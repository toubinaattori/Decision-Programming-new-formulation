using Pkg
using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random
using .Threads

global j = parse(Int64, ARGS[1])
readseed = readdlm("../Revision_1_results/Data/water_cvar_seed_gurobi_$j.csv",',')[1]
seed = Int(readseed)
Random.seed!(seed)

function old_DP_model(diagram::InfluenceDiagram,No_states::Int64)
        dt9 = @elapsed begin
        dt4 = @elapsed begin
            model = Model()
            z = DecisionVariables(model, diagram, names=false)
            variables = PathCompatibilityVariables(model, diagram, z, names=false)
            u = collect(Iterators.flatten(diagram.U(s, diagram.translation) for s in keys(variables)))
            u_sorted = sort(u)
            u_min = u_sorted[1]
            u_max = u_sorted[end]
            M = u_max - u_min + 1
            u_diff = diff(u_sorted)
            if isempty(filter(!iszero, u_diff))
                return u_min    # All utilities are the same, CVaR is equal to that constant utility value
            else
                ϵ = minimum(filter(!iszero, abs.(u_diff))) / 2 
            end

            # Variables and constraints
            η = @variable(model)
            @constraint(model, η ≥ u_min)
            @constraint(model, η ≤ u_max)
            ρ′_s = Dict{Path, VariableRef}()
            for (s, x) in variables
                u_s = diagram.U(s, diagram.translation)
                λ = @variable(model, binary=true)
                λ′ = @variable(model, binary=true)
                ρ = @variable(model)
                ρ′ = @variable(model)
                @constraint(model, η - u_s ≤ M * λ)
                @constraint(model, η - u_s ≥ (M + ϵ) * λ  - M)
                @constraint(model, η - u_s ≤ (M + ϵ) * λ′ - ϵ)
                @constraint(model, η - u_s ≥ M * (λ′ - 1))
                 @constraint(model, 0 ≤ ρ)
                 @constraint(model, 0 ≤ ρ′)
                 @constraint(model, ρ ≤ λ )
                 @constraint(model, ρ′ ≤ λ′)
                 @constraint(model, ρ ≤ ρ′)
                 @constraint(model, ρ′ ≤ x * diagram.P(s) )
                 @constraint(model, (x * diagram.P(s) - (1 - λ)) ≤ ρ)
                ρ′_s[s] = ρ′
            end
            @constraint(model, sum(values(ρ′_s)) == 0.2 )

            # Return CVaR as an expression
            CVaR = @expression(model, sum(ρ_bar * diagram.U(s, diagram.translation) for (s, ρ_bar) in ρ′_s) / (0.2))

            @objective(model, Max, CVaR)
        end
            
        dt7 = @elapsed begin
            @info("Starting the optimization process.")
            optimizer = optimizer_with_attributes(Gurobi.Optimizer)
            set_optimizer(model, optimizer)
        end
        set_time_limit_sec(model, maximum([3600.0-dt7-dt4,0.0001])) 
        println("Model building time with jl-package with $No_states states: ", dt4 + dt7) 
        optimize!(model)
        println(is_solved_and_feasible(model))
    end
    status = primal_status(model)
    if status == JuMP.NO_SOLUTION
        return (dt9, dt4, dt7, 3700, 0.0, status)
    else
        if is_solved_and_feasible(model)
            return (dt9, dt4, dt7, solve_time(model), objective_value(model),"OPTIMAL")
        else
            return (dt9, dt4, dt7, solve_time(model), objective_value(model),"FEASIBLE")
        end
    end
end

for No_states in [2,collect(2:10)...]

    results_old_DP = zeros(50,6)




        SW = 1:No_states
        SH = 1:No_states 
        SPW = 1:No_states
        SM = 1:3
        SF = 1:No_states
        SD = 1:No_states
        SC = 1:No_states
        SPoW = 1:No_states

        diagram = InfluenceDiagram()

        add_node!(diagram, ChanceNode("H", [], ["$i" for i in 1:No_states]))
        add_node!(diagram, ChanceNode("PW", ["H"], ["$i" for i in 1:No_states]))
        add_node!(diagram, DecisionNode("M", [], ["yes", "no","maybe"]))
        add_node!(diagram, ChanceNode("F", ["M","PW"], ["$i" for i in 1:No_states]))
        add_node!(diagram, DecisionNode("D", ["F"], ["$i" for i in 1:No_states]))
        add_node!(diagram, ChanceNode("C", ["D"], ["$i" for i in 1:No_states]))
        add_node!(diagram, ChanceNode("PoW", ["D","PW","H"], ["$i" for i in 1:No_states]))
        add_node!(diagram, ValueNode("U", ["M","C","PoW"]))
        
        generate_arcs!(diagram)
        
        
 SPW = 1:No_states
 SM = 1:3
 SH = 1:No_states 
 SF = 1:No_states
 SD = 1:No_states
 SC = 1:No_states
 SPoW = 1:No_states



PH = rand( No_states)
PH[:] = PH./sum(PH)

PW = zeros(No_states,No_states)
for h in 1:No_states
    # Small H → center near No_states (good quality); large H → center near 1 (poor quality)
    center = No_states + 1 - h
    weights = exp.(-0.4 .* abs.(collect(1:No_states) .- center))
    noise = 0.5 .+ rand( No_states)
    PW[h, :] = (weights .* noise) ./ sum(weights .* noise)
end



F = zeros(3, No_states, No_states)

# m=1 (few samples) → low sharpness (flat); m=3 (many samples) → high sharpness (tight)
for m in 1:3, p in 1:No_states
    sharpness = 0.5 * (m - 1) + 0.3   # m=1→0.3, m=2→0.8, m=3→1.3
    weights = exp.(-sharpness .* abs.(collect(1:No_states) .- p))
    noise = 0.5 .+ rand( No_states)
    F[m, p, :] = (weights .* noise) ./ sum(weights .* noise)
end

C = zeros(No_states, No_states)

for d in 1:No_states
    weights = exp.(-0.8 .* abs.(collect(1:No_states) .- d))
    noise = 0.5 .+ rand( No_states)
    C[d, :] = (weights .* noise) ./ sum(weights .* noise)
end


PoW = zeros(No_states, No_states, No_states, No_states)

# Center = prior quality p, shifted up by dilution d and down by exogenous e (both scaled to [0, No_states-1])
scale = (No_states - 1) / max(1, No_states - 1)
for d in 1:No_states, p in 1:No_states, e in 1:No_states
    dilution_shift   =  (d - 1) * scale   # d=1 → no improvement; d=No_states → max improvement
    exogenous_shift  = -(e - 1) * scale   # e=1 → no degradation;  e=No_states → max degradation
    center = clamp(p + dilution_shift + exogenous_shift, 1, No_states)
    weights = exp.(-0.8 .* abs.(collect(1:No_states) .- center))
    noise = 0.5 .+ rand(No_states)
    PoW[d, p, e, :] = (weights .* noise) ./ sum(weights .* noise)
end


 U = zeros(3, No_states, No_states)


 profits2 = zeros(No_states) 
        profits2 = sort(500 .* rand(No_states))
        costs = zeros(No_states)
        costs = 100 .* rand(No_states)
        costs2 = zeros(1) 
        costs2 = 100 .* rand(No_states)
 scaler = abs(minimum(profits2) - maximum(costs) - maximum(costs2) - 1)
 for k in 1:3, j in 1:No_states, i in 1:No_states
    U[k,i,j] = profits2[i] - costs[j] - costs2[1] + scaler
 end


        
        X_H = ProbabilityMatrix(diagram, "H")
        X_H[:] = PH
        add_probabilities!(diagram, "H", X_H)
        
        X_PW = ProbabilityMatrix(diagram, "PW")
        for i in 1:No_states
            X_PW["$i", :] = PW[i,:]
        end
        add_probabilities!(diagram, "PW", X_PW)
        
        
        X_F = ProbabilityMatrix(diagram, "F")
        for i in 1:No_states
            X_F["yes","$i", :] = F[1,i, :]
            X_F["no","$i", :] = F[2,i, :]
            X_F["maybe","$i", :] = F[3,i, :]
        end
        add_probabilities!(diagram, "F", X_F)
        
        X_C = ProbabilityMatrix(diagram, "C")
        for i in 1:No_states
            X_C["$i", :] = C[i, :]
        end
        add_probabilities!(diagram, "C", X_C)
        
        X_PoW = ProbabilityMatrix(diagram, "PoW")
        for i in 1:No_states, j in 1:No_states, u in 1:No_states
            X_PoW["$i", "$j","$u", :] = PoW[i, j,u, :]
        end
        add_probabilities!(diagram, "PoW", X_PoW)
        
        
        Y_T = UtilityMatrix(diagram, "U")
        for k in ["yes","no","maybe"], j in 1:No_states, i in 1:No_states
            Y_T[k,"$i","$j"] = U[1,i,j]
        end
        add_utilities!(diagram, "U", Y_T)
        
        
        generate_diagram!(diagram)


        println("starting old_DP")
        solutions = old_DP_model(diagram,No_states)


        println("old_DP model building time with jl-package: ",solutions[2])
        println("old_DP solve time with jl-package: ", solutions[4])
        println("old_DP total: ", solutions[1])
        println("old_DP other time: ", solutions[3])
        println("old_DP objective value: ", solutions[5])
        println("No states: ", No_states)
        writedlm("../Revision_1_results/Results_Water_CVaR_DP_old/results_Water_CVaR_DP_old_triton_gurobi_$(j)_$No_states.csv", solutions[1:5])


end

