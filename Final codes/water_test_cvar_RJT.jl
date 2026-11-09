using Pkg
using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random
using .Threads


global j = parse(Int64, ARGS[1])
readseed = readdlm("../Revision_1_results/Data/water_cvar_seed_$j.csv",',')[1]
seed = Int(readseed)
Random.seed!(seed)

function RJT_model(diagram::InfluenceDiagram,No_states::Int64)
    dt9 = @elapsed begin
        dt4 = @elapsed begin
            α = 0.2
            model = Model()
            z = DecisionVariables(model, diagram, names=false)
            variables = RJTVariables(model, diagram, z, names=false)
            M = maximum(diagram.U.Y[1]) - minimum(diagram.U.Y[1]) + 1
            ϵ = minimum(diff(unique(diagram.U.Y[1]))) / 2

            η = @variable(model)
            ρ′_s = Dict{Float64, VariableRef}()

            value_node_name = first(n for (n, node) in diagram.Nodes if isa(node, ValueNode))
            #Assuming only one preceding node
            preceding_node_name = first(filter(tuple -> tuple[2]==collect(keys(diagram.V))[1], diagram.RJT.arcs))[1]

            #Finding the name and index of differing element between value nodes' information set and its preceding nodes rjt cluster. 
            #This is needed in conditional sums for constraints.
            missing_element = setdiff(diagram.RJT.clusters[preceding_node_name], diagram.Nodes[value_node_name].I_j)
            index_to_remove = findall(x -> x in missing_element, diagram.RJT.clusters[preceding_node_name])

            statevars = variables.data[preceding_node_name].statevars
            statevars_dims = collect(size(statevars))
            statevars_dims_ranges = [1:d for d in statevars_dims]

            function remove_index(old_tuple::NTuple{N, Int64}, index::Vector{Int64}) where N
                vector = [old_tuple[i] for i in 1:length(old_tuple) if !(i in index)]
                return collect(ntuple(i -> vector[i], N-length(index)))
            end

            for u in unique(diagram.U.Y[1])
                λ = @variable(model, binary=true)
                λ′ = @variable(model, binary=true)
                ρ = @variable(model)
                ρ′ = @variable(model)
                @constraint(model, η - u ≤ M * λ)
                @constraint(model, η - u ≥ (M + ϵ) * λ - M)
                @constraint(model, η - u ≤ (M + ϵ) * λ′ - ϵ)
                @constraint(model, η - u ≥ M * (λ′ - 1))
                @constraint(model, 0 ≤ ρ)
                @constraint(model, 0 ≤ ρ′)
                @constraint(model, ρ ≤ λ)
                @constraint(model, ρ′ ≤ λ′)
                @constraint(model, ρ ≤ ρ′)

                p = @expression(model, sum(statevars[indices...] for indices in product(statevars_dims_ranges...) if diagram.U.Y[1][remove_index(indices, index_to_remove)...] == u))
                @constraint(model, ρ′ ≤ p)
                @constraint(model, (p - (1 - λ)) ≤ ρ)

                ρ′_s[u] = ρ′
            end

            @constraint(model, sum(values(ρ′_s)) == α)
            CVaR = @expression(model, (sum(ρ_bar * u for (u, ρ_bar) in ρ′_s)/α))
            @objective(model, Max, CVaR)
        end
            
        dt7 = @elapsed begin
            @info("Starting the optimization process.")
            optimizer = optimizer_with_attributes(Gurobi.Optimizer, "Threads" => 1, "Presolve" => 0, "Heuristics" => 0, "NodeMethod" => 1, "Cuts" => 1, "MIPFocus" => 2)
            set_optimizer(model, optimizer)
        end
        println("Model building time with jl-package with $No_states states: ", dt4 + dt7)
        set_time_limit_sec(model, maximum([3600.0-dt7-dt4,0.0001]))  
        optimize!(model)
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


    results_RJT = zeros(50,6)





        SW = 1:No_states
        SH = 1:No_states 
        SPW = 1:No_states
        SM = 1:3
        SF = 1:No_states
        SD = 1:No_states
        SC = 1:No_states
        SPoW = 1:No_states

        seed = 30

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
        
        

        SW = 1:No_states
        SH = 1:No_states 
        SPW = 1:No_states
        SM = 1:3
        SF = 1:No_states
        SD = 1:No_states
        SC = 1:No_states
        SPoW = 1:No_states

        seed = 30

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



        PH = rand(No_states)
        PH[:] = PH./sum(PH)

        PW = zeros(No_states,No_states)
        for h in 1:No_states
            # Small H → center near No_states (good quality); large H → center near 1 (poor quality)
            center = No_states + 1 - h
            weights = exp.(-0.4 .* abs.(collect(1:No_states) .- center))
            noise = 0.5 .+ rand(No_states)
            PW[h, :] = (weights .* noise) ./ sum(weights .* noise)
        end



        F = zeros(3, No_states, No_states)

        # m=1 (few samples) → low sharpness (flat); m=3 (many samples) → high sharpness (tight)
        for m in 1:3, p in 1:No_states
            sharpness = 0.5 * (m - 1) + 0.3   # m=1→0.3, m=2→0.8, m=3→1.3
            weights = exp.(-sharpness .* abs.(collect(1:No_states) .- p))
            noise = 0.5 .+ rand(No_states)
            F[m, p, :] = (weights .* noise) ./ sum(weights .* noise)
        end

        C = zeros(No_states, No_states)

        for d in 1:No_states
            weights = exp.(-0.8 .* abs.(collect(1:No_states) .- d))
            noise = 0.5 .+ rand(No_states)
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


        println("starting RJT")
        solutions = RJT_model(diagram,No_states)


        println("RJT model building time with jl-package: ",solutions[2])
        println("RJT solve time with jl-package: ", solutions[4])
        println("RJT total: ", solutions[1])
        println("RJT other time: ", solutions[2])
        println("RJT objective value: ", solutions[5])
        println("No states: ", No_states)


        writedlm("../Revision_1_results/Results_Water_CVaR_RJT/results_Water_CVaR_RJT_triton_$(j)_$No_states.csv", solutions[1:5])
end


