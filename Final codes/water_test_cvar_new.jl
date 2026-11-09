using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using Dates
using DelimitedFiles
using Dates
using Random

global j = parse(Int64, ARGS[1])
seed = Dates.value(now()) + j
Random.seed!(seed)
writedlm("../Revision_1_results/Data/water_cvar_seed_$j.csv",[seed],",")

function DP_model(No_states::Int64,PH::Vector{Float64},PW::Array{Float64,2},F::Array{Float64,3},PoW::Array{Float64,4},U::Array{Float64,3},C::Array{Float64,2})

    SPW = 1:No_states
    SM = 1:3
    SH = 1:No_states 
    SF = 1:No_states
    SD = 1:No_states
    SC = 1:No_states
    SPoW = 1:No_states
    dt8 = @elapsed begin
        u = collect(Iterators.flatten(U[m,c,pow] for m in SM, c in SC, pow in SPoW))
        u_sorted = sort(u)
        u_min = u_sorted[1]
        u_max = u_sorted[end]
        M = u_max - u_min
        u_diff = diff(u_sorted)
        if isempty(filter(!iszero, u_diff))
            return u_min    # All utilities are the same, CVaR is equal to that constant utility value
        else
            ϵ = minimum(filter(!iszero, abs.(u_diff))) / 2 
        end
        u_unique = unique(u_sorted)
        dt2 = @elapsed begin
            P = zeros(length(u),3, No_states, No_states)
            for m in SM, f in SF, d in SD, ut in 1:length(u_unique)
                for c in SC, p in SPW, pow in SPoW, h in SH
                    if U[m,c,pow] == u_unique[ut]
                        P[ut,m,f,d] += PH[h]*F[m,p,f]*PW[h,p]*C[d,c]*PoW[d,p,h,pow]
                    end
                end
            end
        end

        dt1 = @elapsed begin
            model2 = Model()
            @variable(model2, y[SM, SF, SD], upper_bound = 1.0, lower_bound = 0.0)
            @variable(model2, z1[SM], Bin)
            @variable(model2, z2[SF,SD], Bin)
            
            @constraint(model2, sum(z1) == 1)
            @constraint(model2, [f in SF], sum(z2[f,d] for d in SD) == 1)
            @constraint(model2, [m in SM], sum(y[m,f,d] for f in SF, d in SD ) <= 3*2*No_states^3*z1[m])
            @constraint(model2, [f in SF, d in SD], sum(y[m,f,d] for m in SM) <= 3*2*No_states^3*z2[f,d])
            @constraint(model2, [f in SF], sum(y[m,f,d] for m in SM, d in SD) == 1)

            α = 0.2
        
            # Variables and constraints
            η = @variable(model2)
            @constraint(model2, η ≥ u_min)
            @constraint(model2, η ≤ u_max)
            ρ′_s = Dict{Float64, VariableRef}()
            ρ_s = Dict{Float64, VariableRef}()
            λ′_s = Dict{Float64, VariableRef}()
            λ_s = Dict{Float64, VariableRef}()
            for (j,u) in enumerate(u_unique)
                λ = @variable(model2, binary=true)
                λ′ = @variable(model2, binary=true)
                ρ = @variable(model2)
                ρ′ = @variable(model2)
                @constraint(model2, η - u ≤ M * λ)
                @constraint(model2, η - u ≥ (M + ϵ) * λ - M)
                @constraint(model2, η - u ≤ (M + ϵ) * λ′ - ϵ)
                @constraint(model2, η - u ≥ M * (λ′ - 1))
                @constraint(model2, 0 ≤ ρ)
                @constraint(model2, 0 ≤ ρ′)
                @constraint(model2, ρ ≤ λ)
                @constraint(model2, ρ′ ≤ λ′)
                @constraint(model2, ρ ≤ ρ′)
                @constraint(model2, ρ′ ≤ sum(y[m,f,d] * P[j,m,f,d] for f in SF, d in SD, m in SM))
                @constraint(model2, sum(y[m,f,d] * P[j,m,f,d] for f in SF, d in SD,  m in SM) - (1 - λ) ≤ ρ)
                ρ′_s[u] = ρ′
                ρ_s[u] = ρ
                λ′_s[u] = λ′
                λ_s[u] = λ
            end
            @constraint(model2, sum(values(ρ′_s)) == α)
        
            # Return CVaR as an expression
            CVaRN = @expression(model2, (sum(ρ_bar * u for (u, ρ_bar) in ρ′_s)/α))

            @objective(model2, Max, CVaRN)
            
        end
        dt5 = @elapsed begin
            @info("Starting the optimization process.")
            optimizer = optimizer_with_attributes(Gurobi.Optimizer, "Threads" => 1, "Presolve" => 0, "Heuristics" => 0, "NodeMethod" => 1, "Cuts" => 1, "MIPFocus" => 2)
            set_optimizer(model2, optimizer)
        end  
            println("Model building time with $No_states states: ", dt2 + dt1 + dt5)
            optimize!(model2)
    end
    return (dt2, dt1, dt8, solve_time(model2), objective_value(model2),dt5)

end

for No_states in [2,collect(2:10)...]
    results_DP = zeros(50,6)
    # TODO: Swap these two lines if you want to generate new seeds
    seed_list = zeros(50)
    #seed_list = readdlm("../Results/water_cvar/seeds_$No_states.csv", ',')



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
    noise = 0.5 .+ rand( No_states)
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


println("starting DP")
solutions = DP_model(No_states,PH,PW,F,PoW,U,C)





println("DP parameter calculation time: ",solutions[1])
println("DP model building time: ",solutions[2])
println("DP solve time: ", solutions[4])
println("DP other time: ", solutions[6])
println("DP total: ",  solutions[3])
println("DP objective value: ",  solutions[5])
println("States: ", No_states)

   writedlm("../Revision_1_results/Results_Water_CVaR_DP/results_water_CVaR_DP_triton_$(j)_$No_states.csv", solutions)


end



