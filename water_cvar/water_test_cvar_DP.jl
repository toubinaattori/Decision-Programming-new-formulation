using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using Dates
using DelimitedFiles
using Dates
using Random



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
            optimizer = optimizer_with_attributes(Gurobi.Optimizer,"OptimalityTol" => 1e-09, "MIPGap" => 0)
            set_optimizer(model2, optimizer)
        end  
            optimize!(model2)
    end
    return (dt2, dt1, dt8, solve_time(model2), objective_value(model2),dt5)

end

for No_states in 2:10

    results_DP = zeros(50,6)
    # TODO: Swap these two lines if you want to generate new seeds
    #seed_list = zeros(50)
    seed_list = readdlm("../Results/water_cvar/seeds_$No_states.csv", ',')

    for p in 1:50



 SPW = 1:No_states
 SM = 1:3
 SH = 1:No_states 
 SF = 1:No_states
 SD = 1:No_states
 SC = 1:No_states
 SPoW = 1:No_states


# TODO: uncomment these two lines if you want to generate your own seeds
#seed = Dates.year(now())*Dates.day(now())*Dates.hour(now()) + Dates.minute(now()) + Dates.second(now()) + Dates.millisecond(now())
#seed_list[p] = seed

# TODO:  Use this line if you want to use the same seeds as in the computational experiments
seed = Int(seed_list[i])

PH = rand(Xoshiro(seed-1), No_states)
PH[:] = PH./sum(PH)

PW = zeros(No_states,No_states)
x_V = rand(Xoshiro(seed), No_states,No_states)
for i in 1:No_states
    PW[i,:] = x_V[i,:]./sum(x_V[i,:])
end



F = zeros(3, No_states, No_states)

x_V = zeros(3,No_states,No_states) 
x_V = rand(Xoshiro(seed+1),3,No_states,No_states)
for i in 1:No_states
    F[1,i, :] = x_V[1,i,:]./sum(x_V[1,i,:])
    F[2,i, :] = x_V[2,i,:]./sum(x_V[2,i,:])
    F[3,i, :] = x_V[3,i,:]./sum(x_V[3,i,:])
end

C = zeros(No_states, No_states)

x_V = zeros(No_states,No_states) 
x_V = rand(Xoshiro(seed+2),No_states,No_states)
for i in 1:No_states
    C[i, :] = x_V[i,:]./sum(x_V[i,:])
end


PoW = zeros(No_states, No_states, No_states, No_states)
x_V = rand(Xoshiro(seed+3),No_states,No_states,No_states,No_states)
for i in 1:No_states, j in 1:No_states, u in 1:No_states
    PoW[i, j,u, :] = x_V[i,j,u,:]./sum(x_V[i,j,u,:])
end


 U = zeros(3, No_states, No_states)


 profits2 = zeros(No_states) 
 profits2 = sort(rand(Xoshiro(seed+4), 0:500,No_states))
 costs = zeros(No_states)
 costs = rand(Xoshiro(seed+5),0:100,No_states)
 costs2 = zeros(1) 
 costs2 = rand(Xoshiro(seed+6), 0:100)
 scaler = abs(minimum(profits2) - maximum(costs) - maximum(costs2) - 1)
 for k in 1:3, j in 1:No_states, i in 1:No_states
    U[k,i,j] = profits2[i] - costs[j] - costs2[1] + scaler
 end


println("starting DP")
t = @async DP_model(No_states,PH,PW,F,PoW,U,C)
solutions = [3700,3700,3700,3700,0.0, 3700]
end_time = now() + Dates.Second(3601)
while now() <= end_time
    sleep(0.01)
    if istaskdone(t)
        (dt2, dt1, dt8, dt10, obj, dt5) = fetch(t)
        solutions[1] =dt2
        solutions[2] =dt1
        solutions[3] =dt8
        solutions[4] =dt10
        solutions[5] =obj
        solutions[6] = dt5
        break
    end
end


results_DP[p,1] = solutions[1]
results_DP[p,2] = solutions[2]
results_DP[p,3] = solutions[3]
results_DP[p,4] = solutions[4]
results_DP[p,5] = solutions[5]
results_DP[p,6] = solutions[6]



println("DP parameter calculation time: ",solutions[1])
println("DP model building time: ",solutions[2])
println("DP solve time: ", solutions[4])
println("DP other time: ", solutions[6])
println("DP total: ",  solutions[3])
println("States: ", No_states)
println("Iteration; ", p)

end
writedlm("../Results/water_cvar/results_DP_$No_states.csv", results_DP, ',')
# TODO: Uncomment this if you are creating your own seeds
#writedlm("../Results/water_cvar/seeds_$No_states.csv", seed_list, ',')
end



