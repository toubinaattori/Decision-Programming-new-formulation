using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random

global j = parse(Int64, ARGS[1])
seed = Dates.value(now()) + j
Random.seed!(seed)
writedlm("../Revision_1_results/Data/pig_farm_gurobi_seed_$j.csv",[seed],",")

function DP_model(No_periods::Int64, scaler::Float64, PH0::Vector{Float64}, PH::Array{Float64}, PT1::Matrix{Float64}, PT2::Matrix{Float64}, PT3::Matrix{Float64}, U::Vector{Float64}, C::Vector{Float64})

        dt8 = @elapsed begin
            # Calculating DP-parameters
            p = [[2,2,2,2] for i in 1:No_periods]
            PU = zeros(collect(Iterators.flatten(p))...)
            dt2 = @elapsed begin
                for t1 in product([1:2 for i in 1:No_periods]...), t2 in product([1:2 for i in 1:No_periods]...), t3 in product([1:2 for i in 1:No_periods]...), d in product([1:2 for i in 1:No_periods]...)
                    s = collect(Iterators.flatten(zip(t1, t2, t3, d)))
                    PU[s...] = sum(PH0[h[1]]*prod(PT1[h[i],s[4*i-3]]*PT2[h[i],s[4*i-2]]*PT3[h[i],s[4*i-1]]*PH[h[i],s[4*i],h[i+1]] for i in 1:No_periods)*(sum(C[s[4*j]] for j in 1:No_periods) + U[h[end]] + scaler) for h in product([1:2 for i in 1:No_periods+1]...))
                end
            end
            # Calculating model build time
            dt1 = @elapsed begin
                model = Model()
                set_time_limit_sec(model, 3600.0)
            @variable(model, z[1:No_periods, 1:2, 1:2, 1:2, 1:2], base_name = "z",  binary = true)
            @constraint(model, [t in 1:No_periods, t2 in 1:2, t1 in 1:2, t3 in 1:2], sum(z[t,t2,t1,t3,:]) == 1)

            y = Array{VariableRef}(undef, ntuple(_ -> 2, 4*No_periods))
            for index in CartesianIndices(y)
                y[index] = @variable(model, upper_bound = 1.0, lower_bound = 0.0,base_name = "y_$index")
            end

                a = vcat([[1:2,1:2,1:2,1:2] for _ in 1:No_periods-1]...)
                b = [1:2,1:2,1:2,1:2]


                for index in CartesianIndices(Tuple(length.(b))), period in 1:No_periods
                    ranges = vcat(
                        [[1:2 for _ in 1:4] for _ in 1:period-1]...,
                        [index[k]:index[k] for k in 1:4],
                        [[1:2 for _ in 1:4] for _ in period+1:No_periods]...
                    )
                    @constraint(model, sum(y[i...] for i in Iterators.product(ranges...)) <= length(CartesianIndices(Tuple(length.(a))))*z[period,Tuple(index)...])
                end


                for index in Iterators.product([1:2 for _ in 1:3*No_periods]...)
                    ranges = vcat([[index[3*(t-1)+1]:index[3*(t-1)+1], index[3*(t-1)+2]:index[3*(t-1)+2], index[3*(t-1)+3]:index[3*(t-1)+3], 1:2] for t in 1:No_periods]...)
                    @constraint(model, sum(y[i...] for i in Iterators.product(ranges...)) <= 1)
                end

                @objective(model, Max, sum(y[index]*PU[Tuple(index)...] for index in CartesianIndices(ntuple(_ -> 2, 4*No_periods))))
            end
            println("DP model building time: ", dt1+dt2)
            dt5 = @elapsed begin
                @info("Starting the optimization process.")
                optimizer = optimizer_with_attributes(Gurobi.Optimizer, "Threads" => 1, "Presolve" => 0, "Heuristics" => 0, "NodeMethod" => 1, "Cuts" => 1, "MIPFocus" => 2)
                set_optimizer(model, optimizer)
            end  
            set_time_limit_sec(model, 3610.0)
            optimize!(model)
        end
        if is_solved_and_feasible(model3)
            return (dt2, dt1, dt8, dt5, solve_time(model), objective_value(model))
        else
            return (dt2, dt1, dt8,dt5, 3700.0, 0.0)
        end
    end


for No_periods in 1:7
    

            p_ill = rand()*0.2                  # Initial probability of pig being ill

        spec1 = 0.5 + rand()*0.5             # Specificity of the test1
        sens1 = 0.5 + rand()*0.5             # Sensitivity of the test1
        spec2 = 0.5 + rand()*0.5             # Specificity of the test2
        sens2 = 0.5 + rand()*0.5             # Sensitivity of the test2
        spec3 = 0.5 + rand()*0.5             # Specificity of the test3
        sens4 = 0.5 + rand()*0.5             # Sensitivity of the test3

        p_vals = sort(rand(4))              # Health state update probabilities
        p_hpi = p_vals[1]                   # Probability of ill untreated pig to become healthy, assumed to be the smallest of the random probabilities
        p_hti = p_vals[2]                   # Probability of ill treated pig to become healthy
        p_hph = p_vals[3]                   # Probability of healthy untreated pig to remain healthy
        p_hth = p_vals[4]                   # Probability of healthy treated pig to remain healthy, assumed to be the largest of the random probabilities

        c_treat = -rand()*100               # Cost of treatment
        u_sell = sort(rand(2).*1000).+100   # Utility of selling a pig, ill or healthy (healthy is more valuable)

        PH0 = [p_ill, 1-p_ill]

        PT1 = [1-sens1  sens1; spec1   1-spec1]
        PT2 = [1-sens2  sens2; spec2   1-spec2]
        PT3 = [1-sens4  sens4; spec3   1-spec3]

        PH = zeros(2, 2, 2)
        PH[1, 1, :] = [1-p_hti, p_hti]
        PH[1, 2, :] = [1-p_hpi, p_hpi]
        PH[2, 1, :] = [1-p_hth, p_hth]
        PH[2, 2, :] = [1-p_hph, p_hph]

        U = [u_sell[1], u_sell[2]]
        C = [c_treat, 0.0]


        
        # required to make utilities strictly positive
        scaler = abs(minimum([c_treat*No_periods + u_sell[2] - 1,0]))
        
        
        println("starting DP")
        solutions_DP = DP_model(No_periods, scaler, PH0, PH, PT1, PT2, PT3, U, C)


    println("DP parameter calculation time: ",solutions_DP[1])
    println("DP model building time: ",solutions_DP[2])
    println("DP solve time: ", solutions_DP[5])
    println("DP objective value: ", solutions_DP[6] - scaler)
    println("DP total: ",  solutions_DP[3])

    println("Tests: ", No_periods)
    println(solutions_DP[5])

    writedlm("../Revision_1_results/Results_Pig_Farm_DP_gurobi/results_Pig_Farm_DP_triton_gurobi_$(j)_$No_periods.csv", solutions_DP)



end
