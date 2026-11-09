using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random

global j = parse(Int64, ARGS[1])
seed = Dates.value(now()) + j
Random.seed!(seed)
writedlm("../Revision_1_results/Data/oil_seed_non_seq_gurobi_$j.csv",[seed],",")

function DP_model(SS::UnitRange{Int64}, SO::UnitRange{Int64}, SD::UnitRange{Int64},No_states::Int64, No_tests::Int64,  O::Vector{Float64}, S::Matrix{Float64}, T::Array{Float64,4}, V::Matrix{Float64}, U::Matrix{Float64})

        ST = 1:No_states+1
        STI = 1:2
        dt8 = @elapsed begin
            # Calculating DP-parameters
            PU = zeros([2 for i in 1:No_tests]...,[No_states+1 for i in 1:No_tests]..., 2)
            dt2 = @elapsed begin
                for d in SD, s in product([STI for i in 1:No_tests]...,[ST for i in 1:No_tests]...)
                    PU[s...,d] = sum(O[o]*S[o,ss]*prod(T[i,s[i],ss,s[No_tests + i]] for i in 1:No_tests)*(U[o,d] + sum(V[s[i],i] for i in 1:No_tests)) for ss in SS, o in SO)
                end
            end


            # Calculating model build time
            dt1 = @elapsed begin
                a = [[STI for i in 1:No_tests]...,[ST for i in 1:No_tests]..., SD]
                b = [[ST for i in 1:No_tests]..., SD]
                e =  [ST for i in 1:No_tests]
                model2 = Model()
                set_time_limit_sec(model2, 3600.0)
                y = Array{VariableRef}(undef, Tuple(length.(a)))
                for index in CartesianIndices(y)
                    y[index] = @variable(model2, upper_bound = 1.0, lower_bound = 0.0,base_name = "y_$index")
                end

                z = Array{VariableRef}(undef, Tuple(length.(b)))
                for index in CartesianIndices(z)
                    z[index] = @variable(model2,base_name = "z_$index", binary = true)
                end
                @variable(model2, z2[1:No_tests, STI], base_name = "z2_$index",  binary = true)
                @constraint(model2,[t in 1:No_tests], sum(z2[t,:]) == 1)

                for index in CartesianIndices(Tuple(length.(e)))
                    c = [Tuple(index)..., SD]
                    @constraint(model2, sum(z[index2] for index2 in CartesianIndices(Tuple([x isa Integer ? (x:x) : x for x in c]))) == 1)
                end


                for index in CartesianIndices(Tuple(length.(b)))
                    c = [[STI for i in 1:No_tests]...,Tuple(index)[1:end-1]...,Tuple(index)[end]]
                    @constraint(model2, sum(y[index2] for index2 in CartesianIndices(Tuple([x isa Integer ? (x:x) : x for x in c]))) <= length(CartesianIndices(Tuple(length.(c))))*z[index])
                end
                for i in 1:No_tests
                    for k in STI
                        c = [[STI for j in 1:(i-1)]...,k, [STI for j in (i+1):No_tests]...,[ST for j in 1:No_tests]..., SD]
                        @constraint(model2, sum(y[index2] for index2 in CartesianIndices(Tuple(length.(c)))) <= length(CartesianIndices(Tuple(length.(c))))*(1+No_states)*z2[i,k])
                    end
                end


                for index in CartesianIndices(Tuple(length.(e)))
                    c = [[STI for i in 1:No_tests]...,Tuple(index)..., SD]
                    @constraint(model2, sum(y[index2] for index2 in CartesianIndices(Tuple([x isa Integer ? (x:x) : x for x in c]))) <= 1)
                end

                @objective(model2, Max, sum(y[index]*PU[Tuple(index)...] for index in CartesianIndices(Tuple(length.(a)))))
            end
            println("Model building time with jl-package: ", dt1)
            dt5 = @elapsed begin
                @info("Starting the optimization process.")
                optimizer = optimizer_with_attributes(Gurobi.Optimizer)
                set_optimizer(model2, optimizer)
            end  
            optimize!(model2)
        end

        return (dt2, dt1, dt8,dt5, solve_time(model2), objective_value(model2))
    end


No_states = 3

# Specify Maximum M
No = 8

for No_tests in [1,[i for i in 1:No]...]
    
    NN = No_tests*2+1
    


        SO = 1:No_states
        SS = 1:No_states
        ST = 1:2
        SR = 1:No_states+1
        SD = 1:2

    

        O = zeros(No_states)
        S = zeros(No_states, No_states)
        T = zeros(No_tests, 2, No_states, No_states + 1)
        V = zeros(2,No_tests)
        U = zeros(No_states, 2)

        x_O = rand(No_states)
        O[:] = x_O./sum(x_O)



        for i in 1:No_states
            center = i
            r = 1
            weights = exp.(-r .* abs.(collect(1:No_states) .- center))
            noise = 0.5 .+ rand(No_states)
            S[i, :] = (weights .* noise) ./ sum(weights .* noise)
        end




        for i in 1:No_tests, j in 1:No_states
            center = j
            r = 1
            weights = exp.(-r .* abs.(collect(1:No_states) .- center))
            noise = 0.5 .+ rand(No_states)
            pp = (weights .* noise) ./ sum(weights .* noise)
            T[i, 1, j, :] = [0,pp...]
            T[i, 2, j, :] = [1,zeros(No_states)...]
        end

        for i in 1:No_tests
            x_V = rand(1:50)
            V[1,i] = -x_V + 50*No_tests
            V[2,i] = 50*No_tests
        end

        x_U = sort(rand(200:1100,No_states))
        x_U[1] = 100
        c_D = rand(50:90)

        for i in 1:No_states
            U[i,1] = x_U[i] - c_D
            U[i,2] = 100
        end



    println("starting DP")
    #t = @async DP_model(SS, SO, SD,No_states, No_tests, scaler, O, S, T, V, U)
    solutions_DP = DP_model(SS, SO, SD,No_states, No_tests, O, S, T, V, U)
    #solutions_DP = [3700,3700,3700,3700,3700,0.0]
    #end_time = now() + Dates.Second(3601)
    #while now() <= end_time
    #    sleep(0.01)
    #    if istaskdone(t)
    #        (dt2, dt1, dt8,dt5, solv_tim, obj2) = fetch(t)
    #        solutions_DP[1] =dt2
    #        solutions_DP[2] =dt1
    #        solutions_DP[3] =dt8
    #        solutions_DP[4] =solv_tim
    #        solutions_DP[5] =obj2
    #        solutions_DP[6] = dt5
    #        break
    #    end
    #end



    println("DP parameter calculation time: ",solutions_DP[1])
    println("DP model building time: ",solutions_DP[2])
    println("DP solve time: ", solutions_DP[5])
    println("DP objective value: ", solutions_DP[6])
    println("DP total: ",  solutions_DP[3])

    println("Tests: ", No_tests)
    println(solutions_DP[5])

    writedlm("../Revision_1_results/Results_Oil_DP/results_Oil_DP_triton_non_seq_gurobi_$(j)_$No_tests.csv", solutions_DP)

end
