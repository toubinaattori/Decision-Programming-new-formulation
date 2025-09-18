using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random



function DP_model(SS::UnitRange{Int64}, SO::UnitRange{Int64}, SD::UnitRange{Int64},No_states::Int64, No_tests::Int64, scaler::Float64, O::Vector{Float64}, S::Matrix{Float64}, T::Array{Float64,4}, V::Matrix{Float64}, U::Matrix{Float64})

        ST = 1:No_states
        STI = 1:2
        dt8 = @elapsed begin
            # Calculating DP-parameters
            PU = zeros([2 for i in 1:No_tests]...,[No_states for i in 1:No_tests]..., 2)
            dt2 = @elapsed begin
                for d in SD, s in product([STI for i in 1:No_tests]...,[ST for i in 1:No_tests]...)
                    PU[s...,d] = sum(O[o]*S[o,ss]*prod(T[i,s[i],ss,s[No_tests + i]] for i in 1:No_tests)*(U[o,d] + sum(V[s[i],i] for i in 1:No_tests) + scaler) for ss in SS, o in SO)
                end
            end


            # Calculating model build time
            dt1 = @elapsed begin
                a = [[STI for i in 1:No_tests]...,[ST for i in 1:No_tests]..., SD]
                b = [[ST for i in 1:No_tests]..., SD]
                e =  [ST for i in 1:No_tests]
                f = STI
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
                    @constraint(model2, sum(z[index2] for index2 in CartesianIndices(Tuple(length.(c)))) == 1)
                end


                for index in CartesianIndices(Tuple(length.(b)))
                    c = [Tuple(index)[1:end-1]...,[ST for i in 1:No_tests]...,Tuple(index)[end]]
                    @constraint(model2, sum(y[index2] for index2 in CartesianIndices(Tuple(length.(c)))) <= length(CartesianIndices(Tuple(length.(c))))*(1+No_states)*z[index])
                end

                for i in 1:No_tests
                    for k in STI
                        c = [[STI for j in 1:(i-1)]...,k, [STI for j in (i+1):No_tests]...,[ST for j in 1:No_tests]..., SD]
                        @constraint(model2, sum(y[index2] for index2 in CartesianIndices(Tuple(length.(c)))) <= length(CartesianIndices(Tuple(length.(c))))*(1+No_states)*z2[i,k])
                    end
                end

                for index in CartesianIndices(Tuple(length.(e)))
                    c = [[STI for i in 1:No_tests]...,Tuple(index)..., SD]
                    @constraint(model2, sum(y[index2] for index2 in CartesianIndices(Tuple(length.(c)))) <= 1)
                end

                @objective(model2, Max, sum(y[index]*PU[Tuple(index)...] for index in CartesianIndices(Tuple(length.(a)))))
            end

            dt5 = @elapsed begin
                @info("Starting the optimization process.")
                optimizer = optimizer_with_attributes(Gurobi.Optimizer,"OptimalityTol" => 1e-09, "MIPGap" => 0)
                set_optimizer(model2, optimizer)
            end  
            
            optimize!(model2)
        end

        return (dt2, dt1, dt8,dt5, solve_time(model2), objective_value(model2))
    end


No_states = 3

# Specify Maximum M
No = 3
no_rounds = 50

for No_tests in 1:No
    
    NN = No_tests*2+1
    
    results_DP = zeros(50,6)
    
    # TODO:  Use this if you want to create your own seeds
    #seed_list = zeros(no_rounds)
    # TODO:  Use this if you want to use the same seeds as in the computational experiments
    seed_list = readdlm("../Results/oil/seeds_$No_tests.csv", ',')
    for p in 1:no_rounds

        SO = 1:No_states
        SS = 1:No_states
        ST = 1:2
        SR = 1:No_states+1
        SD = 1:2


        # TODO:  Use these lines if you want to create your own seeds
        #seed = Dates.year(now())*Dates.day(now())*Dates.hour(now()) + Dates.minute(now()) + Dates.second(now()) + Dates.millisecond(now())
        #seed_list[p] = seed

        # TODO:  Use this line if you want to use the same seeds as in the computational experiments
         seed = Int(seed_list[p])

    

        O = zeros(No_states)
        S = zeros(No_states, No_states)
        T = zeros(No_tests, 2, No_states, No_states + 1)
        V = zeros(2,No_tests)
        U = zeros(No_states, 2)

        x_O = rand(Xoshiro(seed),No_states)
        O[:] = x_O./sum(x_O)



        for i in 1:No_states
            x_S = rand(Xoshiro(seed+1),No_states)
            S[i, :] = x_S./sum(x_S)
        end




        for i in 1:No_tests, j in 1:No_states
            x_T = [0,rand(Xoshiro(seed+2),No_states)...]
            T[i, 1, j, :] = x_T./sum(x_T)
            T[i, 2, j, :] = [1,zeros(No_states)...]
        end

        for i in 1:No_tests
            x_V = rand(Xoshiro(seed+3),1:50)
            V[1,i] = -x_V
            V[2,i] = 0
        end

        x_U = sort(rand(Xoshiro(seed+4),100:1000,No_states))
        x_U[1] = 0
        c_D = rand(Xoshiro(seed+5),50:90)

        for i in 1:No_states
            U[i,1] = x_U[i] - c_D
            U[i,2] = 0
        end

        # required to make utilities strictly positive
        scaler = abs(-c_D + sum(V[1,i] for i in 1:No_tests) - 1)



    println("starting DP")
    t = @async DP_model(SS, SO, SD,No_states, No_tests, scaler, O, S, T, V, U)
    solutions_DP = [3700,3700,3700,3700,3700,0.0]
    end_time = now() + Dates.Second(3601)
    while now() <= end_time
        sleep(0.01)
        if istaskdone(t)
            (dt2, dt1, dt8,dt5, solv_tim, obj2) = fetch(t)
            solutions_DP[1] =dt2
            solutions_DP[2] =dt1
            solutions_DP[3] =dt8
            solutions_DP[4] =solv_tim
            solutions_DP[5] =obj2
            solutions_DP[6] = dt5
            break
        end
    end



    println("DP parameter calculation time: ",solutions_DP[1])
    println("DP model building time: ",solutions_DP[2])
    println("DP solve time: ", solutions_DP[4])
    println("DP other time: ", solutions_DP[6])
    println("DP total: ",  solutions_DP[3])

    println("Iteration; ", p)
    println("Tests: ", No_tests)
    println(solutions_DP[5])

    results_DP[p,1] = solutions_DP[1]
    results_DP[p,2] = solutions_DP[2]
    results_DP[p,3] = solutions_DP[3]
    results_DP[p,4] = solutions_DP[4]
    results_DP[p,5] = solutions_DP[5]
    results_DP[p,6] = solutions_DP[6]

end

writedlm("../Results/oil/results_DP_$No_tests.csv", results_DP, ',')
# TODO: # Uncomment this if you are creating your own seeds
#writedlm("../Results/oil/seeds_$No_tests.csv", seed_list, ',')

end
