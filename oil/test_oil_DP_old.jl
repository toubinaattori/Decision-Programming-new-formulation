using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random



function old_DP_model(diagram::InfluenceDiagram)
    dt9 = @elapsed begin
        dt4 = @elapsed begin
            model3 = Model()
            z = DecisionVariables(model3, diagram, names=false)
            variables = PathCompatibilityVariables(model3, diagram, z, names=false)
            EV = expected_value(model3, diagram, variables)
            @objective(model3, Max, EV)
        end
            
        dt7 = @elapsed begin
            @info("Starting the optimization process.")
            @info("Time: ", now())
            optimizer = optimizer_with_attributes(Gurobi.Optimizer,"OptimalityTol" => 1e-09, "MIPGap" => 0)
            set_optimizer(model3, optimizer)
        end  
            if dt4 + dt7 >= 3600
                return(3700,3700,3700,3700,0.0,"NO_SOLUTION")
            end
            set_time_limit_sec(model3, 3600.0-dt7-dt4)
            optimize!(model3)
    end
    status = primal_status(model3)
    if status == JuMP.NO_SOLUTION
        return (dt9, dt4, dt7, solve_time(model3), 0.0, status)
    else
        if is_solved_and_feasible(model3)
            return (dt9, dt4, dt7, solve_time(model3), objective_value(model3),"OPTIMAL")
        else
            return (dt9, dt4, dt7, solve_time(model3), objective_value(model3),"FEASIBLE")
        end
    end
end


No_states = 3
No = 8
no_rounds = 50

for No_tests in 1:No
    
    results_old_DP = zeros(50,6) 

    seeds = readdlm("../Results/oil/seeds_$No_tests.csv", ',')

    NN = No_tests*2+1


    for p in 1:no_rounds

        SO = 1:No_states
        SS = 1:No_states
        ST = 1:2
        SR = 1:No_states+1
        SD = 1:2

        seed = Int(seeds[p])
        @info("Time: ", now())
        @info("Creating the influence diagram.")
        diagram = InfluenceDiagram()

        add_node!(diagram, ChanceNode("O", [], ["$i" for i in 1:No_states]))
        add_node!(diagram, ChanceNode("S", ["O"], ["$i" for i in 1:No_states]))
        for i in 1:No_tests
            add_node!(diagram, DecisionNode("T$i", [], ["yes", "no"]))
            add_node!(diagram, ValueNode("C$i", ["T$i"]))
            add_node!(diagram, ChanceNode("Tr$i", ["S","T$i"], ["N/A", ["$i" for i in 1:No_states]...]))
        end
        add_node!(diagram, DecisionNode("D", ["Tr$i" for i in 1:No_tests], ["yes", "no"]))
        add_node!(diagram, ValueNode("P", ["O","D"]))

        generate_arcs!(diagram)

        O = zeros(No_states)
        S = zeros(No_states, No_states)
        T = zeros(No_tests, 2, No_states, No_states + 1)
        V = zeros(2,No_tests)
        U = zeros(No_states, 2)

        x_O = rand(Xoshiro(seed),No_states)
        O[:] = x_O./sum(x_O)


        X_O = ProbabilityMatrix(diagram, "O")
        X_O[:] = x_O./sum(x_O)
        add_probabilities!(diagram, "O", X_O)

        for i in 1:No_states
            x_S = rand(Xoshiro(seed+1),No_states)
            S[i, :] = x_S./sum(x_S)
        end

        X_V = ProbabilityMatrix(diagram, "S")
        for i in 1:No_states
            X_V["$i", :] = S[i,:]
        end
        add_probabilities!(diagram, "S", X_V)



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


        for i in 1:No_tests
            X_V = ProbabilityMatrix(diagram, "Tr$i")
            for j in 1:No_states
                X_V["$j","yes", :] = T[i, 1, j, :]
                X_V["$j","no", :] = T[i, 2, j, :]
            end
            add_probabilities!(diagram, "Tr$i", X_V)
            Y_T = UtilityMatrix(diagram, "C$i")
            Y_T[:] = [0,V[2,i]]
            add_utilities!(diagram, "C$i", Y_T)
        end



        Y_T = UtilityMatrix(diagram, "P")
        for j in 1:No_states
            Y_T[j,"yes"] = x_U[j]
            Y_T[j,"no"] = c_D
        end
        add_utilities!(diagram, "P", Y_T)

        generate_diagram!(diagram,positive_path_utility=true)
        





        println("starting Old DP model")
        @info("Time: ", now())
        t = @spawn old_DP_model(diagram)   
        solutions = [3700,3700,3700,3700,0.0]
        status = "unknown"
        end_time = now() + Dates.Second(3700)
        while now() <= end_time
            sleep(0.1)
            if istaskdone(t)
                (dt9, dt4, dt7, dt3, obj3,stat) = fetch(t)
                solutions[1] =dt9
                solutions[2] =dt4
                solutions[3] =dt3
                solutions[4] =obj3
                solutions[5] =dt7
                status = stat
                break
            end
        end

        @info("Time: ", now())
        println("Old DP model building time with jl-package: ",solutions[2])
        println("Old DP model solve time with jl-package: ", solutions[3])
        println("Old DP model total: ", solutions[1])
        println("Old DP model other time: ", solutions[5])
        println("Old DP model objective value: ", solutions[4])
        println("iteration: ", p)
        println("No states: ", No_states)
        println("Old DP model status: ", status)

        results_old_DP[p,1] = solutions[1]
        results_old_DP[p,2] = solutions[2]
        results_old_DP[p,3] = solutions[3]
        results_old_DP[p,4] = solutions[4]
        results_old_DP[p,5] = solutions[5]
        if status == "OPTIMAL"
            results_old_DP[p,6] = 1.0
        elseif status == "FEASIBLE"
            results_old_DP[p,6] = 2.0
        else
            results_old_DP[p,6] = 3.0
        end

end

writedlm("../Results/oil/DP_old_$No_tests.csv", results_old_DP, ',')

end
