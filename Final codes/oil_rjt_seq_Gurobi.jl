using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random
using .Threads


j = parse(Int64, ARGS[1])
seed = readdlm("../Revision_1_results/Data/oil_seed_gurobi_$j.csv",',')[1]
seed = Int(seed)
Random.seed!(seed)

function RJT_model(diagram::InfluenceDiagram)
    dt9 = @elapsed begin
        dt11 = @elapsed begin
            dt4 = @elapsed begin
                model3 = Model()
                z = DecisionVariables(model3, diagram, names=false)
                variables = RJTVariables(model3, diagram, z, names=false)
                EV = expected_value(model3, diagram, variables)
                @objective(model3, Max, EV)
            end
            println("Model building time with jl-package: ", dt4)
            dt7 = @elapsed begin
                if dt4 >= 3600
                    return (3700.0, dt4, 3700.0, 3700.0, 0.0, "TIMEOUT")
                end
                @info("Starting the optimization process.")
                @info("Time: ", now())
                optimizer = optimizer_with_attributes(Gurobi.Optimizer)
                set_optimizer(model3, optimizer)
            end
        end
        if dt11 >= 3600
            return (3700.0, dt4, dt7, 3700.0, 0.0, "TIMEOUT")
        end
        set_time_limit_sec(model3, 3600.0)
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
No = 7
no_rounds = 2

for No_tests in [1,[i for i in 1:No]...]

    NN = No_tests*2+1



        SO = 1:No_states
        SS = 1:No_states
        ST = 1:2
        SR = 1:No_states+1
        SD = 1:2


        # TODO:  Use these lines if you want to create your own seeds
        #seed_list[p] = seed

        # TODO:  Use this line if you want to use the same seeds as in the computational experiments
         #seed = Int(seed_list[p])

    

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

        @info("Time: ", now())
        @info("Creating the influence diagram.")
        diagram = InfluenceDiagram()

        add_node!(diagram, ChanceNode("O", [], ["$i" for i in 1:No_states]))
        add_node!(diagram, ChanceNode("S", ["O"], ["$i" for i in 1:No_states]))
        for i in 1:No_tests
            add_node!(diagram, DecisionNode("T$i", [["Tr$j" for j in 1:i-1]...], ["yes", "no"]))
            add_node!(diagram, ValueNode("C$i", ["T$i"]))
            add_node!(diagram, ChanceNode("Tr$i", ["S","T$i"], ["N/A", ["$i" for i in 1:No_states]...]))
        end
        add_node!(diagram, DecisionNode("D", ["Tr$i" for i in 1:No_tests], ["yes", "no"]))
        add_node!(diagram, ValueNode("P", ["O","D"]))

        generate_arcs!(diagram)


        X_O = ProbabilityMatrix(diagram, "O")
        X_O[:] = x_O./sum(x_O)
        add_probabilities!(diagram, "O", X_O)

        X_V = ProbabilityMatrix(diagram, "S")
        for i in 1:No_states
            X_V["$i", :] = S[i,:]
        end
        add_probabilities!(diagram, "S", X_V)
        

        for i in 1:No_tests
            X_V = ProbabilityMatrix(diagram, "Tr$i")
            for j in 1:No_states
                X_V["$j","yes", :] = T[i, 1, j, :]
                X_V["$j","no", :] = T[i, 2, j, :]
            end
            add_probabilities!(diagram, "Tr$i", X_V)
        end

        for i in 1:No_tests
            Y_V = UtilityMatrix(diagram, "C$i")
            Y_V["yes"] = V[1,i]
            Y_V["no"] = V[2,i]
            add_utilities!(diagram, "C$i", Y_V)
        end


        Y_T = UtilityMatrix(diagram, "P")
        for o in 1:No_states, d in 1:2
            Y_T[o,d] = U[o,d]
        end
        add_utilities!(diagram, "P", Y_T)

        generate_diagram!(diagram)
        





        println("starting RJT")
        @info("Time: ", now())
        solutions = RJT_model(diagram)   
        status = "blah"

        @info("Time: ", now())
        println("RJT model building time with jl-package: ",solutions[2])
        println("RJT solve time with jl-package: ", solutions[4])
        println("RJT total: ", solutions[1])
        println("RJT objective value: ", solutions[5])
        println("No states: ", No_tests)
        println("RJT status: ", status)

        writedlm("../Revision_1_results/Results_Oil_RJT/results_Oil_RJT_triton_gurobi_$(j)_$No_tests.csv", solutions)

end
