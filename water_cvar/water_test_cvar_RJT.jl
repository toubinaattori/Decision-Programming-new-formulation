using Pkg
using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random
using .Threads



function RJT_model(diagram::InfluenceDiagram)
    dt9 = @elapsed begin
        dt4 = @elapsed begin
            model3 = Model()
            z = DecisionVariables(model3, diagram, names=false)
            variables = RJTVariables(model3, diagram, z, names=false)
            CVaR = conditional_value_at_risk(model3, diagram, variables,0.2)
            @objective(model3, Max, CVaR)
        end
            
        dt7 = @elapsed begin
            @info("Starting the optimization process.")
            optimizer = optimizer_with_attributes(Gurobi.Optimizer,"OptimalityTol" => 1e-09, "MIPGap" => 0)
            set_optimizer(model3, optimizer)
        end
        set_time_limit_sec(model3, maximum([3600.0-dt7-dt4,0.0001]))  
        optimize!(model3)
    end
    status = primal_status(model3)
    if status == JuMP.NO_SOLUTION
        return (dt9, dt4, dt7, 3700, 0.0, status)
    else
        if is_solved_and_feasible(model3)
            return (dt9, dt4, dt7, solve_time(model3), objective_value(model3),"OPTIMAL")
        else
            return (dt9, dt4, dt7, solve_time(model3), objective_value(model3),"FEASIBLE")
        end
    end
end

for No_states in 2:10

    results_RJT = zeros(50,6)
    seed_list = readdlm("../Results/water_cvar/seeds_$No_states.csv", ',')

    for p in 1:50




        SW = 1:No_states
        SH = 1:No_states 
        SPW = 1:No_states
        SM = 1:3
        SF = 1:No_states
        SD = 1:No_states
        SC = 1:No_states
        SPoW = 1:No_states

        seed = Int(seed_list[p])

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
        
        
        generate_diagram!(diagram,positive_path_utility=true)


        println("starting RJT")
        t = @spawn RJT_model(diagram)
        solutions = [3700,3700,3700,3700,0.0]
        status = "unknown"
        end_time = now() + Dates.Second(3700)
        while now() <= end_time
            sleep(0.01)
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


        println("RJT model building time with jl-package: ",solutions[2])
        println("RJT solve time with jl-package: ", solutions[3])
        println("RJT total: ", solutions[1])
        println("RJT other time: ", solutions[5])
        println("RJT objective value: ", solutions[4])
        println("iteration: ", p)
        println("No states: ", No_states)

        results_RJT[p,1] = solutions[1]
        results_RJT[p,2] = solutions[2]
        results_RJT[p,3] = solutions[3]
        results_RJT[p,4] = solutions[4]
        results_RJT[p,5] = solutions[5]
        if status == "OPTIMAL"
            results_RJT[p,6] = 1.0
        elseif status == "FEASIBLE"
            results_RJT[p,6] = 2.0
        else
            results_RJT[p,6] = 3.0
        end


    end
    writedlm("../Results/water_cvar/results_RJT_$No_states.csv", results_RJT, ',')
end


