using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random
using .Threads



function RJT_model(diagram::InfluenceDiagram)
    dt9 = @elapsed begin
        dt4 = @elapsed begin
            model3 = Model()
            set_time_limit_sec(model3, 3600.0)
            z = DecisionVariables(model3, diagram, names=false)
            variables = RJTVariables(model3, diagram, z, names=false)
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


No = 11

for No_states in 3:No

results_RJT = zeros(50,6)
turbine_seeds = readdlm("../Results/turbine/seeds_$No_states.csv",',')

for p in 1:50
    seed = Int(turbine_seeds[p])
    

         SH = 1:No_states
         SSS = 1:No_states
         STS = 1:No_states
         SSE = 1:No_states
         STE = 1:No_states
         SID = 1:2
         SSR = 1:No_states+1
         STR = 1:No_states+1
         SMD = 1:3
         STF = 1:No_states

                @info("Creating the influence diagram.")
                @info("Time: ", now())
                diagram = InfluenceDiagram()

                add_node!(diagram, ChanceNode("H", [], ["$i" for i in 1:No_states]))
                add_node!(diagram, ChanceNode("SS", ["H"], ["$i" for i in 1:No_states]))
                add_node!(diagram, ChanceNode("TS", ["H"], ["$i" for i in 1:No_states]))
                add_node!(diagram, ChanceNode("SE", ["SS","TS"], ["$i" for i in 1:No_states]))
                add_node!(diagram, ChanceNode("TE", ["TS"], ["$i" for i in 1:No_states]))
                add_node!(diagram, DecisionNode("ID", ["SE","TE"], ["yes", "no"]))
                add_node!(diagram, ChanceNode("SR", ["SE","SS","ID"], ["N/A", ["$i" for i in 1:No_states]...]))
                add_node!(diagram, ChanceNode("TR", ["TE","TS","SR","ID"], ["N/A", ["$i" for i in 1:No_states]...]))
                add_node!(diagram, DecisionNode("MD", ["TR","SR"], ["none", "level 1", "level 2"]))
                add_node!(diagram, ChanceNode("TF", ["MD","TS"], ["$i" for i in 1:No_states]))
                add_node!(diagram, ValueNode("U", ["MD","TF","ID"]))

                generate_arcs!(diagram)
                

                 x_V = rand(Xoshiro(seed),No_states)
                 PSH = x_V./sum(x_V)


                X_V = ProbabilityMatrix(diagram, "H")
                X_V[:] = PSH
                add_probabilities!(diagram, "H", X_V)


                 PSS = zeros(No_states,No_states)
                 x_V = rand(Xoshiro(seed+1),No_states,No_states)
                 for i in 1:No_states
                    PSS[i, :] = x_V[i,:]./sum(x_V[i,:])
                end
                
                X_V = ProbabilityMatrix(diagram, "SS")
                for i in 1:No_states
                    X_V["$i", :] = PSS[i, :] 
                end
                add_probabilities!(diagram, "SS", X_V)

                 PTS = zeros(No_states,No_states)
                 x_V = rand(Xoshiro(seed+2),No_states,No_states)
                 for i in 1:No_states
                    PTS[i, :] = x_V[i,:]./sum(x_V[i,:])
                end

                X_V = ProbabilityMatrix(diagram, "TS")
                for i in 1:No_states
                    X_V["$i", :] = PTS[i, :]
                end
                add_probabilities!(diagram, "TS", X_V)

                 PSE = zeros(No_states,No_states,No_states)
                 x_V = rand(Xoshiro(seed+3),No_states,No_states,No_states)
                 for i in 1:No_states, j in 1:No_states  
                    PSE[i,j, :] = x_V[i,j,:]./sum(x_V[i,j,:])
                end

                X_V = ProbabilityMatrix(diagram, "SE")
                for i in 1:No_states, j in 1:No_states  
                    X_V["$i", "$j", :] = PSE[i,j, :]
                end
                add_probabilities!(diagram, "SE", X_V)

                 PTE = zeros(No_states,No_states)
                 x_V = rand(Xoshiro(seed+4),No_states,No_states)
                 for i in 1:No_states
                    PTE[i, :] = x_V[i,:]./sum(x_V[i,:])
                end


                X_V = ProbabilityMatrix(diagram, "TE")
                for i in 1:No_states
                    X_V["$i", :] = PTE[i, :]
                end
                add_probabilities!(diagram, "TE", X_V)



                 PSR = zeros(No_states,No_states,2,No_states+1)
                 x_V = rand(Xoshiro(seed+5),No_states,No_states,No_states)
                 for i in 1:No_states, j in 1:No_states  
                    PSR[i, j, 1, :] = [0, x_V[i,j,:]./sum(x_V[i,j,:])...]
                    PSR[i, j, 2, :] = [1, zeros(No_states)...]
                end

                X_V = ProbabilityMatrix(diagram, "SR")
                for i in 1:No_states, j in 1:No_states  
                    X_V["$i", "$j", "yes", :] = PSR[i, j, 1, :]
                    X_V["$i", "$j", "no", :] = PSR[i, j, 2, :]
                end
                add_probabilities!(diagram, "SR", X_V)

                 PTR = zeros(No_states,No_states,No_states+1,2,No_states+1)
                 x_V = rand(Xoshiro(seed+6),No_states,No_states,No_states+1,No_states+1)
                 for i in 1:No_states, j in 1:No_states, k in 1:No_states+1
                    PTR[i,j, k, 1, :] = x_V[i,j, k, :]./sum(x_V[i,j, k, :])
                    PTR[i,j, k, 2, :] = [1, zeros(No_states)...]
                end


                X_V = ProbabilityMatrix(diagram, "TR")
                for i in 1:No_states, j in 1:No_states
                    X_V["$i", "$j", "N/A", "yes", :] = PTR[i,j, 1, 1, :]
                    X_V["$i", "$j", "N/A", "no", :] = PTR[i,j, 1, 2, :]
                    for  k in 1:No_states
                        X_V["$i", "$j", "$k", "yes", :] = PTR[i,j, k+1, 1, :]
                        X_V["$i", "$j", "$k", "no", :] = PTR[i,j, k+1, 2, :]
                    end
                end
                add_probabilities!(diagram, "TR", X_V)


                 PTF = zeros(3,No_states,No_states)
                 x_V = rand(Xoshiro(seed+7),3,No_states, No_states)
                 for i in 1:3, j in 1:No_states
                    PTF[i,j, :] = x_V[i,j,:]./sum(x_V[i,j,:])
                end

                X_V = ProbabilityMatrix(diagram, "TF")
                for i in 1:3, j in 1:No_states
                    X_V[i, "$j", :] = PTF[i,j, :]
                end
                add_probabilities!(diagram, "TF", X_V)




                 U = zeros(3,No_states,2)
                 profits = sort(rand(Xoshiro(seed+8),100:1000,No_states))
                 costs = Dict([1, 2, 3] .=> sort(rand(Xoshiro(seed+9),50:90,3)))
                 costs2 = Dict( [1, 2] .=> sort(rand(Xoshiro(seed+10),10:40,2)))
                scaler = abs(minimum(profits) - maximum(collect(values(costs))) - maximum(collect(values(costs2))) - 1)
                 for j in 1:No_states, i in 1:3, k in  1:2
                    U[i,j,k] = profits[j] - costs[i] - costs2[k] + scaler
                end

                Y_T = UtilityMatrix(diagram, "U")
                for j in 1:No_states, i in 1:3, k in  1:2
                    Y_T[i,j,k] = U[i,j,k]
                end
                add_utilities!(diagram, "U", Y_T)



                generate_diagram!(diagram,positive_path_utility=true)


                # dt8 = @elapsed begin
                #     # Calculating DP-parameters
                #     dt2 = @elapsed begin
                #         PU = zeros(No_states, No_states, 2, 1+No_states, 1+No_states, 3)
                #         for se in SSE, te in STE, id in SID, sr in SSR, tr in STR, md in SMD
                #             PU[se,te,id,sr,tr,md] = sum(U[md,tf,id]*PSH[h]*PSS[h, ss]*PTS[h,ts]*PSE[ss,ts,se]*PTE[ts,te]*PSR[se,ss,id,sr]*PTR[te,ts,sr,id,tr]*PTF[md,ts,tf] for h in SH, tf in STF, ss in SSS, ts in STS)
                #         end
                #     end



            


                @info("Time: ", now())
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
        
                @info("Time: ", now())
                println("RJT model building time with jl-package: ",solutions[2])
                println("RJT solve time with jl-package: ", solutions[3])
                println("RJT total: ", solutions[1])
                println("RJT other time: ", solutions[5])
                println("RJT objective value: ", solutions[4])
                println("iteration: ", p)
                println("No states: ", No_states)
                println("RJT status: ", status)
        
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

writedlm("../Results/turbine/RJT_$No_states.csv", results_RJT, ',')

end