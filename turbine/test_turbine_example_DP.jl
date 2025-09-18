using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random


for No_states in  3:11
    results_DP = zeros(50,6)

    # TODO: Swap these two lines if you want to generate new seeds
    turbine_seeds = readdlm("../Results/turbine/seeds_$No_states.csv",',')
    #turbine_seeds = zeros(50)
    
    for i in 1:50
    
        # TODO: uncomment these two lines if you want to generate your own seeds
        #seed = Dates.year(now())*Dates.day(now())*Dates.hour(now()) + Dates.minute(now()) + Dates.second(now())
        #turbine_seeds[i] = seed

        # TODO:  Use this line if you want to use the same seeds as in the computational experiments
        seed = Int(turbine_seeds[i])

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

                

                 x_V = rand(Xoshiro(seed),No_states)
                 PSH = x_V./sum(x_V)



                 PSS = zeros(No_states,No_states)
                 x_V = rand(Xoshiro(seed+1),No_states,No_states)
                 for i in 1:No_states
                    PSS[i, :] = x_V[i,:]./sum(x_V[i,:])
                end
                

                 PTS = zeros(No_states,No_states)
                 x_V = rand(Xoshiro(seed+2),No_states,No_states)
                 for i in 1:No_states
                    PTS[i, :] = x_V[i,:]./sum(x_V[i,:])
                end


                 PSE = zeros(No_states,No_states,No_states)
                 x_V = rand(Xoshiro(seed+3),No_states,No_states,No_states)
                 for i in 1:No_states, j in 1:No_states  
                    PSE[i,j, :] = x_V[i,j,:]./sum(x_V[i,j,:])
                end

                 PTE = zeros(No_states,No_states)
                 x_V = rand(Xoshiro(seed+4),No_states,No_states)
                 for i in 1:No_states
                    PTE[i, :] = x_V[i,:]./sum(x_V[i,:])
                end



                 PSR = zeros(No_states,No_states,2,No_states+1)
                 x_V = rand(Xoshiro(seed+5),No_states,No_states,No_states)
                 for i in 1:No_states, j in 1:No_states  
                    PSR[i, j, 1, :] = [0, x_V[i,j,:]./sum(x_V[i,j,:])...]
                    PSR[i, j, 2, :] = [1, zeros(No_states)...]
                end


                 PTR = zeros(No_states,No_states,No_states+1,2,No_states+1)
                 x_V = rand(Xoshiro(seed+6),No_states,No_states,No_states+1,No_states+1)
                 for i in 1:No_states, j in 1:No_states, k in 1:No_states+1
                    PTR[i,j, k, 1, :] = x_V[i,j, k, :]./sum(x_V[i,j, k, :])
                    PTR[i,j, k, 2, :] = [1, zeros(No_states)...]
                end





                 PTF = zeros(3,No_states,No_states)
                 x_V = rand(Xoshiro(seed+7),3,No_states, No_states)
                 for i in 1:3, j in 1:No_states
                    PTF[i,j, :] = x_V[i,j,:]./sum(x_V[i,j,:])
                end




                 U = zeros(3,No_states,2)
                 profits = sort(rand(Xoshiro(seed+8),100:1000,No_states))
                 costs = Dict([1, 2, 3] .=> sort(rand(Xoshiro(seed+9),50:90,3)))
                 costs2 = Dict( [1, 2] .=> sort(rand(Xoshiro(seed+10),10:40,2)))
                 scaler = abs(minimum(profits) - maximum(collect(values(costs))) - maximum(collect(values(costs2))) - 1)
                 for j in 1:No_states, i in 1:3, k in  1:2
                    U[i,j,k] = profits[j] - costs[i] - costs2[k] + scaler
                end



            function DP_model(SSE::UnitRange{Int64}, STE::UnitRange{Int64}, SID::UnitRange{Int64}, SSR::UnitRange{Int64}, STR::UnitRange{Int64}, SMD::UnitRange{Int64}, No_states::Int64)
                dt8 = @elapsed begin
                    # Calculating DP-parameters
                    dt2 = @elapsed begin
                        PU = zeros(No_states, No_states, 2, 1+No_states, 1+No_states, 3)
                        for se in SSE, te in STE, id in SID, sr in SSR, tr in STR, md in SMD
                            PU[se,te,id,sr,tr,md] = sum(U[md,tf,id]*PSH[h]*PSS[h, ss]*PTS[h,ts]*PSE[ss,ts,se]*PTE[ts,te]*PSR[se,ss,id,sr]*PTR[te,ts,sr,id,tr]*PTF[md,ts,tf] for h in SH, tf in STF, ss in SSS, ts in STS)
                        end
                    end
    
    
                    # Calculating model build time
                    dt1 = @elapsed begin
                        model2 = Model()
                        set_time_limit_sec(model2, 3600.0)
                        @variable(model2, y[SSE, STE, SID, SSR, STR, SMD], upper_bound = 1.0, lower_bound = 0.0)
                        @variable(model2, z1[SSE,STE,SID], Bin)
                        @variable(model2, z2[STR, SSR, SMD], Bin)
                        
                        @constraint(model2,[se in SSE, te in STE], sum(z1[se,te,id] for id in SID) == 1)
                        @constraint(model2,[tr in STR, sr in SSR], sum(z2[tr,sr,md] for md in SMD) == 1)
    
                        @constraint(model2, [se in SSE, te in STE, id in SID], sum(y[se, te, id, sr, tr, md] for sr in SSR, tr in STR, md in SMD) <= 3*(1+No_states)*(1+No_states)*z1[se,te,id])
                        @constraint(model2, [sr in SSR, tr in STR, md in SMD], sum(y[se, te, id, sr, tr, md] for se in SSE, te in STE, id in SID) <= 3*(1+No_states)*(1+No_states)*z2[tr,sr,md])
    
                        @constraint(model2, [se in SSE, te in STE, sr in SSR, tr in STR],  sum(y[se, te, id, sr, tr, md] for id in SID, md in SMD) <= 1)
    
                        @objective(model2, Max, sum(y[se, te, id, sr, tr, md]*PU[se,te,id,sr,tr,md] for se in SSE, te in STE, sr in SSR, tr in STR, id in SID, md in SMD))
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
            

            println("starting DP")
            t = @async DP_model(SSE, STE, SID, SSR, STR, SMD, No_states)
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
            println("Iteration; ", i)
            println("No_states: ", No_states)
            
            # println("RJT model building time: ",dt3)
            # println("RJT solve time: ", solve_time(model))
            
            
            # println(objective_value(model))
            println(solutions_DP[5])
            
            results_DP[i,1] = solutions_DP[1]
            results_DP[i,2] = solutions_DP[2]
            results_DP[i,3] = solutions_DP[3]
            results_DP[i,4] = solutions_DP[4]
            results_DP[i,5] = solutions_DP[5]
            results_DP[i,6] = solutions_DP[6]
end

writedlm("../Results/turbine/DP_NP_$No_states.csv", results_DP, ',')
# TODO: Uncomment this if you are creating your own seeds
#writedlm("../Results/turbine/seeds_$No_states.csv", turbine_seeds, ',')
end