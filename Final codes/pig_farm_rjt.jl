using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics, Random
using Base.Iterators: product
using DelimitedFiles
using Dates
using Random

global j = parse(Int64, ARGS[1])
readseed = readdlm("../Revision_1_results/Data/pig_farm_seed_gurobi_$j.csv",',')[1]
seed = Int(readseed)
Random.seed!(seed)

function DP_model(No_periods::Int64, scaler::Float64, PH0::Vector{Float64}, PH::Array{Float64}, PT1::Matrix{Float64}, PT2::Matrix{Float64}, PT3::Matrix{Float64}, U::Vector{Float64}, C::Vector{Float64})

    dt8 = @elapsed begin
        dt2 = @elapsed begin
            diagram = InfluenceDiagram()

            add_node!(diagram, ChanceNode("H1", [], ["ill", "healthy"]))
            for t in 1:No_periods
                add_node!(diagram, ChanceNode("X1_$t", ["H$t"], ["neg", "pos"]))
                add_node!(diagram, ChanceNode("X2_$t", ["H$t"], ["neg", "pos"]))
                add_node!(diagram, ChanceNode("X3_$t", ["H$t"], ["neg", "pos"]))
                add_node!(diagram, DecisionNode("D$t", ["X1_$t", "X2_$t", "X3_$t"], ["treat", "no_treat"]))
                add_node!(diagram, ValueNode("C$t", ["D$t"]))
                add_node!(diagram, ChanceNode("H$(t+1)", ["H$t", "D$t"], ["ill", "healthy"]))
            end
            add_node!(diagram, ValueNode("U_end", ["H$(No_periods+1)"]))

            generate_arcs!(diagram)

            X_H1 = ProbabilityMatrix(diagram, "H1")
            X_H1[:] = PH0
            add_probabilities!(diagram, "H1", X_H1)

            for t in 1:No_periods
                X_T1 = ProbabilityMatrix(diagram, "X1_$t")
                X_T1["ill", :] = PT1[1, :]
                X_T1["healthy", :] = PT1[2, :]
                add_probabilities!(diagram, "X1_$t", X_T1)

                X_T2 = ProbabilityMatrix(diagram, "X2_$t")
                X_T2["ill", :] = PT2[1, :]
                X_T2["healthy", :] = PT2[2, :]
                add_probabilities!(diagram, "X2_$t", X_T2)

                X_T3 = ProbabilityMatrix(diagram, "X3_$t")
                X_T3["ill", :] = PT3[1, :]
                X_T3["healthy", :] = PT3[2, :]
                add_probabilities!(diagram, "X3_$t", X_T3)

                X_H = ProbabilityMatrix(diagram, "H$(t+1)")
                X_H["ill", "treat", :] = PH[1, 1, :]
                X_H["ill", "no_treat", :] = PH[1, 2, :]
                X_H["healthy", "treat", :] = PH[2, 1, :]
                X_H["healthy", "no_treat", :] = PH[2, 2, :]
                add_probabilities!(diagram, "H$(t+1)", X_H)

                Y_C = UtilityMatrix(diagram, "C$t")
                Y_C["treat"] = C[1]
                Y_C["no_treat"] = C[2]
                add_utilities!(diagram, "C$t", Y_C)
            end

            Y_U = UtilityMatrix(diagram, "U_end")
            Y_U["ill"] = U[1] + scaler
            Y_U["healthy"] = U[2] + scaler
            add_utilities!(diagram, "U_end", Y_U)

            generate_diagram!(diagram)
        end

        dt1 = @elapsed begin
            model = Model()
            z = DecisionVariables(model, diagram, names=false)
            variables = RJTVariables(model, diagram, z, names=false)
            EV = expected_value(model, diagram, variables)
            @objective(model, Max, EV)
        end
        println("DP model building time: ", dt1+dt2)
        dt5 = @elapsed begin
            @info("Starting the optimization process.")
            optimizer = optimizer_with_attributes(Gurobi.Optimizer)
            set_optimizer(model, optimizer)
        end

        optimize!(model)
    end

    return (dt2, dt1, dt8, dt5, solve_time(model), objective_value(model))
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
    println("DP objective value: ", solutions_DP[6])
    println("DP total: ",  solutions_DP[3])

    println("Tests: ", No_periods)

    writedlm("../Revision_1_results/Results_Pig_Farm_RJT/results_Pig_Farm_RJT_triton_gurobi_$(j)_$No_periods.csv", solutions_DP)

end
