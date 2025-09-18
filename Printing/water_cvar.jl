using PrettyTables, DelimitedFiles, Statistics

results_DP = zeros(9,50,6)
for i in 1:9
    results_DP[i,:,:] = readdlm("../Results/water_cvar/results_DP_"*string(i+1)*".csv", ',')
end


results_DP_old = zeros(9,50,6)
for i in 1:9
    if isfile("../Results/water_cvar/results_DP_old_"*string(i+1)*".csv")
        results_DP_old[i,:,:] = readdlm("../Results/water_cvar/results_DP_old_"*string(i+1)*".csv", ',')
    end
end

results_RJT = zeros(9,50,6)
for i in 1:9
    if isfile("../Results/water_cvar/results_RJT_"*string(i+1)*".csv")
        results_RJT[i,:,:] = readdlm("../Results/water_cvar/results_RJT_"*string(i+1)*".csv", ',')
    end
end



mean_total_DP = zeros(9)
mean_total_DP_old = zeros(9)
mean_total_RJT = zeros(9)
mean_sol_time_DP = zeros(9)
mean_sol_time_DP_old = zeros(9)
mean_sol_time_RJT = zeros(9)
no_optimal_RJT = zeros(9)
no_feasible_RJT = zeros(9)
no_optimal_DP_old = zeros(9)
no_feasible_DP_old = zeros(9)
for i in 1:9
    mean_total_DP[i] = round.(mean(results_DP[i,2:50,3]) , digits = 2) 
    res = [results_RJT[i,j,:] for j in 1:50 if results_RJT[i,j,6] <= 1.5 && results_RJT[i,j,6] >= 0.5]
    if !isempty(res)
        mean_total_RJT[i] = round.(mean([r[1] for r in res]), digits = 2)
        mean_sol_time_RJT[i] = round.(mean([r[3] for r in res]), digits = 2)
    end
    res_o = [results_DP_old[i,j,:] for j in 1:50 if results_DP_old[i,j,6] <= 1.5 && results_DP_old[i,j,6] >= 0.5]
    if !isempty(res_o)
        mean_total_DP_old[i] = round.(mean([r[1] for r in res_o]), digits = 2)
        mean_sol_time_DP_old[i] = round.(mean([r[3] for r in res_o]), digits = 2)
    end
    mean_sol_time_DP[i] = round.(mean(results_DP[i,2:50,4]), digits = 2)
    mean_total_RJT[i] = round.(mean(results_RJT[i,2:50,1]), digits = 2)
    mean_sol_time_RJT[i] = round.(mean(results_RJT[i,2:50,3]), digits = 2)
    no_optimal_RJT[i] = count(results_RJT[i,1:50,6] .<= 1.5 .&& results_RJT[i,1:50,6] .> 0.0)
    no_feasible_RJT[i] = count(results_RJT[i,1:50,6] .<= 2.5 .&& results_RJT[i,1:50,6] .> 0.0)
    no_optimal_DP_old[i] = count(results_DP_old[i,1:50,6] .<= 1.5 .&& results_DP_old[i,1:50,6] .> 0.0)
    no_feasible_DP_old[i] = count(results_DP_old[i,1:50,6] .<= 2.5 .&& results_DP_old[i,1:50,6] .> 0.0)
end

blah = [Int(i) for i in 2:10]

column_labels = ["N", "M. Tot DP", "M. Sol DP", "M. Tot DP old", "M. Sol DP old", "No Opt DP old", "No Feas DP old", "M. Tot RJT", "M. Sol RJT", "No Opt RJT", "No Feas RJT"]
results = hcat(blah, mean_total_DP, mean_sol_time_DP, mean_total_DP_old, mean_sol_time_DP_old, no_optimal_DP_old, no_feasible_DP_old, mean_total_RJT, mean_sol_time_RJT, no_optimal_RJT, no_feasible_RJT)
pretty_table(results, header = column_labels)
