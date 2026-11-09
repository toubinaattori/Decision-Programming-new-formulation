using DecisionProgramming
using JuMP, Gurobi, DataStructures, Statistics
using Base.Iterators: product
using DelimitedFiles
using DataFrames
using Dates
using Random
using Distributions


j = parse(Int64, ARGS[1])
seed = readdlm("../Revision_1_results/Data/turbine_seed_$j.csv",',')[1]
Random.seed!(Int(seed))
States_all = 1:7

PTL = Dict{Int64, Vector{Float64}}()
PFH = Dict{Int64, Vector{Float64}}()
PSS = Dict{Int64, Array{Float64,2}}()
PTS = Dict{Int64, Array{Float64,3}}()
PSE = Dict{Int64, Array{Float64,3}}()
PTE = Dict{Int64, Array{Float64,3}}()
PSR = Dict{Int64, Array{Float64,4}}()
PTR = Dict{Int64, Array{Float64,5}}()
PTF = Dict{Int64, Array{Float64,4}}()
U = Dict{Int64, Array{Float64,3}}()
U_shit = Dict{Int64, Dict{Tuple{Int64,Int64,Int64}, Array{Float64}}}()

for No_states in States_all
    PTL[No_states] = zeros(No_states)
    PFH[No_states] = zeros(No_states)
    PSS[No_states] = zeros(No_states,No_states)
    PTS[No_states] = zeros(No_states,No_states,No_states)
    PSE[No_states] = zeros(No_states,No_states,No_states)
    PTE[No_states] = zeros(No_states,No_states,No_states)
    PSR[No_states] = zeros(3,No_states,No_states,No_states)
    PTR[No_states] = zeros(3,No_states,No_states,No_states,No_states)
    PTF[No_states] = zeros(3,No_states,No_states,No_states)
    U_shit[No_states] = Dict{Tuple{Int64,Int64,Int64}, Array{Float64}}()
    for i in 1:3, md in 1:3, tf in 1:No_states
        U_shit[No_states][(i, md, tf)] = []
    end
end


results_RJT = zeros(7,5)

for No_states in States_all
     TL = rand(No_states)
    FH = rand(No_states)
    PTL[No_states] = TL./sum(TL)
    PFH[No_states] = FH./sum(FH)

    # SS: 1=worst, N=best. High FH (worse history) → center near 1; low FH → center near N.
    for i in 1:No_states
        center = No_states + 1 - i
        weights = exp.(-1.0 .* abs.(collect(1:No_states) .- center))
        noise = 0.5 .+ rand(No_states)
        PSS[No_states][i, :] = (weights .* noise) ./ sum(weights .* noise)
    end
    # TS: 1=worst, N=best. Equal contribution from FH and TL.
    # High FH + high TL → center near 1 (worst); low FH + low TL → center near N (best).
    for i in 1:No_states, j in 1:No_states
        center = No_states + 1 - round(Int, (i + j) / 2)
        weights = exp.(-1.0 .* abs.(collect(1:No_states) .- center))
        noise = 0.5 .+ rand(No_states)
        PTS[No_states][i,j, :] = (weights .* noise) ./ sum(weights .* noise)
    end

    # SE peaks at true turbine state (TS=j, 1=worst N=best).
    # Sharpness increases with sensor quality (SS=i: 1=worst, N=best)
    # and slightly increases for better turbine state (TS=j: 1=worst, N=best).
    for i in 1:No_states, j in 1:No_states
        ss_factor = (i - 1) / max(1, No_states - 1)               # 0.0 (worst SS) → 1.0 (best SS)
        ts_factor = 0.7 + 0.3 * (j - 1) / max(1, No_states - 1)   # 0.7 (worst TS) → 1.0 (best TS)
        sharpness = 2.0 * ss_factor * ts_factor + 0.3
        weights = exp.(-sharpness .* abs.(collect(1:No_states) .- i))
        noise = 0.5 .+ rand(No_states)
        PSE[No_states][i,j, :] = (weights .* noise) ./ sum(weights .* noise)
    end

    # TE peaks at average of TS and SE: fuses true turbine state and sensor estimate.
    for i in 1:No_states, j in 1:No_states
        sharpness = 0.5 + 1.5 * (j - 1) / max(1, No_states - 1)   # 0.5 (worst TS) → 2.0 (best TS)
        weights = exp.(-sharpness .* abs.(collect(1:No_states) .- i))
        noise = 0.5 .+ rand(No_states)
        PTE[No_states][i,j, :] = (weights .* noise) ./ sum(weights .* noise)
    end

    # SR: depends on ID(i), SS(j), SE(l) → outcome SR
    # ID=1 (no inspection):        SR = SE deterministically (no new information gained)
    # ID=2 (sensor inspection):    sharp peak at true SS — accurate sensor reading
    # ID=3 (condition inspection): moderate peak at true SS — incidental sensor information
    for j in 1:No_states, l in 1:No_states
        PSR[No_states][1, j, l, :] .= 0.0
        PSR[No_states][1, j, l, l]  = 1.0
        weights = exp.(-2.0 .* abs.(collect(1:No_states) .- j))
        noise = 0.5 .+ rand(No_states)
        PSR[No_states][2, j, l, :] = (weights .* noise) ./ sum(weights .* noise)
        weights = exp.(-1.0 .* abs.(collect(1:No_states) .- j))
        noise = 0.5 .+ rand(No_states)
        PSR[No_states][3, j, l, :] = (weights .* noise) ./ sum(weights .* noise)
    end

    # TR: depends on ID(i), TS(j), SR(k), TE(l) → outcome TR
    # ID=1 (no inspection):        TR = TE deterministically (no new information gained)
    # ID=2 (sensor inspection):    moderate peak at TS; SR has slight influence (25% weight)
    # ID=3 (condition inspection): sharp peak at TS; SR has slight influence (25% weight)
    for j in 1:No_states, k in 1:No_states, l in 1:No_states
        PTR[No_states][1, j, k, l, :] .= 0.0
        PTR[No_states][1, j, k, l, l]  = 1.0
        sharpness1 = 0.5 + 0.5 * (k - 1) / max(1, No_states - 1)   # 0.5 (worst TS) → 1 (best TS)
        sharpness2 = 0.5 + 2.5 * (k - 1) / max(1, No_states - 1)   # 0.5 (worst TS) → 3 (best TS)
        weights = exp.(-sharpness1 .* abs.(collect(1:No_states) .- j))
        noise = 0.5 .+ rand(No_states)
        PTR[No_states][2, j, k, l, :] = (weights .* noise) ./ sum(weights .* noise)
        weights = exp.(-sharpness2 .* abs.(collect(1:No_states) .- j))
        noise = 0.5 .+ rand(No_states)
        PTR[No_states][3, j, k, l, :] = (weights .* noise) ./ sum(weights .* noise)
    end

    # TF: maintenance level shifts center toward best state; higher load broadens the distribution.
    # States: 1=worst, N=best. Maintenance moves center toward higher (better) states.
    # MD=1 (none):    no improvement; high load → uncertain outcome
    # MD=2 (level 1): moderate improvement ~(N-1)/3 states toward best
    # MD=3 (level 2): strong improvement  ~2(N-1)/3 states toward best
    step = max(1, div(No_states - 1, 3))
    for i in 1:3, j in 1:No_states, k in 1:No_states
        improvement = (i - 1) * step
        center      = clamp(k + improvement, 1, No_states)
        load_factor = 1.0 - 0.6 * (j - 1) / max(1, No_states - 1)  # 1.0 (low load) → 0.4 (high load)
        sharpness   = 2.0 * load_factor + 0.3
        weights = exp.(-sharpness .* abs.(collect(1:No_states) .- center))
        noise = 0.5 .+ rand(No_states)
        PTF[No_states][i,j,k,:] = (weights .* noise) ./ sum(weights .* noise)
    end


    U_temp = zeros(3,3,No_states)
    Rewards = sort(rand(100:10000,No_states))
    costs = sort(rand(10:100,3))
    costs2 = sort(rand(500:3000,3))
    scaler = costs[end] + costs2[end] + 1
    for j in 1:No_states, i in 1:3, k in  1:3
        U_temp[i,k,j] = Rewards[j] - costs[i] - costs2[k] + scaler
    end
    U[No_states] = U_temp
end


for No_states in States_all

@info("Calculating parameters")

    STL = 1:No_states
    SH = 1:No_states
    SSS = 1:No_states
    STS = 1:No_states
    SSE = 1:No_states
    STE = 1:No_states
    SID = 1:3
    SSR = 1:No_states
    STR = 1:No_states
    SMD = 1:3
    STF = 1:No_states

    @info("Creating the influence diagram.")
    @info("Time: ", now())
    diagram = InfluenceDiagram()

    add_node!(diagram, ChanceNode("FH", [], ["$i" for i in 1:No_states]))
    add_node!(diagram, ChanceNode("TL", [], ["$i" for i in 1:No_states]))
    add_node!(diagram, ChanceNode("SS", ["FH"], ["$i" for i in 1:No_states]))
    add_node!(diagram, ChanceNode("TS", ["FH","TL"], ["$i" for i in 1:No_states]))
    add_node!(diagram, ChanceNode("SE", ["SS","TS"], ["$i" for i in 1:No_states]))
    add_node!(diagram, ChanceNode("TE", ["TS","SE"], ["$i" for i in 1:No_states]))
    add_node!(diagram, DecisionNode("ID", ["FH","SE","TE"], ["none", "sensor","condition"]))
    add_node!(diagram, ChanceNode("SR", ["ID","SS","SE"], ["$i" for i in 1:No_states]))
    add_node!(diagram, ChanceNode("TR", ["ID","TS","SR","TE"], ["$i" for i in 1:No_states]))
    add_node!(diagram, DecisionNode("MD", ["FH", "TR","SR"], ["none","level 1","level 2"]))
    add_node!(diagram, ChanceNode("TF", ["MD","TL","TS"], ["$i" for i in 1:No_states]))
    add_node!(diagram, ValueNode("U", ["ID","MD","TF"]))

    generate_arcs!(diagram)

     X_V = ProbabilityMatrix(diagram, "TL")
        X_V = PTL[No_states]
        add_probabilities!(diagram, "TL", X_V)

    X_V = ProbabilityMatrix(diagram, "FH")
        X_V = PFH[No_states]
        add_probabilities!(diagram, "FH", X_V)

        X_V = ProbabilityMatrix(diagram, "SS")
        X_V = PSS[No_states]
        add_probabilities!(diagram, "SS", X_V)

        X_V = ProbabilityMatrix(diagram, "TS")
        X_V = PTS[No_states]
        add_probabilities!(diagram, "TS", X_V)

        X_V = ProbabilityMatrix(diagram, "SE")
        X_V = PSE[No_states]
        add_probabilities!(diagram, "SE", X_V)

        X_V = ProbabilityMatrix(diagram, "TE")
        X_V = PTE[No_states]
        add_probabilities!(diagram, "TE", X_V)

        X_V = ProbabilityMatrix(diagram, "SR")
        X_V = PSR[No_states]
        add_probabilities!(diagram, "SR", X_V)

        X_V = ProbabilityMatrix(diagram, "TR")
        X_V = PTR[No_states]
        add_probabilities!(diagram, "TR", X_V)

        X_V = ProbabilityMatrix(diagram, "TF")
        X_V = PTF[No_states]
        add_probabilities!(diagram, "TF", X_V)

        Y_T = UtilityMatrix(diagram, "U")
        Y_T = U[No_states]
        add_utilities!(diagram, "U", Y_T)


    generate_diagram!(diagram,positive_path_utility=true)

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
           #optimizer = optimizer_with_attributes(Gurobi.Optimizer)
            optimizer = optimizer_with_attributes(Gurobi.Optimizer, "Threads" => 1, "Presolve" => 0, "Heuristics" => 0, "NodeMethod" => 1, "Cuts" => 1, "MIPFocus" => 2)
            set_optimizer(model3, optimizer)
        end  
            optimize!(model3)
    end

    results_RJT[No_states,1] = dt9
    results_RJT[No_states,2] = dt7
    results_RJT[No_states,3] = dt4
    results_RJT[No_states,4] = solve_time(model3)
    writedlm("../Revision_1_results/Results_Turbine_RJT/results_turbine_RJT_triton_$(j)_$No_states.csv", results_RJT[No_states, :])
    println("Total time for ", No_states, " states: ", results_RJT[No_states,1], " seconds.")
    println("Model construction time for ", No_states, " states: ", results_RJT[No_states,3], " seconds.")
    println("Optimization time for ", No_states, " states: ", results_RJT[No_states,4], " seconds.")
end




