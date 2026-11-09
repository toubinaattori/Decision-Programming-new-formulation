# Value of discretization for the continuous turbine inspection and maintenance problem,
# implemented as specified in Section 2.3 of the supplementary material
# ("An Efficient MILP Formulation for Solving Influence Diagrams").
#
# One run = one random continuous instance. For each discretization size N in Table 10:
#   1. discretize the influence diagram using sample paths (conditional frequencies),
#   2. solve the discretized diagram with the RDP formulation,
#   3. evaluate the resulting strategy on the continuous problem by Monte Carlo.
#
# Usage:  julia value_of_discretization_paper.jl <instance_index> [n_estimation_samples] [n_evaluation_samples]
# Output: Results/turbine_paper/{profits,parameters,seed}_<instance_index>.csv
#         profits_<i>.csv holds the average evaluated utility for N = 2,...,6 (one per line).


# Wrapped in a module so that globals left in a REPL session by other scripts
# (e.g. a variable named `keys` or `edges`) cannot shadow the names used here.
module TurbineValueOfDiscretizationPaper

using JuMP, Gurobi, Distributions, Random, DelimitedFiles, Statistics

# ---------------------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------------------
const INSTANCE   = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
const N_ESTIMATE = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100_000   # sample paths for discretization
const N_EVALUATE = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 100_000   # sample paths for strategy evaluation
const SEED       = 20_250_000 + INSTANCE
const OUT_DIR    = joinpath(@__DIR__, "..", "Results", "turbine_paper")

# Table 9 writes N(mu, sigma^2), i.e. the second argument is a VARIANCE. The earlier scripts
# passed these numbers to Normal() directly, i.e. used them as standard deviations.
# false = follow the paper (sigma = sqrt(variance)); true = reproduce the earlier scripts.
const SPREAD_AS_STD = false

# Table 10 breakpoints. State n of every chance node is the interval [edges[n], edges[n+1]].
const BREAKPOINTS = Dict(
    2 => [50.0],
    3 => [50.0, 75.0],
    4 => [25.0, 50.0, 75.0],
    5 => [25.0, 50.0, 75.0, 87.5],
    6 => [25.0, 50.0, 62.5, 75.0, 87.5],
)
const STATES_ALL = sort(collect(keys(BREAKPOINTS)))
edges(N) = [0.0; BREAKPOINTS[N]; 100.0]
midpoints(N) = (e = edges(N); [(e[k] + e[k+1]) / 2 for k in 1:N])
bin(e::Vector{Float64}, x) = clamp(searchsortedfirst(e, x) - 1, 1, length(e) - 1)

const GUROBI_PARAMS = ("OptimalityTol" => 1e-9, "MIPGap" => 0.0)

# ---------------------------------------------------------------------------------------
# Random instance (Section 2.3)
# ---------------------------------------------------------------------------------------
function generate_instance(rng)
    var_ss = rand(rng, Uniform(10, 100))
    var_ts = rand(rng, Uniform(10, 100))
    var_se = rand(rng, Uniform(10, 100))
    tr = sort(rand(rng, Uniform(0.2, 10), 2), rev = true)   # larger variance for IN = sensor check
    tf = sort(rand(rng, Uniform(0.2, 10), 3), rev = true)   # largest for M = none, smallest for M = level 2
    return (
        var_ss = var_ss, var_ts = var_ts, var_se = var_se,
        var_sr = 1.0,                          # SR | IN in {1,2} ~ N(s_SS, 1)
        var_tr = tr,                           # [IN = sensor check, IN = turbine check]
        var_tf = tf,                           # [M = none, M = level 1, M = level 2]
        cost_in = [0.0, rand(rng, Uniform(0, 20)), rand(rng, Uniform(50, 300))],
        cost_m  = [0.0, rand(rng, Uniform(100, 2000)), rand(rng, Uniform(3000, 8000))],
        reward_multiplier = 1 / 100,
    )
end

reward(P, s_tf) = P.reward_multiplier * s_tf^3
utility(P, id, md, s_tf) = reward(P, s_tf) - P.cost_in[id] - P.cost_m[md]

# ---------------------------------------------------------------------------------------
# Continuous conditional distributions (Table 9)
# Decision indices: IN 1 = none, 2 = sensor check, 3 = turbine check
#                   M  1 = none, 2 = level 1,      3 = level 2
# ---------------------------------------------------------------------------------------
sd(v) = SPREAD_AS_STD ? v : sqrt(v)

function tnormal(rng, mu, v, a, b)
    b - a < 1e-9 && return a                       # degenerate interval, e.g. s_TS = 0
    return rand(rng, truncated(Normal(mu, max(sd(v), 1e-6)), a, b))
end

# FH, W, SS, TS, SE, TE (do not depend on decisions)
function sample_upstream(rng, P)
    fh = rand(rng, Uniform(0, 100))
    w  = rand(rng, Uniform(0, 100))
    ss = tnormal(rng, 100 - fh, P.var_ss, 0, 100)
    ts = tnormal(rng, 100 - (fh + w) / 2, P.var_ts, 0, 100)
    se = rand(rng) < ts / 100 ? tnormal(rng, ss, P.var_se, 0, 100) : rand(rng, Uniform(0, 100))
    te = tnormal(rng, ts, (100 - se) / 5, 0, 100)
    return fh, w, ss, ts, se, te
end

# SR and TR given IN
function sample_sr_tr(rng, P, id, ss, ts, se, te)
    id == 1 && return se, te                                    # delta(s_SE), delta(s_TE)
    sr = tnormal(rng, ss, P.var_sr, 0, 100)
    if id == 2
        tr = rand(rng) < sr / 100 ? tnormal(rng, ts, P.var_tr[1], 0, 100) : rand(rng, Uniform(0, 100))
    else
        tr = tnormal(rng, ts, P.var_tr[2], 0, 100)
    end
    return sr, tr
end

# TF given M
function sample_tf(rng, P, md, ts, w)
    md == 1 && return tnormal(rng, ts - 0.02w, P.var_tf[1], 0, ts)
    md == 2 && return tnormal(rng, ts + (100 - ts) / 2 - 0.02w, P.var_tf[2], ts, 100)
    return tnormal(rng, 100.0, P.var_tf[3], ts, 100)
end

# ---------------------------------------------------------------------------------------
# Discretization from sample paths (Section 2.3)
#   P(s_c | s_I(c)) = #{paths with s_c, s_I(c)} / #{paths with s_I(c)}
#   U(s_I(v))       = average utility of the paths with s_I(v)
# Decisions are drawn uniformly at random on the sample paths so that every
# information state of the chance nodes downstream of the decisions is visited.
# ---------------------------------------------------------------------------------------
mutable struct Counts
    FH::Array{Float64,1}; W::Array{Float64,1}
    SS::Array{Float64,2}; TS::Array{Float64,3}
    SE::Array{Float64,3}; TE::Array{Float64,3}
    SR::Array{Float64,4}; TR::Array{Float64,5}
    TF::Array{Float64,4}
    Usum::Array{Float64,3}; Ucnt::Array{Float64,3}
end
Counts(N) = Counts(zeros(N), zeros(N), zeros(N, N), zeros(N, N, N), zeros(N, N, N), zeros(N, N, N),
                   zeros(3, N, N, N), zeros(3, N, N, N, N), zeros(3, N, N, N), zeros(3, 3, N), zeros(3, 3, N))

# normalize along the last dimension; unseen information states get a uniform distribution
function normalize_last(A)
    B = similar(A)
    for I in CartesianIndices(size(A)[1:end-1])
        row = A[I, :]
        s = sum(row)
        B[I, :] = s > 0 ? row ./ s : fill(1 / length(row), length(row))
    end
    return B
end

function discretize(rng, P)
    C = Dict(N => Counts(N) for N in STATES_ALL)
    E = Dict(N => edges(N) for N in STATES_ALL)
    for _ in 1:N_ESTIMATE
        fh, w, ss, ts, se, te = sample_upstream(rng, P)
        id = rand(rng, 1:3)
        sr, tr = sample_sr_tr(rng, P, id, ss, ts, se, te)
        md = rand(rng, 1:3)
        tf = sample_tf(rng, P, md, ts, w)
        u = utility(P, id, md, tf)
        for N in STATES_ALL
            e, c = E[N], C[N]
            bfh, bw, bss, bts, bse, bte, bsr, btr, btf =
                bin(e, fh), bin(e, w), bin(e, ss), bin(e, ts), bin(e, se), bin(e, te), bin(e, sr), bin(e, tr), bin(e, tf)
            c.FH[bfh] += 1; c.W[bw] += 1
            c.SS[bfh, bss] += 1
            c.TS[bfh, bw, bts] += 1
            c.SE[bss, bts, bse] += 1
            c.TE[bts, bse, bte] += 1
            c.SR[id, bss, bse, bsr] += 1
            c.TR[id, bts, bte, bsr, btr] += 1
            c.TF[md, bw, bts, btf] += 1
            c.Usum[id, md, btf] += u; c.Ucnt[id, md, btf] += 1
        end
    end

    D = Dict{Int,Any}()
    for N in STATES_ALL
        c = C[N]
        PSR = normalize_last(c.SR)
        PTR = normalize_last(c.TR)
        # IN = none: SR = SE and TR = TE exactly (Dirac deltas)
        for ss in 1:N, se in 1:N, sr in 1:N
            PSR[1, ss, se, sr] = sr == se ? 1.0 : 0.0
        end
        for ts in 1:N, te in 1:N, sr in 1:N, tr in 1:N
            PTR[1, ts, te, sr, tr] = tr == te ? 1.0 : 0.0
        end
        mid = midpoints(N)
        U = [c.Ucnt[id, md, tf] > 0 ? c.Usum[id, md, tf] / c.Ucnt[id, md, tf] : utility(P, id, md, mid[tf])
             for id in 1:3, md in 1:3, tf in 1:N]
        D[N] = (FH = c.FH ./ sum(c.FH), W = c.W ./ sum(c.W),
                SS = normalize_last(c.SS), TS = normalize_last(c.TS),
                SE = normalize_last(c.SE), TE = normalize_last(c.TE),
                SR = PSR, TR = PTR, TF = normalize_last(c.TF), U = U)
    end
    return D
end

# ---------------------------------------------------------------------------------------
# RDP formulation of the discretized influence diagram
# Information sets: IN <- {FH, SE, TE},  M <- {FH, SR, TR}
# y[fh,se,te,id,sr,tr,md] = probability-free indicator of the observed path,
# PU[...] = sum over unobserved nodes (W, SS, TS, TF) of probability * utility.
# ---------------------------------------------------------------------------------------
function solve_rdp(N, D, env)
    shift = 1.0 - minimum(D.U)                 # utilities must be positive in this formulation
    Upos = D.U .+ shift

    # EU[id, md, w, ts] = sum_tf P(tf | md, w, ts) * U(id, md, tf)
    EU = [sum(D.TF[md, w, ts, tf] * Upos[id, md, tf] for tf in 1:N) for id in 1:3, md in 1:3, w in 1:N, ts in 1:N]

    PU = zeros(N, N, N, 3, N, N, 3)
    for fh in 1:N, w in 1:N, ss in 1:N, ts in 1:N
        p0 = D.FH[fh] * D.W[w] * D.SS[fh, ss] * D.TS[fh, w, ts]
        p0 == 0 && continue
        for se in 1:N, te in 1:N
            p1 = p0 * D.SE[ss, ts, se] * D.TE[ts, se, te]
            p1 == 0 && continue
            for id in 1:3, sr in 1:N, tr in 1:N
                p2 = p1 * D.SR[id, ss, se, sr] * D.TR[id, ts, te, sr, tr]
                p2 == 0 && continue
                for md in 1:3
                    PU[fh, se, te, id, sr, tr, md] += p2 * EU[id, md, w, ts]
                end
            end
        end
    end

    S = Base.OneTo(N)          # OneTo => JuMP builds plain Arrays (linear indexing below)
    bigM = 3 * N^2
    model = direct_model(Gurobi.Optimizer(env))
    for (k, v) in GUROBI_PARAMS
        set_attribute(model, k, v)
    end
    @variable(model, 0 <= y[S, S, S, 1:3, S, S, 1:3] <= 1)
    @variable(model, z_in[S, S, S, 1:3], Bin)      # z_in[fh, se, te, id]
    @variable(model, z_m[S, S, S, 1:3], Bin)       # z_m[fh, sr, tr, md]
    @constraint(model, [fh in S, se in S, te in S], sum(z_in[fh, se, te, :]) == 1)
    @constraint(model, [fh in S, sr in S, tr in S], sum(z_m[fh, sr, tr, :]) == 1)
    @constraint(model, [fh in S, se in S, te in S, id in 1:3],
        sum(y[fh, se, te, id, sr, tr, md] for sr in S, tr in S, md in 1:3) <= bigM * z_in[fh, se, te, id])
    @constraint(model, [fh in S, sr in S, tr in S, md in 1:3],
        sum(y[fh, se, te, id, sr, tr, md] for se in S, te in S, id in 1:3) <= bigM * z_m[fh, sr, tr, md])
    @constraint(model, [fh in S, se in S, te in S, sr in S, tr in S],
        sum(y[fh, se, te, id, sr, tr, md] for id in 1:3, md in 1:3) <= 1)
    @objective(model, Max, sum(PU[i] * y[i] for i in eachindex(PU) if PU[i] != 0))
    optimize!(model)
    termination_status(model) == MOI.OPTIMAL || @warn "N = $N not solved to optimality: $(termination_status(model))"

    strat_in = [argmax(value.(z_in[fh, se, te, :])) for fh in S, se in S, te in S]
    strat_m  = [argmax(value.(z_m[fh, sr, tr, :])) for fh in S, sr in S, tr in S]
    return (strat_in = strat_in, strat_m = strat_m,
            expected_utility = objective_value(model) - shift, solve_time = solve_time(model))
end

# ---------------------------------------------------------------------------------------
# Evaluation on the continuous problem. Common random numbers: on each path, SR/TR are
# sampled for every IN alternative and TF for every M alternative, and every strategy
# then picks its own branch.
# ---------------------------------------------------------------------------------------
function evaluate(rng, P, strategies)
    total = Dict(N => 0.0 for N in STATES_ALL)
    E = Dict(N => edges(N) for N in STATES_ALL)
    for _ in 1:N_EVALUATE
        fh, w, ss, ts, se, te = sample_upstream(rng, P)
        srtr = [sample_sr_tr(rng, P, id, ss, ts, se, te) for id in 1:3]
        tfs = [sample_tf(rng, P, md, ts, w) for md in 1:3]
        for N in STATES_ALL
            e, st = E[N], strategies[N]
            id = st.strat_in[bin(e, fh), bin(e, se), bin(e, te)]
            sr, tr = srtr[id]
            md = st.strat_m[bin(e, fh), bin(e, sr), bin(e, tr)]
            total[N] += utility(P, id, md, tfs[md])
        end
    end
    return Dict(N => total[N] / N_EVALUATE for N in STATES_ALL)
end

# ---------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------
function main()
    rng = Xoshiro(SEED)
    P = generate_instance(rng)
    @info "Instance $INSTANCE (seed $SEED)" P

    @info "Discretizing with $N_ESTIMATE sample paths"
    D = discretize(rng, P)

    env = Gurobi.Env()
    strategies = Dict{Int,Any}()
    for N in STATES_ALL
        @info "Solving the RDP model with N = $N"
        strategies[N] = solve_rdp(N, D[N], env)
        @info "N = $N: discretized expected utility $(strategies[N].expected_utility), solve time $(strategies[N].solve_time) s"
    end

    @info "Evaluating strategies with $N_EVALUATE sample paths"
    profits = evaluate(rng, P, strategies)
    for N in STATES_ALL
        println("N = $N: average utility $(profits[N]), improvement over N = 2: ",
                round(100 * (profits[N] - profits[2]) / abs(profits[2]), digits = 2), " %")
    end

    mkpath(OUT_DIR)
    writedlm(joinpath(OUT_DIR, "profits_$INSTANCE.csv"), [profits[N] for N in STATES_ALL], ',')
    writedlm(joinpath(OUT_DIR, "parameters_$INSTANCE.csv"),
        [P.var_ss, P.var_ts, P.var_se, P.var_sr, P.var_tr[1], P.var_tr[2], P.var_tf[1], P.var_tf[2], P.var_tf[3],
         P.cost_in[2], P.cost_in[3], P.cost_m[2], P.cost_m[3], P.reward_multiplier], ',')
    writedlm(joinpath(OUT_DIR, "seed_$INSTANCE.csv"), [SEED], ',')
end

main()

end # module
