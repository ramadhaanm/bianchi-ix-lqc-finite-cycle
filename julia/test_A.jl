# Checks of the A-quantization implementation (BianchiLQC.rhs_A!) against Wilson-Ewing's effective
# constraint (1005.5565, eq. 4.2) evaluated directly.
#   julia --threads=auto --project=julia julia/test_A.jl

include(joinpath(@__DIR__, "BianchiLQC.jl"))
include(joinpath(@__DIR__, "CornishLevin.jl"))
using .BianchiLQC, .CornishLevin, DoubleFloats, OrdinaryDiffEqVerner, Printf, Random

failures = 0
report(name, ok, detail = "") = (global failures += !ok;
    println(ok ? "[PASS] " : "[FAIL] ", name, isempty(detail) ? "" : "\n       " * detail))

# --- 1. RHS against finite differences of C_A(p, K) (G = 1, l_o = 1, r_o = 2, lapse sqrt(q)) --------
const Tq = Double64
g = Tq(2375) / 10000; Dl = 4sqrt(Tq(3)) * Tq(pi) * g
function C_A(p, K)
    q = p[1] * p[2] * p[3]
    Gam = [(p[2] * p[3] / p[1]^2 - p[2] / p[3] - p[3] / p[2]) / 2,
           (p[1] * p[3] / p[2]^2 - p[1] / p[3] - p[3] / p[1]) / 2,
           (p[1] * p[2] / p[3]^2 - p[1] / p[2] - p[2] / p[1]) / 2]
    mub = [sqrt(Dl * p[1] / (p[2] * p[3])), sqrt(Dl * p[2] / (p[1] * p[3])), sqrt(Dl * p[3] / (p[1] * p[2]))]
    s = [sin(mub[i] * (g * K[i] + Gam[i])) for i in 1:3]
    pot = 2(p[1]^2 + p[2]^2 + p[3]^2) - (p[1] * p[2] / p[3])^2 - (p[2] * p[3] / p[1])^2 - (p[3] * p[1] / p[2])^2
    return -q / (8pi * g^2 * Dl) * (s[1] * s[2] + s[2] * s[3] + s[3] * s[1]) -
           1 / (8pi * g^2 * sqrt(Dl)) * ((p[1] * p[2])^Tq(1.5) / sqrt(p[3]) * s[3] + (p[2] * p[3])^Tq(1.5) / sqrt(p[1]) * s[1] +
                                         (p[3] * p[1])^Tq(1.5) / sqrt(p[2]) * s[2]) -
           (1 + g^2) / (32pi * g^2) * pot
end
function canonical_xth(p, K)   # (x', th') in the time T with lapse gamma sqrt(Delta)
    h = Tq(1e-13)
    dCdK = [(C_A(p, K .+ h .* (1:3 .== i)) - C_A(p, K .- h .* (1:3 .== i))) / 2h for i in 1:3]
    dCdp = [(C_A(p .* (1 .+ h .* (1:3 .== i)), K) - C_A(p .* (1 .- h .* (1:3 .== i)), K)) / (2h * p[i]) for i in 1:3]
    pd = -8pi .* dCdK; Kd = 8pi .* dCdp
    lap = g * sqrt(Dl) / sqrt(p[1] * p[2] * p[3])
    mub = [sqrt(Dl * p[1] / (p[2] * p[3])), sqrt(Dl * p[2] / (p[1] * p[3])), sqrt(Dl * p[3] / (p[1] * p[2]))]
    xd = pd ./ p .* lap
    # th_i = g mub_i K_i;  dth_i/dT = g (mub_i K_i' + K_i dmub_i/dT), dln mub_i = (x_i' - x_j' - x_k')/2
    thd = [g * (mub[i] * Kd[i] * lap + K[i] * mub[i] * (xd[i] - sum(xd) + xd[i]) / 2) for i in 1:3]
    return [xd; thd]
end
m = Model{Tq}(:IX; quant = :A)
rng = MersenneTwister(5); worst = Tq(0)
for _ in 1:20
    x = Tq.(rand(rng, 3) .* 8); th = Tq.(rand(rng, 3) .* 6 .- 3)
    p = exp.(x)
    mub = [sqrt(Dl * p[1] / (p[2] * p[3])), sqrt(Dl * p[2] / (p[1] * p[3])), sqrt(Dl * p[3] / (p[1] * p[2]))]
    K = th ./ (g .* mub)
    ref = canonical_xth(p, K)
    m.n .= round.(Int, th ./ pi); u = Tq[x; th .- m.n .* m.piT]
    du = zeros(Tq, 6); BianchiLQC.rhs!(du, u, m, Tq(0))
    global worst = max(worst, maximum(abs.(du .- ref) ./ (abs.(ref) .+ Tq(1e-20))))
end
report("A-model RHS == Hamilton's equations of Wilson-Ewing eq. (4.2) (20 random points, Double64 finite differences)",
       worst < 1e-12, @sprintf("max relative difference %.1e", Float64(worst)))

# --- 2. constraint along a trajectory through bounces and recollapses (Float64) ------------------------
m = Model{Float64}(:IX; quant = :A)
u = from_cl(m, cl_initial(1.5, 1.5), 20.0)
@printf("   initial constraint error %.1e\n", constraint_error(m, u))
rows = String[]; worstc = 0.0; ok = true
for c in 1:2, d in (+1, -1)
    global u, worstc, ok
    try
        u, e = evolve!(m, u, d; maxiters = 10^7)
        worstc = max(worstc, e)
        la = loga(u)
        push!(rows, @sprintf("%s%d lnV=%7.3f  ln a = %7.3f %7.3f %7.3f", d > 0 ? "B" : "R", c, lnV(u), la...))
    catch err
        ok = false; push!(rows, "failed: " * sprint(showerror, err)); break
    end
end
foreach(r -> println("   ", r), rows)
report("A-model: 2 bounce/recollapse cycles, constraint conserved", ok && worstc < 1e-6, @sprintf("max E_r %.1e", worstc))

# --- 3. far from the Planck scale the A model is classical GR (Cornish-Levin mapping check) -------------
ucl0 = cl_initial(1.5, 1.5)
target = sum(ucl0[1:3]) - 12.0
cb = ContinuousCallback((u, t, i) -> u[1] + u[2] + u[3] - target, nothing, terminate!)
ref = solve(ODEProblem(CornishLevin.cl_rhs!, ucl0, (0.0, 1e6)), Vern9(); reltol = 1e-13, abstol = 1e-14,
            callback = cb, save_everystep = false).u[end]
m = Model{Float64}(:IX; quant = :A)
u = from_cl(m, ucl0, 20.0)
u, _ = evolve_until!(m, u, (m, v) -> lnV(v) - (target + 60.0), -1; reltol = 1e-13)
d = maximum(abs.(loga(u) .- 20.0 .- Tuple(ref[1:3])))
report("A model at Lambda = 20 follows classical GR for 12 e-folds of contraction", d < 1e-8,
       @sprintf("max |d ln a| = %.1e", d))

println(failures == 0 ? "\nAll A-model checks passed." : "\n$failures check(s) FAILED.")
exit(failures != 0)
