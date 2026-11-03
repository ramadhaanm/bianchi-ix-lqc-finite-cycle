# Baseline checks for the Julia port, against the Python results in derivations/DERIVATIONS.md.
#   julia --project=julia julia/test_baseline.jl
#
#   A. LQC Bianchi I bounce reproduces Fig. 1 of 2303.16747 (k -> 2/3 - k).
#   C. LQC Bianchi IX: ln V at successive volume extrema in Float64 agrees with the Python
#      values; Double64 runs at two tolerances show how many cycles extended precision resolves.

include(joinpath(@__DIR__, "BianchiLQC.jl"))
using .BianchiLQC
using DoubleFloats, Printf

failures = 0
report(name, ok, detail = "") = (global failures += !ok;
    println(ok ? "[PASS] " : "[FAIL] ", name, isempty(detail) ? "" : "\n       " * detail))

# --- A. Bianchi I bounce ------------------------------------------------------
m = Model{Float64}(:I)
k = (0.954498, 0.232394, -0.186891); k = k ./ sum(k)
u0 = initial_data(m, k, 1e6, -1e-3)
kin = kasner(m, u0)                               # before evolving: kasner() reads the offsets m.n
uA, _ = evolve!(m, u0, +1)                        # to the bounce...
v0 = lnV(u0)
# ...then on until the volume is back to its initial value (reuse evolve! machinery via a loop)
using OrdinaryDiffEqVerner
back = ContinuousCallback((u, t, i) -> lnV(u) - v0, terminate!, nothing)
sol = solve(ODEProblem(BianchiLQC.rhs!, uA, (0.0, 1e40), m), Vern9(); reltol = 1e-12,
            abstol = [fill(1e-14, 3); fill(1e-26, 3)], callback = back, save_everystep = false)
kout = kasner(m, sol.u[end])
report("A. Bianchi I LQC bounce: k -> 2/3 - k (paper Fig. 1)",
       # the rule is exact only asymptotically; at theta ~ 1e-3 O(theta^2) corrections remain (sum k^2 = 1 + 4e-8)
       all(abs.(kout .- (2 / 3 .- kin)) .< 2e-7) && all(abs.(kout .- (-0.287832, 0.434273, 0.853558)) .< 2e-6),
       @sprintf("k_out = (%.6f, %.6f, %.6f)", kout...))

# --- C. Bianchi IX cycles -----------------------------------------------------
function extrema_lnV(T, ncycles, reltol)
    m = Model{T}(:IX)
    u = initial_data(m, kasner_from_u(T(34) / 10), T(10)^12, -T(1) / 1000)
    out, worst = T[], zero(T)
    for c in 1:ncycles, d in (+1, -1)
        u, e = evolve!(m, u, d; reltol = T(reltol))
        worst = max(worst, e)
        push!(out, lnV(u))
    end
    return out, worst
end

python_ref = [33.534, 57.552, 33.932, 57.527, 33.891, 57.236]  # derivations/02_numerical_checks.py, rtol 1e-13
t64 = @elapsed (f64, w64) = extrema_lnV(Float64, 3, 1e-13)
report("C. Float64 extrema match the Python run (first 3 cycles, 3 decimals)",
       all(abs.(f64 .- python_ref) .< 2e-3), "Julia: " * join((@sprintf("%.3f", v) for v in f64), " ") *
       @sprintf("   (%.1f s, max local E_r %.1e)", t64, w64))

ncyc = 8
t1 = @elapsed (d1, _) = extrema_lnV(Double64, ncyc, 1e-24)
t2 = @elapsed (d2, _) = extrema_lnV(Double64, ncyc, 1e-26)
dev = abs.(Float64.(d1 .- d2))
println("\n   Double64, |ln V(reltol 1e-24) - ln V(reltol 1e-26)| at successive extrema (B = bounce, R = recollapse):")
for (i, dv) in enumerate(dev)
    @printf("      %s%-2d  ln V = %8.4f   deviation %.1e\n", isodd(i) ? "B" : "R", (i + 1) ÷ 2, Float64(d2[i]), dv)
end
@printf("   timing: %.1f s and %.1f s for %d cycles\n", t1, t2, ncyc)
report("C'. Double64 agrees with Float64 where Float64 is converged (first 4 extrema, < 1e-8)",
       all(abs.(Float64.(d2[1:4]) .- f64[1:4]) .< 1e-8))

println(failures == 0 ? "\nAll Julia baseline checks passed." : "\n$failures check(s) FAILED.")
exit(failures != 0)
