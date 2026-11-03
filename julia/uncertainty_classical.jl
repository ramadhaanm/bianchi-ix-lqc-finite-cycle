# Uncertainty-exponent dimension of the classical Mixmaster basin boundary (Grebogi, McDonald, Ott &
# Yorke, Phys. Lett. A 99 (1983) 415), on the Cornish-Levin Fig. 3 region. Target: D = 1.86 +- 0.01.
#   julia --threads=auto --project=julia julia/uncertainty_classical.jl [N] [lnVmax]
#
# For random x in the region and each eps, x is "uncertain" if the outcome at x differs from the
# outcome at x + eps*e (e a random unit vector). f(eps) ~ eps^alpha and D = 2 - alpha.
# Pairs with an undetermined outcome (label 0) are excluded; their fraction is reported so the
# possible bias can be judged.

include(joinpath(@__DIR__, "CornishLevin.jl"))
using .CornishLevin, DoubleFloats, Printf, Random

N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 3000
lnVmax = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 2.5e5
# precision: "f64" (reltol 1e-14, trustworthy to |ln V| ~ 2.5e5) or "d64" (Double64, reltol 1e-26)
T, rtol = (length(ARGS) >= 3 && ARGS[3] == "d64") ? (Double64, 1e-26) : (Float64, 1e-14)
ur, vr = (1.32, 1.36), (1.28, 1.32)
epss = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 1e-6, 1e-7, 1e-8]

rng = MersenneTwister(2026)
pts = [(ur[1] + rand(rng) * (ur[2] - ur[1]), vr[1] + rand(rng) * (vr[2] - vr[1])) for _ in 1:N]
dirs = [(c = cos(2pi * rand(rng)); s = sqrt(1 - c^2) * rand(rng, (-1, 1)); (c, s)) for _ in 1:N]

lab(u, v) = cl_outcome(u, v; lnVmax = lnVmax, T = T, reltol = rtol)[1]
lab(1.33, 1.3)  # compile
base = zeros(Int, N)
Threads.@threads for k in 1:N
    base[k] = lab(pts[k]...)
end
@printf("N = %d, lnVmax = %.1e, %s (reltol %.0e): undetermined at base points %.3f\n", N, lnVmax, T, rtol,
        count(==(0), base) / N)
flush(stdout)

logs, logf = Float64[], Float64[]
println("     eps       f(eps)   uncertain/valid   excluded")
for eps in epss
    pert = zeros(Int, N)
    Threads.@threads for k in 1:N
        pert[k] = lab(pts[k][1] + eps * dirs[k][1], pts[k][2] + eps * dirs[k][2])
    end
    valid = (base .!= 0) .& (pert .!= 0)
    unc = count(valid .& (base .!= pert))
    f = unc / count(valid)
    @printf("  %8.0e   %8.4f   %6d/%-6d     %.3f\n", eps, f, unc, count(valid), 1 - count(valid) / N)
    unc > 0 && (push!(logs, log(eps)); push!(logf, log(f)))
end
# least-squares slope, with a simple standard error
n = length(logs); mx, my = sum(logs) / n, sum(logf) / n
alpha = sum((logs .- mx) .* (logf .- my)) / sum((logs .- mx) .^ 2)
res = logf .- (my .+ alpha .* (logs .- mx))
se = sqrt(sum(res .^ 2) / (n - 2) / sum((logs .- mx) .^ 2))
@printf("alpha = %.3f +- %.3f   =>   D = 2 - alpha = %.3f +- %.3f   (Cornish & Levin: 1.86 +- 0.01)\n",
        alpha, se, 2 - alpha, se)
