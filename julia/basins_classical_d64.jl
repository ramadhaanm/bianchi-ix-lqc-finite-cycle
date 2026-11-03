# Classical Cornish-Levin basins in Double64, deep enough that almost every point gets an outcome
# (Float64 leaves ~25% undetermined at |ln V| <= 1e5; see results/RESULTS.md section 1).
#   julia --threads=auto --project=julia julia/basins_classical_d64.jl <n_grid> [region] [lnVmax]
#     region: fig2 (u, v in [1,2]^2, default) or fig3 (Cornish-Levin zoom); lnVmax default 1e13
# Rows "i,j,label" (label 1/2/3 = axis a/b/c, 0 = undetermined) are appended as points finish, so the
# run is resumable. Output: results/cl_basins_<region>_<n>_d64.csv; the header matches cl_basins_*.csv
# and plot_basins.py / make_figures.py read it after `assemble` (written at the end).

include(joinpath(@__DIR__, "CornishLevin.jl"))
using .CornishLevin, DoubleFloats, Printf

n = parse(Int, ARGS[1])
region = length(ARGS) >= 2 ? ARGS[2] : "fig2"
lnVmax = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e13
ur, vr = region == "fig3" ? ((1.32, 1.36), (1.28, 1.32)) : ((1.0, 2.0), (1.0, 2.0))
rows_path = joinpath(@__DIR__, "..", "results", "cl_basins_$(region)_$(n)_d64.rows")
grid_path = joinpath(@__DIR__, "..", "results", "cl_basins_$(region)_$(n)_d64.csv")

done = Set{Tuple{Int,Int}}()
if isfile(rows_path)
    for l in eachline(rows_path)
        f = split(l, ","); push!(done, (parse(Int, f[1]), parse(Int, f[2])))
    end
end
us = range(ur..., length = n); vs = range(vr..., length = n)
todo = [(i, j) for i in 1:n for j in 1:n if (i, j) ∉ done]
cl_outcome(1.5, 1.5; T = Double64, reltol = 1e-26, lnVmax = lnVmax)   # compile
@printf("classical %s %dx%d Double64, |ln V| <= %.0e: %d points to do (%d done), %d threads\n",
        region, n, n, lnVmax, length(todo), length(done), Threads.nthreads()); flush(stdout)

q = Channel{Tuple{Int,Int}}(length(todo)); foreach(p -> put!(q, p), todo); close(q)
io = open(rows_path, "a"); lk = ReentrantLock(); t0 = time(); cnt = Threads.Atomic{Int}(0)
step = max(1, length(todo) ÷ 50)
@sync for _ in 1:Threads.nthreads()
    Threads.@spawn for (i, j) in q
        lab = cl_outcome(us[j], vs[i]; T = Double64, reltol = 1e-26, lnVmax = lnVmax)[1]
        k = Threads.atomic_add!(cnt, 1) + 1
        lock(lk) do
            println(io, "$i,$j,$lab"); flush(io)
            if k % step == 0 || k == length(todo)
                el = time() - t0
                @printf("  %6d/%d (%5.1f%%)  elapsed %6.1f min  ETA %6.1f min\n", k, length(todo), 100k / length(todo),
                        el / 60, el / 60 * (length(todo) - k) / k); flush(stdout)
            end
        end
    end
end
close(io)

# assemble the label matrix in the format of cl_basins_*.csv
L = fill(-1, n, n)
for l in eachline(rows_path)
    f = parse.(Int, split(l, ",")); L[f[1], f[2]] = f[3]
end
open(grid_path, "w") do out
    println(out, "# u in $(ur), v in $(vr); rows = v ascending, cols = u ascending; 1=a 2=b 3=c 0=undetermined; Double64, |lnV| <= $(lnVmax)")
    for i in 1:n
        println(out, join(L[i, :], ","))
    end
end
@printf("finished in %.1f min; undetermined %.3f; wrote %s\n", (time() - t0) / 60, count(==(0), L) / n^2, grid_path)
