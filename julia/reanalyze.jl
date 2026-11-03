# Reanalysis v3: raw, resumable trajectory records with pathwise validity and two epsilon controls.
# julia --threads=8 --project=julia julia/reanalyze.jl <phase|labels|long|basin> <L> <K|A> <N> [onset] [budget]
# N is number of pairs (grid side length for basin). Output lives exclusively in results/revised/.
include(joinpath(@__DIR__, "BianchiLQC.jl"))
include(joinpath(@__DIR__, "CornishLevin.jl"))
using .BianchiLQC, .CornishLevin, DoubleFloats, Printf, Random, SHA
include(joinpath(@__DIR__, "Cycles.jl"))

mode, Lam, quant, N = ARGS[1], parse(Int, ARGS[2]), Symbol(ARGS[3]), parse(Int, ARGS[4])
onset = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.5
budget = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 300_000
mode in ("phase", "labels", "long", "basin") || error("unknown mode")
const T = mode == "long" ? Double64 : Float64
ncyc = mode == "long" ? 6 : mode == "basin" ? 2 : 3
rtol = mode == "long" ? T(1e-28) : T(1e-12)
epss = mode == "long" ? [1e-24, 1e-25] : mode == "labels" ?
    [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6,
     3e-7, 1e-7, 3e-8, 1e-8, 3e-9, 1e-9, 3e-10, 1e-10] : [1e-9, 1e-10]
roles = mode == "basin" ? ["base"] : ["base", "ctrl", ["e$(j)" for j in eachindex(epss)]...]
root = joinpath(@__DIR__, "..")
outdir = joinpath(root, "results", "revised"); mkpath(outdir)
tag = "$(mode)_L$(Lam)_$(quant)_N$(N)_on$(onset)_b$(budget)"
path = joinpath(outdir, tag * ".csv")
sourcehash = bytes2hex(sha256(join(read(joinpath(@__DIR__, f), String)
    for f in ("BianchiLQC.jl", "Cycles.jl", "CornishLevin.jl", "reanalyze.jl"))))
meta = "# schema=3 mode=$mode Lambda=$Lam quant=$quant N=$N cycles=$ncyc onset=$onset budget=$budget rtol=$rtol eps=$(join(epss,';')) source=$sourcehash"

# Sampling each point with its own seed makes extension/resumption independent of ensemble size.
function point(k)
    if mode == "basin"
        i, j = divrem(k - 1, N)
        return (T(1) + T(j) / (N-1), T(1) + T(i) / (N-1), zero(T), zero(T))
    end
    rng = MersenneTwister(110000 + k)
    a, b, angle = rand(rng), rand(rng), 2pi * rand(rng)
    (T(1) + T(a), T(1) + T(b), T(cos(angle)), T(sin(angle)))
end

# Each event: ln(a1,a2,a3) in hi/lo pairs, lnV, proxy, directional rates,
# cumulative sampled min ln(a), cumulative constraint residual, cutoff mask.
fields = ["a1_hi", "a1_lo", "a2_hi", "a2_lo", "a3_hi", "a3_lo", "lnV",
    "proxy", "h1", "h2", "h3", "minpath", "error", "cutmask"]
events = ["$(ev)$(c)" for c in 1:ncyc for ev in ("On", "B", "Off", "R", "P")]
header = join(["id", "role", "u", "v", "epsilon", "completed", "reason", "initial_lnV",
    "final_min", "final_error", "reentries", "skipped", ["$(ev)_$(f)" for ev in events for f in fields]...], ',')
hilo(x::Float64) = (x, 0.0)
hilo(x::Double64) = (x.hi, x.lo)

function trajectory(k, role)
    a, b, da, db = point(k)
    epsilon = startswith(role, "e") ? T(epss[parse(Int, role[2:end])]) : zero(T)
    aa, bb = a + epsilon * da, b + epsilon * db
    rt = role == "ctrl" ? rtol / T(10) : rtol
    m = Model{T}(:IX; quant=quant)
    u = from_cl(m, cl_initial(aa, bb; T=T), T(Lam))
    initial_lnV = Float64(lnV(u))
    mon = ValidityMonitor(u)
    rec = fill(NaN, ncyc * 5, length(fields))
    completed, reentries, skipped = 0, 0, 0
    reason = "complete"
    function record!(slot, v)
        la = loga(v); hh = directional_rates(m, v)
        rec[slot, :] .= [hilo(la[1])..., hilo(la[2])..., hilo(la[3])..., Float64(lnV(v)),
            Float64(lqc_strength(m, v)), Float64.(hh)..., Float64(mon.minimum), Float64(mon.worst), mon.cutmask]
    end
    try
        for c in 1:ncyc
            slots = Dict(:On=>5c-4, :B=>5c-3, :Off=>5c-2, :R=>5c-1)
            u, info = run_cycle!(m, u; record=(name,v)->record!(slots[name],v),
                monitor=mon, onset=onset, reltol=rt, maxiters=budget)
            reentries += info.reentries; skipped += info.skipped_max
            completed = c  # R has been reached; P is a separate subsequent segment.
            if mode == "phase"
                lv = lnV(u)
                u, _ = evolve_until!(m, u, (m,v)->lnV(v) - (lv-T(5)), -1;
                    monitor=mon, reltol=rt, maxiters=budget)
                record!(5c, u)
            end
        end
    catch err
        err isa ErrorException || rethrow()
        reason = mon.reason !== :none ? string(mon.reason) : occursin("MaxIters", err.msg) ? "budget" : "integration"
    end
    join([k; role; Float64(a); Float64(b); Float64(epsilon); completed; reason; initial_lnV;
        Float64(mon.minimum); Float64(mon.worst); reentries; skipped; vec(permutedims(rec))], ',')
end

done = Set{Tuple{Int,String}}()
if isfile(path)
    open(path) do io
        readline(io) == meta || error("Refusing to mix incompatible settings or source revisions: $path")
        readline(io) == header || error("CSV schema mismatch")
        for line in eachline(io)
            f = split(line, ',')
            length(f) == length(split(header, ',')) || error("Incomplete CSV row; preserve and repair before resuming")
            push!(done, (parse(Int,f[1]), f[2]))
        end
    end
else
    open(path, "w") do io
        println(io, meta); println(io, header)
    end
end
np = mode == "basin" ? N^2 : N
jobs = [(k,r) for k in 1:np for r in roles if (k,r) ∉ done]
println("$tag: $(length(jobs)) jobs, $(Threads.nthreads()) threads; raw data -> $path"); flush(stdout)
trajectory(1, "base") # compile before concurrent jobs
queue = Channel{Tuple{Int,String}}(length(jobs)); foreach(j->put!(queue,j), jobs); close(queue)
lk = ReentrantLock(); countdone = Threads.Atomic{Int}(0); t0 = time()
open(path, "a") do io
    @sync for _ in 1:Threads.nthreads()
        Threads.@spawn for (k, role) in queue
            line = trajectory(k,role)
            lock(lk) do
                println(io,line); flush(io)
                n = Threads.atomic_add!(countdone,1) + 1
                if n % max(1, length(jobs) ÷ 40) == 0 || n == length(jobs)
                    @printf("%s %d/%d elapsed %.1f min\n", tag,n,length(jobs),(time()-t0)/60); flush(stdout)
                end
            end
        end
    end
end
println("DONE $tag"); flush(stdout)
