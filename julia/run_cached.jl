# Same production runner and algorithm, with bit-identical immutable coefficient memoization.
# Arguments are exactly those of reanalyze.jl. The simulation-source hash remains the
# reference algorithm's hash; the execution cache hash is recorded explicitly in the log.
include(joinpath(@__DIR__, "VernerDoubleCache.jl"))
using SHA
println("Execution cache: bit-identical Vern9 Double64 tables; SHA256=" *
    bytes2hex(sha256(read(joinpath(@__DIR__, "VernerDoubleCache.jl")))))
flush(stdout)
include(joinpath(@__DIR__, "reanalyze.jl"))
