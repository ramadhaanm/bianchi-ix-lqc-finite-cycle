include(joinpath(@__DIR__, "BianchiLQC.jl"))
include(joinpath(@__DIR__, "CornishLevin.jl"))
using .BianchiLQC, .CornishLevin, Test
include(joinpath(@__DIR__, "Cycles.jl"))

@testset "Pathwise validity regression" begin
    for tol in (1e-12, 1e-13)
        m = Model{Float64}(:IX)
        u = from_cl(m, cl_initial(1.9223751031077165, 1.4036725164488069), 20.0)
        mon = ValidityMonitor(u)
        u, info = run_cycle!(m,u; monitor=mon, reltol=tol)
        @test mon.reason === :none
        @test info.lna_min_B > 8
        @test mon.cutmask & 4 == 0  # first cycle remains above 10 Planck lengths
        @test_throws ErrorException run_cycle!(m,u; monitor=mon, reltol=tol)
        @test mon.reason === :physical
        @test mon.cutmask & 1 == 1
        @test mon.cutmask == 7
        @test mon.worst < 1e-3
        @test abs(mon.minimum) < 1e-6 # termination at the crossing, not after the excursion
    end
    m = Model{Float64}(:IX)
    u = from_cl(m, cl_initial(1.5,1.5), 20.0)
    u[1:3] .-= 100
    mon = ValidityMonitor(u)
    @test_throws ErrorException evolve!(m,u,+1; monitor=mon)
    @test mon.reason === :physical
    m = Model{Float64}(:IX; quant=:A)
    u = from_cl(m, cl_initial(1.5469022866032291, 1.9149896759345253), 20.0)
    mon = ValidityMonitor(u)
    u, _ = run_cycle!(m,u; monitor=mon)
    @test mon.reason === :none
    @test_throws ErrorException run_cycle!(m,u; monitor=mon)
    @test mon.reason === :physical
    @test mon.cutmask == 7
    # A genuine constraint violation: at the LQC onset all holonomy sines are O(0.5); shifting one
    # angle by 0.3 violates S = W by ~10% of the kinetic scale and must stop the run. (A former
    # version shifted an angle by 0.1 in a Kasner state with all sines ~1e-9; the shifted state has
    # two tiny sines and lies almost on the constraint surface, so it is not a meaningful violation.)
    m = Model{Float64}(:IX)
    u = from_cl(m, cl_initial(1.5,1.5),20.0)
    u, _ = to_onset!(m, u)
    u[6] += 0.3
    @test constraint_error(m, u) > 1e-2
    mon = ValidityMonitor(u)
    @test_throws ErrorException evolve!(m,u,+1; monitor=mon)
    @test mon.reason === :constraint
end

@testset "Constraint residual normalization (false-alarm regression)" begin
    # At this K bounce (Lambda = 30) two holonomy sines pass through zero together. The absolute
    # residual stays ~6e-11, but the former normalization (sum of |pairwise products|) collapsed
    # to ~9e-11 and reported a relative error 0.67, terminating a healthy trajectory.
    m = Model{Float64}(:IX)
    u = from_cl(m, cl_initial(1.0606966018464825, 1.7682456879825685), 30.0)
    mon = ValidityMonitor(u)
    u, info = run_cycle!(m, u; monitor = mon)
    @test mon.reason === :none
    @test mon.worst < 1e-6
    @test info.lna_min_B > 8
end
