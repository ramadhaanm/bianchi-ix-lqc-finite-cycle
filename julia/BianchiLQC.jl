"""
Effective LQC ('K' quantization) and classical Bianchi I / II / IX dynamics, generic in the
number type (Float64, Double64, Float128, BigFloat).

Port of derivations/model.py (WrappedModel); the equations are those verified in
derivations/01_symbolic_checks.py, section 5:

    x_i'  = cos th_i (sin th_j + sin th_k)
    th_i' = -S + 1/2 sum_{l != i} x_l' (th_l - th_i) + F_i(x),    constraint S = W(x)

with x_i = ln p_i, th_i = gamma mubar_i K_i and lapse N = gamma sqrt(Delta).
Angles are stored exactly as th_i = n_i pi + d_i (integers n_i, |d_i| <~ pi/2).
State vector u = (x1, x2, x3, d1, d2, d3); the integers n live in the parameter object.
"""
module BianchiLQC

using OrdinaryDiffEqVerner

export Model, initial_data, kasner_from_u, evolve!, evolve_until!, lqc_strength, lnV, dlnV, kasner,
       constraint_error, theta, from_cl, loga, ValidityMonitor, directional_rates

"""Cumulative trajectory monitor. cutmask records crossings of a_i=1,3,10.
Physical and numerical exits are separate; the minimum is sampled at accepted steps
and at threshold roots. Passing this monitor to every segment prevents re-admission.
"""
mutable struct ValidityMonitor{T}
    minimum::T
    worst::T
    cutmask::Int
    reason::Symbol
end
ValidityMonitor(u::AbstractVector{T}) where {T} = ValidityMonitor{T}(minimum(loga(u)), zero(T), 0, :none)

const WRAP = 1.6707963267948966  # pi/2 + 0.1: re-wrap threshold (hysteresis avoids event chatter)

mutable struct Model{T}
    c::T                    # gamma^2 Delta / 4
    nwall::NTuple{3,T}      # Bianchi class A structure constants
    ix::Bool                # Bianchi IX cross terms
    lqc::Bool
    n::Vector{Int}          # current integer offsets of theta / pi
    piT::T                  # pi in the working precision (never use Float64 `pi` in products)
    quant::Symbol           # :K (extrinsic-curvature holonomies) or :A (connection holonomies, Bianchi IX)
    cA::T                   # (1 + gamma^2) / gamma^2: potential factor of the A quantization
    sqD::T                  # sqrt(Delta)
end

function Model{T}(kind::Symbol = :IX; lqc::Bool = true, quant::Symbol = :K) where {T}
    gamma = T(2375) / T(10000)
    delta = 4 * sqrt(T(3)) * T(pi) * gamma
    nwall = kind === :I ? (zero(T), zero(T), zero(T)) :
            kind === :II ? (one(T), zero(T), zero(T)) :
            kind === :IX ? (one(T), one(T), one(T)) : error("unknown Bianchi type $kind")
    quant in (:K, :A) || error("quant must be :K or :A")
    quant === :A && !(lqc && kind === :IX) && error("the A quantization is implemented for LQC Bianchi IX only")
    Model{T}(gamma^2 * delta / 4, nwall, kind === :IX, lqc, zeros(Int, 3), T(pi), quant,
             (1 + gamma^2) / gamma^2, sqrt(delta))
end

# --- A quantization (Wilson-Ewing 2010, eq. 4.2), verified in derivations/03_A_quantization.py ------
# Holonomies of c_i = gamma K_i + Gamma_i(p): sin(mubar_i c_i) = sin(psi_i), psi_i = th_i + phi_i(x),
#   phi_l = (sqrt(Delta)/2) (a_l/(a_j a_k) - a_j/(a_l a_k) - a_k/(a_l a_j)),
#   eta_i = sqrt(Delta) e^{(x_j + x_k - 3 x_i)/2},
#   x_i'  = cos psi_i (sig_i + eta_i),   sig_i = sin psi_j + sin psi_k,
#   th_i' = -S_A + 1/2 sum_{l!=i} x_l' (th_l - th_i) - sum_l x_l' dphi_l/dx_i
#           - sum_l omega_li eta_l sin psi_l + cA F^K_i,     omega_li = -1/2 (l = i), 3/2 (l != i),
#   constraint  S_A + sum_l eta_l sin psi_l = cA W^K.
# Only common shifts th -> th + 2 pi are symmetries (th + pi flips the eta sin psi terms).

"""phi_l(x) and the Jacobian dphi[l, i] = d phi_l / d x_i."""
function phis(m::Model{T}, x) where {T}
    L = loga(x)
    phi = zeros(T, 3); dphi = zeros(T, 3, 3)
    for l in 1:3
        j, k = cyc(l)
        # monomials a_l/(a_j a_k), a_j/(a_l a_k), a_k/(a_l a_j) with signs +, -, -;
        # exponents c on ln a; d ln a_m / d x_i = 1/2 - delta_mi, so dM/dx_i = M (sum(c)/2 - c_i) = M (-1/2 - c_i)
        for (sgn, top) in ((1, l), (-1, j), (-1, k))
            c = (-1, -1, -1) .+ 2 .* ((1, 2, 3) .== top)
            M = exp(c[1] * L[1] + c[2] * L[2] + c[3] * L[3])
            phi[l] += sgn * M
            for i in 1:3
                dphi[l, i] += sgn * M * (-T(1) / 2 - c[i])
            end
        end
    end
    return phi .* (m.sqD / 2), dphi .* (m.sqD / 2)
end

etas(m::Model, x) = ntuple(i -> begin (j, k) = cyc(i); m.sqD * exp((x[j] + x[k] - 3x[i]) / 2) end, 3)

# sin / cos of psi_i = n_i pi + d_i + phi_i without cancellation
@inline function scpsi(m::Model{T}, d, i, ph) where {T}
    sgn = isodd(m.n[i]) ? -one(T) : one(T)
    s, c = sincos(d + ph)
    return sgn * s, sgn * c
end

cyc(i) = (mod1(i + 1, 3), mod1(i + 2, 3))

# sin / cos of th_i = n_i pi + d_i without cancellation
@inline function sc(m::Model{T}, d, i) where {T}
    if m.lqc
        sgn = isodd(m.n[i]) ? -one(T) : one(T)
        s, c = sincos(d)
        return sgn * s, sgn * c
    end
    return m.n[i] * m.piT + d, one(T)
end

function walls(m::Model{T}, x) where {T}
    E = ntuple(i -> begin (j, k) = cyc(i); m.nwall[i] * exp(x[j] + x[k] - 3x[i]) end, 3)
    e = m.ix ? ntuple(i -> begin (j, k) = cyc(i); exp(x[i] - x[j] - x[k]) end, 3) : (zero(T), zero(T), zero(T))
    return E, e
end

function W(m, x)
    E, e = walls(m, x)
    return m.c * (sum(E) - 2sum(e))
end

function rhs!(du, u, m::Model{T}, t) where {T}
    m.quant === :A && return rhs_A!(du, u, m, t)
    x = (u[1], u[2], u[3])
    s1, c1 = sc(m, u[4], 1); s2, c2 = sc(m, u[5], 2); s3, c3 = sc(m, u[6], 3)
    xd = (c1 * (s2 + s3), c2 * (s3 + s1), c3 * (s1 + s2))
    S = s1 * s2 + s2 * s3 + s3 * s1
    E, e = walls(m, x)
    @inbounds for i in 1:3
        j, k = cyc(i)
        mix = zero(T)
        for l in (j, k)
            mix += xd[l] * ((m.n[l] - m.n[i]) * m.piT + (u[3 + l] - u[3 + i]))
        end
        F = m.c * (-2E[i] + 2E[j] + 2E[k] - 4e[i])
        du[i] = xd[i]
        du[3 + i] = -S + mix / 2 + F
    end
    return nothing
end

function rhs_A!(du, u, m::Model{T}, t) where {T}
    x = (u[1], u[2], u[3])
    phi, dphi = phis(m, x)
    eta = etas(m, x)
    sc_ = ntuple(i -> scpsi(m, u[3 + i], i, phi[i]), 3)
    s = ntuple(i -> sc_[i][1], 3)
    xd = ntuple(i -> begin (j, k) = cyc(i); sc_[i][2] * (s[j] + s[k] + eta[i]) end, 3)
    S = s[1] * s[2] + s[2] * s[3] + s[3] * s[1]
    E, e = walls(m, x)
    @inbounds for i in 1:3
        j, k = cyc(i)
        acc = -S + m.cA * m.c * (-2E[i] + 2E[j] + 2E[k] - 4e[i])
        for l in (j, k)
            acc += xd[l] * ((m.n[l] - m.n[i]) * m.piT + (u[3 + l] - u[3 + i])) / 2
        end
        for l in 1:3
            acc -= xd[l] * dphi[l, i] + (l == i ? -T(1) / 2 : T(3) / 2) * eta[l] * s[l]
        end
        du[i] = xd[i]
        du[3 + i] = acc
    end
    return nothing
end

theta(m, u) = [m.n[i] * m.piT + u[3 + i] for i in 1:3]
lnV(u) = (u[1] + u[2] + u[3]) / 2

function xdot(m, u)
    if m.quant === :A
        x = (u[1], u[2], u[3]); phi, _ = phis(m, x); eta = etas(m, x)
        sA = ntuple(i -> scpsi(m, u[3 + i], i, phi[i]), 3)
        return ntuple(i -> begin (j, k) = cyc(i); sA[i][2] * (sA[j][1] + sA[k][1] + eta[i]) end, 3)
    end
    s = ntuple(i -> sc(m, u[3 + i], i), 3)
    ntuple(i -> begin (j, k) = cyc(i); s[i][2] * (s[j][1] + s[k][1]) end, 3)
end
dlnV(m, u) = sum(xdot(m, u)) / 2
directional_rates(m, u) = (xd = xdot(m, u); ntuple(i -> sum(xd) / 2 - xd[i], 3))

function kasner(m, u)
    xd = xdot(m, u)
    H = ntuple(i -> (sum(xd) - 2xd[i]) / 2, 3)
    return H ./ sum(H)
end

# Relative Hamiltonian-constraint residual. The kinetic part is normalized by sum_i sin^2,
# which vanishes only if all three holonomy sines vanish. (An earlier version used
# sum |sin_i sin_j|, which collapses whenever two sines pass through zero together -- this
# happens at bounces -- and reported O(1) "errors" for absolute residuals ~1e-11.)
function constraint_error(m::Model{T}, u) where {T}
    x = (u[1], u[2], u[3])
    E, e = walls(m, x)
    if m.quant === :A
        phi, _ = phis(m, x); eta = etas(m, x)
        s = ntuple(i -> scpsi(m, u[3 + i], i, phi[i])[1], 3)
        S = s[1] * s[2] + s[2] * s[3] + s[3] * s[1]
        lin = sum(eta[l] * s[l] for l in 1:3)
        scale = sum(abs2, s) + sum(abs(eta[l] * s[l]) for l in 1:3) + m.cA * m.c * (sum(E) + 2sum(e))
        return abs(S + lin - m.cA * W(m, x)) / scale
    end
    s = ntuple(i -> sc(m, u[3 + i], i)[1], 3)
    S = s[1] * s[2] + s[2] * s[3] + s[3] * s[1]
    scale = sum(abs2, s) + m.c * (sum(E) + 2sum(e))
    return abs(S - W(m, x)) / scale
end

"""Kasner exponents k at equal p_i = P with th_i = s k_i (s < 0: contracting).
th_3 is fixed in closed form by the constraint, which is linear in sin th_3."""
function initial_data(m::Model{T}, k, P, s) where {T}
    x = fill(log(T(P)), 3)
    th = T[s * T(k[1]), s * T(k[2]), 0]
    m.n .= 0
    if m.quant === :A
        th[3] = theta3_A(m, x, th)
    else
        w = W(m, x)
        th[3] = m.lqc ? asin((w - sin(th[1]) * sin(th[2])) / (sin(th[1]) + sin(th[2]))) :
                        (w - th[1] * th[2]) / (th[1] + th[2])
    end
    return T[x..., th...]
end

"""A quantization: the constraint is linear in sin psi_3, so theta_3 follows in closed form."""
function theta3_A(m::Model, x, th)
    phi, _ = phis(m, x); eta = etas(m, x)
    s1, s2 = sin(th[1] + phi[1]), sin(th[2] + phi[2])
    s3 = (m.cA * W(m, x) - s1 * s2 - eta[1] * s1 - eta[2] * s2) / (s1 + s2 + eta[3])
    return asin(s3) - phi[3]
end

kasner_from_u(u) = (d = 1 + u + u^2; (-u / d, (1 + u) / d, u * (1 + u) / d))

"""
    from_cl(m, ucl, Lam) -> u

Map a classical state in Cornish-Levin variables (ln a, ln b, ln c, d/dtau of each; dt = abc dtau)
to (x, d) after rescaling all scale factors by e^Lam (velocities scale by e^{2 Lam}):
    x1 = beta + gamma, x2 = alpha + gamma, x3 = alpha + beta,   theta_i = gamma sqrt(Delta) alpha_i' / (abc).
Classically Lam is pure gauge; in LQC it sets the distance to the Planck scale (theta ~ e^{-Lam}).
For the LQC model theta_3 is then re-fixed on the exact LQC constraint surface. Resets m.n.
"""
function from_cl(m::Model{T}, ucl, Lam) where {T}
    al = T.(ucl[1:3]) .+ T(Lam)
    vel = T.(ucl[4:6]) .* exp(2 * T(Lam))
    x = T[al[2] + al[3], al[1] + al[3], al[1] + al[2]]
    th = 2sqrt(m.c) .* vel ./ exp(sum(al))          # gamma sqrt(Delta) = 2 sqrt(c)
    m.n .= 0
    if m.quant === :A
        th[3] = theta3_A(m, x, th)
    elseif m.lqc
        s1, s2 = sin(th[1]), sin(th[2])
        th[3] = asin((W(m, x) - s1 * s2) / (s1 + s2))
    end
    return T[x; th]
end

"""ln of the scale factors: ln a_i = (x_j + x_k - x_i)/2."""
loga(u) = ((u[2] + u[3] - u[1]) / 2, (u[1] + u[3] - u[2]) / 2, (u[1] + u[2] - u[3]) / 2)

function rewrap!(integrator)
    m, u = integrator.p, integrator.u
    for i in 1:3
        k = round(Int, u[3 + i] / m.piT)
        u[3 + i] -= k * m.piT
        m.n[i] += k
    end
    # exact common shift th -> th - m pi (K quantization); only even m (2 pi) for the A quantization
    m.n .-= m.quant === :A ? 2round(Int, sum(m.n) / 6) : round(Int, sum(m.n) / 3)
    return nothing
end

"""
    evolve!(m, u, direction; reltol, abstol, tmax) -> (u_new, worst_constraint_error)

Integrate from state `u` (with offsets m.n) until the volume reaches an extremum:
direction = +1 for a minimum (bounce), -1 for a maximum (recollapse).
"""
evolve!(m::Model, u, direction::Int; kw...) = evolve_until!(m, u, dlnV, direction; kw...)

"""Holonomy-based phase marker, not a directional Hubble rate or a classicality proof."""
lqc_strength(m, u) = maximum(abs(sc(m, u[3 + i], i)[1]) for i in 1:3)

"""
    evolve_until!(m, u, g, direction; ...) -> (u_new, worst_constraint_error)

Integrate until the scalar g(m, u) crosses zero upwards (direction = +1) or downwards (-1).
g must be invariant under the re-wrapping of the angles (use sin/cos of theta, not d itself).
"""
function evolve_until!(m::Model{T}, u, g, direction::Int; reltol = T(1e-12), abstol = nothing,
                       tmax = T(1e40), alg = Vern9(), maxiters = 10^6,
                       monitor = nothing) where {T}
    # x_i ~ 10-100: absolute tolerance reltol/100. The angles range from ~1e-10 (near a recollapse)
    # to O(1) (bounce), and single components cross zero; a fixed tiny absolute tolerance then forces
    # needless micro-steps. Their absolute tolerance therefore tracks the size of the angle vector:
    # atol_d = reltol * max_i |d_i| (updated after every step), unless `abstol` is given explicitly.
    adaptive_atol = abstol === nothing
    atol = adaptive_atol ? T[fill(reltol / 100, 3); fill(reltol * maximum(abs, u[4:6]), 3)] : abstol
    worst = Ref(zero(T))
    function sample!(v)
        er = constraint_error(m, v)
        worst[] = max(worst[], er)
        if monitor !== nothing
            mn = minimum(loga(v))
            monitor.minimum = min(monitor.minimum, mn)
            monitor.worst = max(monitor.worst, er)
            for (bit, cutoff) in enumerate((1, 3, 10))
                mn < log(T(cutoff)) && (monitor.cutmask |= 1 << (bit - 1))
            end
            mn < 0 && (monitor.reason = :physical)
            monitor.reason === :none && er >= T(1e-3) && (monitor.reason = :constraint)
        end
    end
    sample!(u)
    monitor !== nothing && monitor.reason !== :none && error("validity exit: $(monitor.reason)")
    stop = ContinuousCallback((u, t, integ) -> g(integ.p, u),
                              direction > 0 ? terminate! : nothing,   # upcrossing (dlnV: volume minimum)
                              direction < 0 ? terminate! : nothing;   # downcrossing (dlnV: maximum)
                              save_positions = (false, false))
    wrap = ContinuousCallback((u, t, integ) -> max(u[4]^2, u[5]^2, u[6]^2) - T(WRAP)^2,
                              rewrap!, nothing; save_positions = (false, false))
    function after_step!(integ)
        sample!(integ.u)
        if monitor !== nothing && monitor.reason !== :none
            terminate!(integ)
        end
        if adaptive_atol
            a = reltol * maximum(abs, view(integ.u, 4:6))
            integ.opts.abstol[4] = integ.opts.abstol[5] = integ.opts.abstol[6] = a
        end
        u_modified!(integ, false)
    end
    track = DiscreteCallback((u, t, integ) -> true, after_step!; save_positions = (false, false))
    callbacks = if monitor === nothing
        CallbackSet(stop, wrap, track)
    else
        # Root searches catch excursions below the cutoffs between accepted steps.
        guards = map(enumerate((1, 3, 10))) do (bit, cutoff)
            function crossed!(integ)
                sample!(integ.u)
                monitor.cutmask |= 1 << (bit - 1)
                if cutoff == 1
                    monitor.reason = :physical
                    terminate!(integ)
                else
                    u_modified!(integ, false)
                end
            end
            ContinuousCallback((v, t, integ) -> minimum(loga(v)) - log(T(cutoff)),
                nothing, crossed!; save_positions=(false, false), interp_points=20,
                abstol=zero(T), reltol=zero(T))
        end
        CallbackSet(stop, wrap, guards..., track)
    end
    prob = ODEProblem(rhs!, copy(u), (zero(T), tmax), m)
    sol = solve(prob, alg; reltol = reltol, abstol = atol, callback = callbacks,
                save_everystep = false, save_start = false, maxiters = maxiters)
    sample!(sol.u[end])
    monitor !== nothing && monitor.reason !== :none && error("validity exit: $(monitor.reason)")
    sol.retcode == ReturnCode.Terminated || error("event not reached: $(sol.retcode)")
    return sol.u[end], worst[]
end

end # module
