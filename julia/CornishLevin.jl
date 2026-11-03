"""
Classical Mixmaster in the variables of Cornish & Levin, PRL 78 (1997) 998 [gr-qc/9605029].

    alpha = ln a, beta = ln b, gamma = ln c;  ' = d/dtau with dt = abc dtau (singularity at tau -> infinity)
    (ln a^2)'' = (b^2 - c^2)^2 - a^4   et cyc.                                   (their eq. 1)
    H = a'b' + a'c' + b'c' - (1/4)[a^4 + b^4 + c^4 - 2(a^2b^2 + b^2c^2 + a^2c^2)] = 0   (eq. 2)

(eq. 1 is the Euler-Lagrange equation of eq. 2; checked by hand: 2 alpha'' = (b^2-c^2)^2 - a^4.)

Initial data (eq. 4) on alpha = 0:  beta = Sigma/D, gamma = Sigma (v + uv)/D, D = 1 + v + uv,
velocities proportional to the Kasner indices p = (1+u, -u, u+u^2)/(1+u+u^2) (s = 1). Since
sum_{i<j} p_i p_j = 0, the constraint forces the potential to vanish: with a = 1 >= b, c this is
b + c = 1, which fixes Sigma < 0 uniquely. We integrate towards the singularity (tau forward with
velocities -p, using time-reversal symmetry).

Outcome: once in a Kasner phase with u > umax (=7 in the paper), the axis with the largest
exponent, i.e. the one collapsing fastest towards the singularity: 1 = a, 2 = b, 3 = c.
"""
module CornishLevin

using OrdinaryDiffEqVerner

export cl_initial, cl_outcome, kasner_u

function cl_rhs!(du, u, p, t)
    a2, b2, c2 = exp(2u[1]), exp(2u[2]), exp(2u[3])
    du[1], du[2], du[3] = u[4], u[5], u[6]
    du[4] = ((b2 - c2)^2 - a2^2) / 2
    du[5] = ((c2 - a2)^2 - b2^2) / 2
    du[6] = ((a2 - b2)^2 - c2^2) / 2
    return nothing
end

# The same flow, written for long runs towards the singularity. In Misner time tau the velocities
# decay exponentially (they underflow after ~10 eras) and the epochs grow without bound. Variables:
#   z = (alpha_1..3, w_1..3, l),  w_i = d alpha_i / d lnV (the Kasner exponents, sum w = 1),
#   l = ln|d lnV/dtau|.  With L = ln V, f_i = alpha_i'' (tau), F = sum f, Q = e^{2l}:
#   d alpha_i/dL = w_i,   d w_i/dL = (f_i - w_i F)/Q,   d l/dL = F/Q
# and f_i/Q is evaluated with combined exponents (e^{4 alpha - 2l}, ...), so nothing underflows.
# Clock: sigma = ln(-L), contracting towards L -> -infinity (the intent of the Cornish-Levin
# T = ln ln(1/t)): d/dsigma = -max(-L, 1) d/dL. A pure reparametrization: outcomes unchanged.
function fQ(z)
    a, b, c, l = z[1], z[2], z[3], z[7]
    e(x) = exp(x - 2l)
    f1 = (e(4b) - 2e(2b + 2c) + e(4c) - e(4a)) / 2
    f2 = (e(4c) - 2e(2c + 2a) + e(4a) - e(4b)) / 2
    f3 = (e(4a) - 2e(2a + 2b) + e(4b) - e(4c)) / 2
    return f1, f2, f3
end

function cl_rhs_L!(dz, z, p, t)
    f = fQ(z)
    F = f[1] + f[2] + f[3]
    k = -max(-(z[1] + z[2] + z[3]), one(eltype(z)))
    for i in 1:3
        dz[i] = k * z[3 + i]
        dz[3 + i] = k * (f[i] - z[3 + i] * F)
    end
    dz[7] = k * F
    return nothing
end

"""Constraint in (alpha, w, l): sum_{i<j} w_i w_j - (1/4) V e^{-2l} = 0."""
function cl_constraint_L(z)
    a, b, c, l = z[1], z[2], z[3], z[7]
    e(x) = exp(x - 2l)
    V = e(4a) + e(4b) + e(4c) - 2(e(2a + 2b) + e(2b + 2c) + e(2a + 2c))
    return z[4] * z[5] + z[4] * z[6] + z[5] * z[6] - V / 4
end

potential(u) =(a2 = exp(2u[1]); b2 = exp(2u[2]); c2 = exp(2u[3]);
                a2^2 + b2^2 + c2^2 - 2(a2 * b2 + b2 * c2 + a2 * c2))
cl_constraint(u) = u[4] * u[5] + u[4] * u[6] + u[5] * u[6] - potential(u) / 4

function cl_initial(u, v; s = 1.0, T = Float64)
    u, v, s = T(u), T(v), T(s)
    D = 1 + v + u * v
    w = v + u * v
    # e^{Sigma/D} + e^{Sigma w/D} = 1, unique root Sigma < 0 (bisection, monotone)
    f(S) = exp(S / D) + exp(S * w / D) - 1
    lo, hi = -one(T), zero(T)
    while f(lo) > 0
        lo *= 2
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    S = (lo + hi) / 2
    den = 1 + u + u^2
    p = ((1 + u) / den, -u / den, (u + u^2) / den)
    return T[0, S / D, S * w / D, -s * p[1], -s * p[2], -s * p[3]]   # contracting
end

"""u-parameter of a Kasner triple (any order): u = -p_min / (p_mid + p_min)."""
function kasner_u(pk)
    q = sort(collect(pk))
    return -q[1] / (q[2] + q[1])
end

"""
    cl_outcome(u, v; umax=7, kasner_tol=1e-6, lnVmax=1e5) -> (label, u_final, |constraint|)

label 1/2/3 = axis a/b/c collapsing fastest when the first Kasner phase with u > umax is reached;
0 if not reached before |ln V| = lnVmax. (Long runs accumulate round-off: in Float64 at reltol
1e-14 the constraint stays < 1e-8 down to |ln V| ~ 4e4 and < 1e-6 down to ~2.5e5; use
T = Double64 with reltol ~1e-26 to go deeper.)
"""
function cl_outcome(u0, v0; umax = 7.0, kasner_tol = 1e-6, lnVmax = 1e5, reltol = 1e-14, T = Float64)
    label, ufin, cerr = Ref(0), Ref(NaN), Ref(NaN)
    function check(z, t, integ)
        maximum(abs, fQ(z)) < kasner_tol || return false        # potential negligible: Kasner phase
        pk = (z[4], z[5], z[6])                                  # Kasner exponents
        ku = kasner_u(pk)
        ku > umax || return false
        label[] = argmax(pk); ufin[] = Float64(ku); cerr[] = Float64(abs(cl_constraint_L(z)))
        return true
    end
    cb = DiscreteCallback(check, terminate!; save_positions = (false, false))
    u = cl_initial(u0, v0; T = T)
    vsum = u[4] + u[5] + u[6]
    z = [u[1:3]; u[4:6] ./ vsum; log(abs(vsum))]
    prob = ODEProblem(cl_rhs_L!, z, (zero(T), log(T(lnVmax))))
    solve(prob, Vern9(); reltol = T(reltol), abstol = T(reltol) / 100, callback = cb, save_everystep = false,
          save_start = false, maxiters = 10^7)
    return label[], ufin[], cerr[]
end

end # module
