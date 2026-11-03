"""Symbolic checks for the classical and effective-LQC Bianchi IX Hamiltonian.

Run:  .venv/Scripts/python derivations/01_symbolic_checks.py

Conventions (Blackmore & Wilson-Ewing 2023, arXiv:2303.16747, 'K' quantization):
    p1 = a2 a3, p2 = a1 a3, p3 = a1 a2,   q = det(spatial metric) = p1 p2 p3
    {K_i, p_j} = 8 pi G delta_ij
    C_H = N * H,  lapse N = sqrt(q) for the "tau" time used in the numerics.
"""
import sympy as sp

G, gam, Dl = sp.symbols("G gamma Delta", positive=True)
a1, a2, a3 = sp.symbols("a1 a2 a3", positive=True)
p1, p2, p3 = sp.symbols("p1 p2 p3", positive=True)
K1, K2, K3 = sp.symbols("K1 K2 K3", real=True)
p = (p1, p2, p3)
K = (K1, K2, K3)
q = p1 * p2 * p3
idx = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]  # (i, j, k) cyclic

ok = lambda e: sp.simplify(e) == 0


def report(name, passed):
    print(f"[{'PASS' if passed else 'FAIL'}] {name}")
    if not passed:
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# 1. Spatial curvature of the Bianchi IX metric  a_i^2 sigma_i^2  on S^3,
#    computed from scratch in Euler coordinates (no structure-constant formula).
# ---------------------------------------------------------------------------
th, ph, ps = sp.symbols("theta phi psi", real=True)
X = (th, ph, ps)
sig = [  # left-invariant 1-forms of SU(2), components along (dtheta, dphi, dpsi)
    (-sp.sin(ps), sp.cos(ps) * sp.sin(th), 0),
    (sp.cos(ps), sp.sin(ps) * sp.sin(th), 0),
    (0, sp.cos(th), 1),
]
A = (a1, a2, a3)
g = sp.Matrix(3, 3, lambda m, n: sum(A[i] ** 2 * sig[i][m] * sig[i][n] for i in range(3)))
ginv = g.inv()
Gam = [[[sum(ginv[l, s] * (sp.diff(g[s, m], X[n]) + sp.diff(g[s, n], X[m]) - sp.diff(g[m, n], X[s]))
             for s in range(3)) / 2 for n in range(3)] for m in range(3)] for l in range(3)]
Ric = sp.zeros(3, 3)
for m in range(3):
    for n in range(3):
        Ric[m, n] = sum(sp.diff(Gam[l][m][n], X[l]) - sp.diff(Gam[l][m][l], X[n])
                        + sum(Gam[l][l][s] * Gam[s][m][n] - Gam[l][n][s] * Gam[s][m][l] for s in range(3))
                        for l in range(3))
R3 = sum(ginv[m, n] * Ric[m, n] for m in range(3) for n in range(3))
R3_formula = -(a1**4 + a2**4 + a3**4 - 2 * (a1**2 * a2**2 + a2**2 * a3**2 + a3**2 * a1**2)) / (2 * a1**2 * a2**2 * a3**2)
# trig factors do not simplify cleanly, so compare at random points (also tests homogeneity)
import random
random.seed(1)
diffs = []
for _ in range(6):
    r = lambda lo, hi: sp.Rational(random.randint(int(lo * 1000), int(hi * 1000)), 1000)  # exact
    pt = {th: r(0.2, 2.9), ps: r(0, 6), ph: r(0, 6), a1: r(0.3, 3), a2: r(0.3, 3), a3: r(0.3, 3)}
    diffs.append(abs(sp.N((R3 - R3_formula).subs(pt), 30)))
report("3R(Bianchi IX) = -(sum a_i^4 - 2 sum a_i^2 a_j^2) / (2 a1^2 a2^2 a3^2)  [6 random points]",
       max(diffs) < 1e-25)
report("isotropic limit: 3R = 3/(2a^2)  (S^3 of radius 2a)",
       ok(R3_formula.subs({a2: a1, a3: a1}) - sp.Rational(3, 2) / a1**2))

# Potential U_IX and the relation 3R = -16 pi G U / q
U_IX = (p2**2 * p3**2 / p1**2 + p1**2 * p3**2 / p2**2 + p1**2 * p2**2 / p3**2
        - 2 * (p1**2 + p2**2 + p3**2)) / (32 * sp.pi * G)
to_a = {p1: a2 * a3, p2: a1 * a3, p3: a1 * a2}
report("3R = -16 pi G U_IX / q", ok(R3_formula - (-16 * sp.pi * G * U_IX / q).subs(to_a)))

# ---------------------------------------------------------------------------
# 2. Classical Hamiltonian constraint (lapse N):
#    C = N [ -(1/(8 pi G sqrt q)) sum_{i} p_j p_k K_j K_k  +  U / sqrt q ]
#    from  N sqrt(q)/(16 pi G) [ (K_ab K^ab - K^2) - 3R ].
# ---------------------------------------------------------------------------
N = sp.Symbol("N", positive=True)
Kin_cl = -sum(p[j] * p[k] * K[j] * K[k] for (i, j, k) in idx) / (8 * sp.pi * G)
C_cl = N * (Kin_cl + U_IX) / sp.sqrt(q)

# ADM check: with K^i_j = diag(K_i / a_i) the kinetic term is sqrt(q)(K^a_b K^b_a - K^2)/(16 pi G)
Kmix = [K[i] / A[i] for i in range(3)]
adm_kin = a1 * a2 * a3 * (sum(x**2 for x in Kmix) - sum(Kmix) ** 2) / (16 * sp.pi * G)
report("kinetic term matches ADM  sqrt(q)(K_ab K^ab - K^2)/16piG",
       ok(adm_kin - (Kin_cl / sp.sqrt(q)).subs(to_a)))
adm_pot = -a1 * a2 * a3 * R3_formula / (16 * sp.pi * G)
report("potential term matches ADM  -sqrt(q) 3R / 16piG", ok(adm_pot - (U_IX / sp.sqrt(q)).subs(to_a)))


def hamilton(C):
    """dp_i/dt = -8piG dC/dK_i,  dK_i/dt = +8piG dC/dp_i."""
    return ([-8 * sp.pi * G * sp.diff(C, K[i]) for i in range(3)],
            [8 * sp.pi * G * sp.diff(C, p[i]) for i in range(3)])


# With N = 1 the classical equations give K_i = da_i/dt: check dp_1/dt = d(a2 a3)/dt
t = sp.Symbol("t", positive=True)
af = [sp.Function(f"a{i+1}", positive=True)(t) for i in range(3)]
dp, dK = hamilton(C_cl.subs(N, 1))
sub_t = {p1: af[1] * af[2], p2: af[0] * af[2], p3: af[0] * af[1],
         K1: af[0].diff(t), K2: af[1].diff(t), K3: af[2].diff(t)}
report("dp1/dt = d(a2 a3)/dt when K_i = da_i/dt (N=1)",
       ok(dp[0].subs(sub_t) - (af[1] * af[2]).diff(t)))

# Kasner solution a_i = t^{k_i} solves the Bianchi I equations (U = 0) for a valid Kasner triple
k1v, k2v, k3v = sp.Rational(-1, 3), sp.Rational(2, 3), sp.Rational(2, 3)
assert k1v + k2v + k3v == 1 and k1v**2 + k2v**2 + k3v**2 == 1
kas = {af[0]: t**k1v, af[1]: t**k2v, af[2]: t**k3v}
C_I = (Kin_cl / sp.sqrt(q)).subs(N, 1)
dpI, dKI = hamilton(C_I)
res = [sp.simplify((dKI[i].subs(sub_t).subs(kas).doit() - af[i].diff(t, 2).subs(kas).doit())) for i in range(3)]
cons = sp.simplify(C_I.subs(sub_t).subs(kas).doit())
report("Kasner (-1/3, 2/3, 2/3) solves Bianchi I: K-equations and constraint",
       all(r == 0 for r in res) and cons == 0)

# Geometric Kasner exponent  k_i = (dln a_i) / sum_j (dln a_j); with K_i = da_i/dt:
#   k_i = p_i K_i / sum_j p_j K_j   (used to set initial data)
kexpr = [(af[i].diff(t) / af[i]) / sum(af[j].diff(t) / af[j] for j in range(3)) for i in range(3)]
pK = [p[i] * K[i] for i in range(3)]
report("k_i = p_i K_i / sum_j p_j K_j  (classical)",
       all(ok(kexpr[i] - (pK[i] / sum(pK)).subs(sub_t)) for i in range(3)))

# ---------------------------------------------------------------------------
# 3. Effective LQC ('K' quantization), lapse N = sqrt(q):
#    C_eff = -(q / (8 pi G gamma^2 Delta)) sum_i sin(th_j) sin(th_k) + U
#    th_i = mubar_i * gamma * K_i,  mubar_i = sqrt(Delta p_i / (p_j p_k))
# ---------------------------------------------------------------------------
mubar = [sp.sqrt(Dl * p[i] / (p[j] * p[k])) for (i, j, k) in idx]
thv = [mubar[i] * gam * K[i] for i in range(3)]
C_eff = -q / (8 * sp.pi * G * gam**2 * Dl) * sum(sp.sin(thv[j]) * sp.sin(thv[k]) for (i, j, k) in idx) + U_IX
C_cl_tau = C_cl.subs(N, sp.sqrt(q))

eps = sp.Symbol("epsilon", positive=True)  # Delta = eps^2 * Delta0
lead = sp.series(C_eff.subs(Dl, eps**2), eps, 0, 2).removeO()
report("Delta -> 0:  C_eff -> C_classical (lapse sqrt q)", ok(lead - C_cl_tau))

# each pair term: mubar_j mubar_k = Delta / p_i  (so gamma and Delta cancel classically)
report("mubar_j mubar_k = Delta / p_i",
       all(ok(mubar[j] * mubar[k] - Dl / p[i]) for (i, j, k) in idx))

# ---------------------------------------------------------------------------
# 4. Discrete symmetry used after each bounce (2303.16747, eq. 3.8):
#    th_i -> th_i + pi  at fixed p, i.e. K_i -> K_i + pi/(gamma mubar_i).
#    It is canonical because pi/(gamma mubar_i) = dF/dp_i with F = 2 pi sqrt(q) / (gamma sqrt Delta).
# ---------------------------------------------------------------------------
F = 2 * sp.pi * sp.sqrt(q) / (gam * sp.sqrt(Dl))
shift = [sp.pi / (gam * mubar[i]) for i in range(3)]
report("shift is a gradient: pi/(gamma mubar_i) = dF/dp_i", all(ok(shift[i] - sp.diff(F, p[i])) for i in range(3)))
Ks = {K[i]: K[i] + shift[i] for i in range(3)}
report("C_eff invariant under the shift", ok(sp.expand_trig(C_eff.subs(Ks, simultaneous=True) - C_eff)))
dp_e, dK_e = hamilton(C_eff)
# equations of motion for the shifted variables: K'_i = K_i + dF/dp_i, so dK'_i/dt = dK_i/dt + sum_j F_ij dp_j/dt
inv = all(ok(sp.expand_trig(dp_e[i].subs(Ks, simultaneous=True) - dp_e[i])) for i in range(3))
inv &= all(ok(sp.expand_trig(dK_e[i].subs(Ks, simultaneous=True)
                             - (dK_e[i] + sum(sp.diff(F, p[i], p[j]) * dp_e[j] for j in range(3)))))
           for i in range(3))
report("equations of motion covariant under the shift (exact symmetry)", inv)

# ---------------------------------------------------------------------------
# 5. Reformulation used in the numerics: x_i = ln p_i, theta_i = gamma mubar_i K_i,
#    constant lapse N_T = gamma sqrt(Delta)  (T = proper time / (gamma sqrt Delta)).
#      x_i'     = cos th_i (sin th_j + sin th_k)
#      th_i'    = -S + (1/2) sum_{l != i} x_l' (th_l - th_i) + F_i(x)
#      S        = sum_pairs sin th_j sin th_k ,   constraint:  S = W(x)
#      W        = (gamma^2 Delta / 4) [ sum_i E_i - 2 sum_i e^{x_i - x_j - x_k} ],  E_i = e^{x_j + x_k - 3 x_i}
#      F_i      = (gamma^2 Delta / 4) [ -2 E_i + 2 E_j + 2 E_k - 4 e^{x_i - x_j - x_k} ]
#    (G = 1 in W and F.)  Checked against Hamilton's equations of C_eff.
# ---------------------------------------------------------------------------
x = sp.symbols("x1 x2 x3", real=True)
T_ = sp.symbols("theta1 theta2 theta3", real=True)
dtau_dT = gam * sp.sqrt(Dl) / sp.sqrt(q)  # lapse ratio N_T / N_tau
xdot = [sp.diff(sp.log(p[i]), p[i]) * dp_e[i] * dtau_dT for i in range(3)]
thdot = [gam * (sp.diff(thv[i], p[0]) * dp_e[0] + sp.diff(thv[i], p[1]) * dp_e[1] + sp.diff(thv[i], p[2]) * dp_e[2]
                + sp.diff(thv[i], K[i]) * dK_e[i]) * dtau_dT / gam for i in range(3)]
to_xth = {p[i]: sp.exp(x[i]) for i in range(3)}
to_xth.update({K[i]: T_[i] / (gam * mubar[i].subs({p[l]: sp.exp(x[l]) for l in range(3)})) for i in range(3)})
S_ = sum(sp.sin(T_[j]) * sp.sin(T_[k]) for (i, j, k) in idx)
E_ = [sp.exp(x[j] + x[k] - 3 * x[i]) for (i, j, k) in idx]
e_ = [sp.exp(x[i] - x[j] - x[k]) for (i, j, k) in idx]
W_ = gam**2 * Dl / 4 * (sum(E_) - 2 * sum(e_))
F_ = [gam**2 * Dl / 4 * (-2 * E_[i] + 2 * E_[j] + 2 * E_[k] - 4 * e_[i]) for (i, j, k) in idx]
xd_claim = [sp.cos(T_[i]) * (sp.sin(T_[j]) + sp.sin(T_[k])) for (i, j, k) in idx]
thd_claim = [-S_ + sum(xd_claim[l] * (T_[l] - T_[i]) for l in range(3) if l != i) / 2 + F_[i] for i in range(3)]


def same(a, b):
    d = sp.simplify(sp.powsimp(sp.expand(sp.powdenest((a - b).subs(to_xth, simultaneous=True).subs(G, 1), force=True)), force=True))
    return d == 0


report("x_i' = cos th_i (sin th_j + sin th_k)", all(same(xdot[i], xd_claim[i]) for i in range(3)))
report("th_i' = -S + (1/2) sum_{l!=i} x_l'(th_l - th_i) + F_i", all(same(thdot[i], thd_claim[i]) for i in range(3)))
report("constraint: (8 pi G gamma^2 Delta / q) C_eff = -(S - W)",
       same(8 * sp.pi * G * gam**2 * Dl / q * C_eff, -(S_ - W_)))

print("\nAll symbolic checks passed.")
