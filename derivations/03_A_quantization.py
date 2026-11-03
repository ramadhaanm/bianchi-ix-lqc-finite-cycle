"""RQ3: the 'A' (connection-based) loop quantization of Bianchi IX, Wilson-Ewing PRD 82 (2010)
043508 [arXiv:1005.5565], eqs. (2.13), (2.16), (4.2); compared with the 'K' quantization used so far.

Run:  .venv/Scripts/python derivations/03_A_quantization.py

Wilson-Ewing's variables: (c_i, p_i), {c_i, p_j} = 8 pi G gamma delta_ij, lapse N = sqrt(p1 p2 p3),
r_o = 2, positive octant (epsilon = +1). Classical constraint (2.16), vacuum:
    C = -1/(8 pi G gamma^2) [ sum p_i p_j c_i c_j + l_o (p1 p2 c3 + p2 p3 c1 + p3 p1 c2)
                              + (l_o^2/4)(1+gamma^2) (2 sum p_i^2 - sum (p_i p_j / p_k)^2) ]
Spin connection (2.13): Gamma_1 = (eps/r_o)(p2 p3/p1^2 - p2/p3 - p3/p2) (+ cyclic), in units where
A^i = (c^i / l_o) omega^i, so the c-component of Gamma is  l_o * Gamma_i.
Plan: c_i = gamma K_i + l_o Gamma_i(p), with K_i the variables of the K quantization.
"""
import sympy as sp

G, gam, Dl, lo = sp.symbols("G gamma Delta l_o", positive=True)
p = sp.symbols("p1 p2 p3", positive=True)
K = sp.symbols("K1 K2 K3", real=True)
idx = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
q = p[0] * p[1] * p[2]
ro, eps = 2, 1


def report(name, ok):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        raise SystemExit(1)


zero = lambda e: sp.simplify(e) == 0

Gam = [sp.Rational(eps, ro) * (p[j] * p[k] / p[i] ** 2 - p[j] / p[k] - p[k] / p[j]) for (i, j, k) in idx]
c = [gam * K[i] + lo * Gam[i] for i in range(3)]

# 1. canonical: {c_i, c_j} = 0 requires d(l_o Gamma_j)/dp_i symmetric (Gamma is a gradient)
report("spin connection is a gradient: dGamma_j/dp_i = dGamma_i/dp_j",
       all(zero(sp.diff(Gam[j], p[i]) - sp.diff(Gam[i], p[j])) for i in range(3) for j in range(3)))

# 2. classical (2.16) in K variables == K-model classical constraint (lapse sqrt q)
pot = 2 * sum(x**2 for x in p) - sum((p[i] * p[j] / p[k]) ** 2 for (k, i, j) in idx)
C_WE = -1 / (8 * sp.pi * G * gam**2) * (sum(p[i] * p[j] * c[i] * c[j] for (k, i, j) in idx)
                                       + lo * sum(p[i] * p[j] * c[k] for (k, i, j) in idx)
                                       + lo**2 / 4 * (1 + gam**2) * pot)
U_K = (sum((p[j] * p[k] / p[i]) ** 2 for (i, j, k) in idx) - 2 * sum(x**2 for x in p)) / (32 * sp.pi * G)
C_Kcl = -sum(p[j] * p[k] * K[j] * K[k] for (i, j, k) in idx) / (8 * sp.pi * G) + U_K
diff = sp.expand(C_WE - C_Kcl)
lo_val = sp.solve(sp.Eq(diff.subs({K[0]: 0, K[1]: 0, K[2]: 0, p[0]: 1, p[1]: 2, p[2]: 3}), 0), lo)
print("   l_o values matching the K-model normalization:", lo_val)
lo1 = [v for v in lo_val if v.is_positive][0]
report("cross terms in K cancel (only potential terms differ)",
       all(zero(sp.diff(diff, K[i])) for i in range(3)))
report(f"with l_o = {lo1}: classical (2.16) == K-model classical constraint exactly", zero(diff.subs(lo, lo1)))

# 3. effective A-quantization (4.2), written in K variables: sin(mubar_i c_i) = sin(th_i + phi_i)
mub = [sp.sqrt(Dl * p[i] / (p[j] * p[k])) for (i, j, k) in idx]
th = [gam * mub[i] * K[i] for i in range(3)]
phi = [mub[i] * lo1 * Gam[i] for i in range(3)]           # mubar_i * (l_o Gamma_i)
sA = [sp.sin(th[i] + phi[i]) for i in range(3)]
C_A = (-q / (8 * sp.pi * G * gam**2 * Dl) * sum(sA[i] * sA[j] for (k, i, j) in idx)
       - lo1 / (8 * sp.pi * G * gam**2 * sp.sqrt(Dl)) * sum((p[i] * p[j]) ** sp.Rational(3, 2) / sp.sqrt(p[k]) * sA[k]
                                                            for (k, i, j) in idx)
       - lo1**2 / (32 * sp.pi * G * gam**2) * (1 + gam**2) * pot)
e = sp.Symbol("e", positive=True)
lead = sp.series(C_A.subs(Dl, e**2), e, 0, 1).removeO()
report("A-model effective constraint -> classical constraint as Delta -> 0", zero(lead - C_Kcl))

# phase offsets phi_i in terms of scale factors: size of the spin-connection correction
a = sp.symbols("a1:4", positive=True)
to_a = {p[0]: a[1] * a[2], p[1]: a[0] * a[2], p[2]: a[0] * a[1]}
print("   phi_1 =", sp.simplify(phi[0].subs(to_a)))

# 4. symmetries: common shift th -> th + n pi, generated canonically by K_i -> K_i + n pi/(gamma mubar_i)
for n_, name in ((1, "pi"), (2, "2 pi")):
    sh = {K[i]: K[i] + n_ * sp.pi / (gam * mub[i]) for i in range(3)}
    inv = zero(sp.expand_trig(sp.expand(C_A.subs(sh, simultaneous=True) - C_A)))
    print(f"   C_A invariant under common shift th -> th + {name}: {inv}")
    if n_ == 2:
        report("C_A invariant under th -> th + 2 pi (so only 2-pi common shifts may be used numerically)", inv)

# 5. (x, theta) form with lapse N_T = gamma sqrt(Delta), derived by hand, checked against Hamilton's
#    equations of C_A:   psi_i = th_i + phi_i,  s_i = sin psi_i,  sig_i = s_j + s_k,
#      eta_i = sqrt(Delta) e^{(x_j + x_k - 3 x_i)/2},  omega_li = -1/2 (l = i), 3/2 (l != i)
#      x_i'  = cos psi_i (sig_i + eta_i)
#      th_i' = -S_A + 1/2 sum_{l!=i} x_l' (th_l - th_i) - sum_l x_l' dphi_l/dx_i
#              - sum_l omega_li eta_l s_l + ((1+gamma^2)/gamma^2) F^K_i(x)
#      constraint: S_A + sum_l eta_l s_l = ((1+gamma^2)/gamma^2) W^K(x)
#    with F^K, W^K the K-model potential force / constraint function (01_symbolic_checks.py, sec. 5).
dp = [-8 * sp.pi * G * sp.diff(C_A, K[i]) for i in range(3)]
dK = [8 * sp.pi * G * sp.diff(C_A, p[i]) for i in range(3)]
lap = gam * sp.sqrt(Dl) / sp.sqrt(q)
xdot = [dp[i] / p[i] * lap for i in range(3)]
thdot = [(sum(sp.diff(th[i], p[l]) * dp[l] for l in range(3)) + sp.diff(th[i], K[i]) * dK[i]) * lap for i in range(3)]
x = sp.symbols("x1:4", real=True)
T_ = sp.symbols("t1:4", real=True)
sub = {p[i]: sp.exp(x[i]) for i in range(3)}
sub.update({K[i]: T_[i] / (gam * mub[i].subs({p[l]: sp.exp(x[l]) for l in range(3)})) for i in range(3)})
phx = [ph.subs({p[l]: sp.exp(x[l]) for l in range(3)}) for ph in phi]
psi = [T_[i] + phx[i] for i in range(3)]
s_ = [sp.sin(v) for v in psi]
eta = [sp.sqrt(Dl) * sp.exp((x[j] + x[k] - 3 * x[i]) / 2) for (i, j, k) in idx]
E_ = [sp.exp(x[j] + x[k] - 3 * x[i]) for (i, j, k) in idx]
e_ = [sp.exp(x[i] - x[j] - x[k]) for (i, j, k) in idx]
WK = gam**2 * Dl / 4 * (sum(E_) - 2 * sum(e_))
FK = [gam**2 * Dl / 4 * (-2 * E_[i] + 2 * E_[j] + 2 * E_[k] - 4 * e_[i]) for (i, j, k) in idx]
xd = [sp.cos(psi[i]) * (s_[j] + s_[k] + eta[i]) for (i, j, k) in idx]
SA = sum(s_[i] * s_[j] for (k, i, j) in idx)
om = lambda l, i: sp.Rational(-1, 2) if l == i else sp.Rational(3, 2)
thd = [-SA + sum(xd[l] * (T_[l] - T_[i]) for l in range(3) if l != i) / 2
       - sum(xd[l] * sp.diff(phx[l], x[i]) for l in range(3))
       - sum(om(l, i) * eta[l] * s_[l] for l in range(3))
       + (1 + gam**2) / gam**2 * FK[i] for i in range(3)]
import random
random.seed(3)
worst = 0
for _ in range(8):
    vals = {x[i]: sp.Rational(random.randint(0, 900), 100) for i in range(3)}
    vals.update({T_[i]: sp.Rational(random.randint(-300, 300), 100) for i in range(3)})
    vals.update({G: 1, gam: sp.Rational(2375, 10000), Dl: sp.Rational(517, 100)})
    for a_, b_ in zip(xdot + thdot, xd + thd):
        ref = sp.N(a_.subs(sub, simultaneous=True).subs(vals), 40)
        got = sp.N(b_.subs(vals), 40)
        worst = max(worst, abs(ref - got) / (abs(ref) + sp.Float("1e-30")))
report(f"(x, theta) A-model equations == Hamilton's equations (8 random points, max rel. diff {float(worst):.1e})",
       worst < 1e-25)
cons = sp.N(((8 * sp.pi * G * gam**2 * Dl / q) * C_A).subs(sub, simultaneous=True).subs(vals)
            + (SA + sum(eta[l] * s_[l] for l in range(3)) - (1 + gam**2) / gam**2 * WK).subs(vals), 40)
report("constraint: (8 pi G gamma^2 Delta / q) C_A = -(S_A + sum eta s - ((1+g^2)/g^2) W^K)", abs(cons) < 1e-25)

print("\nAll A-quantization checks passed.")
