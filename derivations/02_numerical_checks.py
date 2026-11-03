"""Numerical baseline: reproduce known results before doing anything new.

  0. The hand-coded (x, theta) right-hand side equals Hamilton's equations of the
     canonical effective Hamiltonian (sympy-generated), at random phase-space points.
  A. Effective LQC Bianchi I bounce: k_i -> 2/3 - k_i, with the exact numbers of
     Fig. 1 of Blackmore & Wilson-Ewing (arXiv:2303.16747).
  B. Classical Bianchi II wall bounce (BKL rule, eq. 2.6 of the same paper).
  C. Effective LQC Bianchi IX: several bounce/recollapse cycles with control of
     the Hamiltonian constraint (the paper needed 128-bit arithmetic for this).

Run:  .venv/Scripts/python derivations/02_numerical_checks.py
"""
import os
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import DELTA, GAMMA, Model, State, WrappedModel, kasner_from_u

RTOL = 1e-12
failures = 0


def report(name, passed, detail=""):
    global failures
    failures += not passed
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f"\n       {detail}" if detail else ""), flush=True)


def integrate(m, y, T0, event, T_max=1e12):
    """Integrate until `event` fires; return final state, time and worst constraint error."""
    event.terminal = True
    sol = solve_ivp(m.rhs, (T0, T0 + T_max), y, method="DOP853", rtol=RTOL, atol=1e-14, events=event)
    if sol.status != 1:
        raise RuntimeError(f"event not reached: {sol.message}")
    err = max(m.constraint_error(sol.y[:, i]) for i in range(sol.y.shape[1]))
    return sol.y[:, -1], sol.t[-1], err


# ---------------------------------------------------------------------------
# 0. Hand-coded RHS vs. canonical Hamilton equations
# ---------------------------------------------------------------------------
P = sp.symbols("p1 p2 p3", positive=True)
Ks = sp.symbols("K1 K2 K3", real=True)
idx = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
q = P[0] * P[1] * P[2]
mub = [sp.sqrt(DELTA * P[i] / (P[j] * P[k])) for (i, j, k) in idx]
th = [GAMMA * mub[i] * Ks[i] for i in range(3)]
U = (sum((P[j] * P[k] / P[i]) ** 2 for (i, j, k) in idx) - 2 * sum(x**2 for x in P)) / (32 * sp.pi)
C = -q / (8 * sp.pi * GAMMA**2 * DELTA) * sum(sp.sin(th[j]) * sp.sin(th[k]) for (i, j, k) in idx) + U
lapse = GAMMA * np.sqrt(DELTA) / sp.sqrt(q)  # N_T / N_tau, with C written for N_tau = sqrt(q)
dp = [-8 * sp.pi * sp.diff(C, Ks[i]) * lapse for i in range(3)]
dK = [8 * sp.pi * sp.diff(C, P[i]) * lapse for i in range(3)]
xd = [dp[i] / P[i] for i in range(3)]
thd = [sum(sp.diff(th[i], P[l]) * dp[l] for l in range(3)) + sp.diff(th[i], Ks[i]) * dK[i] for i in range(3)]
canon = sp.lambdify([P, Ks], xd + thd, "mpmath")
m = Model("IX", lqc=True)
rng = np.random.default_rng(0)
worst = 0.0
for _ in range(20):
    xv = rng.uniform(0, 12, 3)
    tv = rng.uniform(-3, 3, 3)
    pv = np.exp(xv)
    Kv = tv / (GAMMA * np.sqrt(DELTA * pv / (np.roll(pv, -1) * np.roll(pv, -2))))
    ref = np.array([float(v) for v in canon(list(pv), list(Kv))])
    got = m.rhs(0, np.concatenate([xv, tv]))
    worst = max(worst, np.max(np.abs(got - ref) / (np.abs(ref) + 1e-300)))
report("0. (x, theta) RHS == canonical Hamilton equations (20 random points)", worst < 1e-9,
       f"max relative difference {worst:.1e}")
wm, worst = WrappedModel("IX", lqc=True), 0.0
for _ in range(20):
    yv = np.concatenate([rng.uniform(0, 12, 3), rng.uniform(-7, 7, 3)])
    st = State(yv[:3], yv[3:])
    worst = max(worst, np.max(np.abs(wm.rhs_nd(0, np.concatenate([st.x, st.d]), st.n) - m.rhs(0, yv))
                              / np.abs(m.rhs(0, yv))))
report("0'. exact (x, n, d) representation gives the same RHS", worst < 1e-12, f"max relative difference {worst:.1e}")

# ---------------------------------------------------------------------------
# A. LQC Bianchi I bounce
# ---------------------------------------------------------------------------
m = Model("I", lqc=True)
k_in = np.array([0.954498, 0.232394, -0.186891])
y0 = m.initial_data(k_in / k_in.sum(), P=1e6, s=-1e-3)  # contracting towards the bounce
v0 = m.lnV(y0)
back = lambda T, y: m.lnV(y) - v0  # volume back to its initial value after the bounce
back.direction = +1
y1, T1, err = integrate(m, y0, 0.0, back)
k_out = m.kasner(y1)
expected = np.array([-0.287832, 0.434273, 0.853558])
report("A. Bianchi I LQC bounce: k_in -> 2/3 - k_in  (paper Fig. 1)",
       # the rule is exact only asymptotically; at theta ~ 1e-3 O(theta^2) corrections remain (sum k^2 = 1 + 4e-8)
       np.allclose(k_out, 2 / 3 - m.kasner(y0), rtol=0, atol=2e-7) and np.allclose(k_out, expected, rtol=0, atol=2e-6),
       f"k_in={np.round(m.kasner(y0), 6)}  k_out={np.round(k_out, 6)}  paper={expected}  max E_r={err:.1e}")

# ---------------------------------------------------------------------------
# B. Classical Bianchi II: one Kasner transition off the a1^4 wall
# ---------------------------------------------------------------------------
m = Model("II", lqc=False)
y0 = m.initial_data(np.array([6 / 7, -2 / 7, 3 / 7]), P=1e14, s=+1.0)  # expanding, wall negligible
k_in = m.kasner(y0)
k1 = k_in[0]
k_exp = np.array([-k1, k_in[1] + 2 * k1, k_in[2] + 2 * k1]) / (1 + 2 * k1)
v0 = m.lnV(y0)
far = lambda T, y: m.lnV(y) - v0 - 45  # wall hit after ~10 e-folds, negligible again (<1e-15) after ~37
far.direction = +1
y1, T1, err = integrate(m, y0, 0.0, far, T_max=1e60)  # classical: ln V ~ ln t
k_out = m.kasner(y1)
report("B. Bianchi II classical: BKL transition rule k1 -> -k1/(1+2k1), ...",
       np.allclose(k_out, k_exp, rtol=0, atol=1e-9),
       f"k_in={np.round(k_in, 10)}  k_out={np.round(k_out, 10)}  rule={np.round(k_exp, 10)}  max E_r={err:.1e}")

# ---------------------------------------------------------------------------
# C. LQC Bianchi IX: cycles of bounce (min volume) and recollapse (max volume)
# ---------------------------------------------------------------------------
# The local constraint error is a poor criterion here: some bounces coincide with a curvature-wall
# bounce in a nearly frozen configuration (two sines ~ 1e-6), which inflates the relative error.
# The meaningful test is convergence of the trajectory (volume at every extremum) under rtol.
m = WrappedModel("IX", lqc=True)
extremum = lambda T, z, n: m.dlnV_nd(z, n)
n_cycles, runs = 6, {}
for rtol in (1e-12, 1e-13):
    st = State(*np.split(m.initial_data(kasner_from_u(3.4), P=1e12, s=-1e-3), 2))  # contracting Kasner epoch
    rows, worst = [], 0.0
    for c in range(n_cycles):
        for direction, kind in ((+1, "bounce"), (-1, "recollapse")):
            worst = max(worst, m.evolve(st, extremum, direction=direction, rtol=rtol))
            rows.append((kind, m.lnV(st.y), st.y[3:]))
    runs[rtol] = (rows, worst)

rows, worst = runs[1e-13]
dev = np.array([abs(a[1] - b[1]) for a, b in zip(runs[1e-12][0], rows)])
print("\n   C. LQC Bianchi IX cycle log (volume extrema, rtol=1e-13) and |ln V(1e-12) - ln V(1e-13)|:")
print("      type         ln V      deviation   theta_i")
for (kind, lv, th), dv in zip(rows, dev):
    print(f"      {kind:10s} {lv:8.3f}   {dv:9.1e}   {np.array2string(th, precision=4)}")
growth = dev[2:8] / np.maximum(dev[1:7], 1e-16)
report(f"C. Bianchi IX LQC: {n_cycles} bounce/recollapse cycles; first 4 extrema converged to < 1e-8",
       len(rows) == 2 * n_cycles and dev[:4].max() < 1e-8 and worst < 1e-3,
       f"max local E_r = {worst:.1e}; deviation grows ~x{np.exp(np.mean(np.log(growth))):.0f} per extremum "
       f"(sensitivity to initial data: double precision resolves ~3 cycles)")

print("\nAll numerical checks passed." if failures == 0 else f"\n{failures} check(s) FAILED.")
sys.exit(failures != 0)
