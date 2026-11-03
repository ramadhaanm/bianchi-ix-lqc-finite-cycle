"""Effective LQC ('K' quantization) and classical Bianchi I / II / IX dynamics.

Variables (verified against the canonical Hamiltonian in 01_symbolic_checks.py, section 5):
    x_i     = ln p_i                      (p1 = a2 a3, p2 = a1 a3, p3 = a1 a2)
    theta_i = gamma * mubar_i * K_i,      mubar_i = sqrt(Delta p_i / (p_j p_k))
    time T  with constant lapse N = gamma sqrt(Delta), i.e. proper time t = gamma sqrt(Delta) T

Equations of motion (LQC):
    x_i'     = cos th_i (sin th_j + sin th_k)
    th_i'    = -S + 1/2 sum_{l != i} x_l' (th_l - th_i) + F_i(x)
    S        = sum_pairs sin th_j sin th_k ,   Hamiltonian constraint:  S = W(x)
Classical GR is the same system with sin th -> th, cos th -> 1.

theta enters only through sin, cos and differences, so th_i -> th_i + n pi (all i together)
is an exact symmetry; `normalize` uses it to keep theta near 0 after each bounce, which
avoids the loss of precision that otherwise needs 128-bit arithmetic (cf. 2303.16747).

Units G = hbar = l_Pl = 1, gamma = 0.2375, Delta = 4 sqrt(3) pi gamma.
"""
import numpy as np

GAMMA = 0.2375
DELTA = 4 * np.sqrt(3) * np.pi * GAMMA
NWALL = {"I": (0, 0, 0), "II": (1, 0, 0), "IX": (1, 1, 1)}  # Bianchi class A structure constants


class Model:
    def __init__(self, kind="IX", lqc=True, gamma=GAMMA, delta=DELTA):
        self.kind, self.lqc = kind, lqc
        self.c = gamma**2 * delta / 4
        self.n = np.array(NWALL[kind], float)
        self.ix = kind == "IX"

    # --- potential: W (constraint) and F (force) ---------------------------
    def _walls(self, x):
        j, k = np.roll(x, -1), np.roll(x, -2)
        E = self.n * np.exp(j + k - 3 * x)  # E_i = e^{x_j + x_k - 3 x_i}   (a_i^4 walls)
        e = np.exp(x - j - k) if self.ix else np.zeros(3)  # e^{x_i - x_j - x_k}  (Bianchi IX cross terms)
        return E, e

    def W(self, x):
        E, e = self._walls(x)
        return self.c * (E.sum() - 2 * e.sum())

    def F(self, x):
        E, e = self._walls(x)
        return self.c * (-2 * E + 2 * np.roll(E, -1) + 2 * np.roll(E, -2) - 4 * e)

    # --- kinetic part --------------------------------------------------------
    def _trig(self, th):
        return (np.sin(th), np.cos(th)) if self.lqc else (th, np.ones(3))

    def S(self, th):
        s, _ = self._trig(th)
        return s[0] * s[1] + s[1] * s[2] + s[2] * s[0]

    def xdot(self, th):
        s, c = self._trig(th)
        return c * (s.sum() - s)

    def rhs(self, T, y):
        x, th = y[:3], y[3:]
        xd = self.xdot(th)
        S = self.S(th)
        # 1/2 sum_{l != i} x_l' (th_l - th_i) = 1/2 (sum_l x_l' th_l - th_i sum_l x_l')
        mix = 0.5 * (xd @ th - xd.sum() * th)
        thd = -S + mix + self.F(x)
        return np.concatenate([xd, thd])

    # --- diagnostics ---------------------------------------------------------
    def constraint_error(self, y):
        """|S - W| relative to the size of the individual terms (local version of E_r)."""
        x, th = y[:3], y[3:]
        s, _ = self._trig(th)
        E, e = self._walls(x)
        scale = abs(s[0] * s[1]) + abs(s[1] * s[2]) + abs(s[2] * s[0]) + self.c * (E.sum() + 2 * e.sum())
        return abs(self.S(th) - self.W(x)) / scale

    def dlnV(self, y):
        return 0.5 * self.xdot(y[3:]).sum()

    @staticmethod
    def lnV(y):
        return 0.5 * y[:3].sum()

    def kasner(self, y):
        """Geometric quasi-Kasner exponents k_i = H_i / sum_j H_j, with ln a_i = (x_j + x_k - x_i)/2."""
        xd = self.xdot(y[3:])
        H = 0.5 * (xd.sum() - 2 * xd)
        return H / H.sum()

    @staticmethod
    def normalize(y):
        """Exact symmetry theta_i -> theta_i - n pi (common n) bringing mean(theta) into [-pi/2, pi/2]."""
        y = y.copy()
        y[3:] -= np.pi * np.round(y[3:].mean() / np.pi)
        return y

    # --- initial data --------------------------------------------------------
    def initial_data(self, k, P, s):
        """Kasner exponents k at equal p_i = P with theta_i = s k_i (s < 0 contracting).

        theta_3 is then re-solved so that the constraint S = W holds to machine precision.
        Needs |s| << 1 (LQC negligible) and s^2 >> gamma^2 Delta / P (curvature negligible).
        """
        from scipy.optimize import brentq
        x = np.full(3, np.log(P))
        th = s * np.asarray(k, float)
        f = lambda t3: self.S(np.array([th[0], th[1], t3])) - self.W(x)
        lo, hi = th[2] - 0.5 * abs(th[2]), th[2] + 0.5 * abs(th[2])
        if f(lo) * f(hi) < 0:
            th[2] = brentq(f, lo, hi, xtol=1e-300, rtol=1e-15, maxiter=500)
        return np.concatenate([x, th])


class State:
    """Phase-space point with theta_i = n_i pi + d_i stored exactly (integers n_i, |d_i| <~ pi/2).

    Near a bounce one direction can sit at theta_i ~ -pi while the others are ~ 0; storing
    theta_i directly then loses the small physical offset to rounding. In this representation
    sin(theta_i) = (-1)^n_i sin(d_i) and theta_l - theta_i = (n_l - n_i) pi + (d_l - d_i)
    are evaluated without cancellation.
    """

    WRAP = np.pi / 2 + 0.1  # re-wrap d_i when |d_i| exceeds this (hysteresis avoids event chatter)

    def __init__(self, x, theta, T=0.0):
        theta = np.asarray(theta, float)
        self.x = np.asarray(x, float).copy()
        self.n = np.round(theta / np.pi).astype(int)
        self.d = theta - np.pi * self.n
        self.T = float(T)

    @property
    def y(self):
        """Lossy (x, theta) view, for diagnostics only."""
        return np.concatenate([self.x, self.n * np.pi + self.d])

    def rewrap(self):
        k = np.round(self.d / np.pi).astype(int)
        self.d -= np.pi * k
        self.n += k
        self.n -= int(np.round(self.n.mean()))  # exact common shift theta -> theta - m pi


class WrappedModel(Model):
    """Same dynamics as Model, evaluated on the exact (x, n, d) representation."""

    def _sc(self, d, n):
        sign = 1.0 - 2.0 * (n % 2)
        if self.lqc:
            return sign * np.sin(d), sign * np.cos(d)
        return n * np.pi + d, np.ones(3)

    def rhs_nd(self, T, z, n):
        x, d = z[:3], z[3:]
        s, c = self._sc(d, n)
        xd = c * (s.sum() - s)
        S = s[0] * s[1] + s[1] * s[2] + s[2] * s[0]
        dth = (n[None, :] - n[:, None]) * np.pi + (d[None, :] - d[:, None])  # dth[i, l] = th_l - th_i
        mix = 0.5 * (dth @ xd)  # 1/2 sum_l x_l' (th_l - th_i); the l = i term is zero
        return np.concatenate([xd, -S + mix + self.F(x)])

    def constraint_error_nd(self, z, n):
        x, d = z[:3], z[3:]
        s, _ = self._sc(d, n)
        E, e = self._walls(x)
        S = s[0] * s[1] + s[1] * s[2] + s[2] * s[0]
        scale = abs(s[0] * s[1]) + abs(s[1] * s[2]) + abs(s[2] * s[0]) + self.c * (E.sum() + 2 * e.sum())
        return abs(S - self.W(x)) / scale

    def dlnV_nd(self, z, n):
        s, c = self._sc(z[3:], n)
        return 0.5 * (c * (s.sum() - s)).sum()

    def evolve(self, st, stop, direction=0, T_max=1e30, rtol=1e-12):
        """Advance State `st` in place until stop(T, z, n) crosses zero in `direction`.

        Returns the worst relative constraint error seen along the way.
        """
        from scipy.integrate import solve_ivp
        atol = np.array([1e-14] * 3 + [1e-26] * 3)
        worst, T_end = 0.0, st.T + T_max
        while True:
            n = st.n.copy()
            ev_stop = lambda T, z: stop(T, z, n)
            ev_stop.terminal, ev_stop.direction = True, direction
            wraps = []
            for i in range(3):
                w = lambda T, z, i=i: z[3 + i] ** 2 - State.WRAP**2
                w.terminal, w.direction = True, +1
                wraps.append(w)
            sol = solve_ivp(lambda T, z: self.rhs_nd(T, z, n), (st.T, T_end), np.concatenate([st.x, st.d]),
                            method="DOP853", rtol=rtol, atol=atol, events=[ev_stop] + wraps)
            worst = max(worst, max(self.constraint_error_nd(sol.y[:, j], n) for j in range(sol.y.shape[1])))
            st.x, st.d, st.T = sol.y[:3, -1].copy(), sol.y[3:, -1].copy(), sol.t[-1]
            if sol.status == 0:
                raise RuntimeError("stop condition not reached before T_max")
            if sol.status < 0:
                raise RuntimeError(sol.message)
            if len(sol.t_events[0]):
                return worst
            st.rewrap()


def kasner_from_u(u):
    """BKL parametrization: k = (-u, 1+u, u(1+u)) / (1+u+u^2)."""
    d = 1 + u + u * u
    return np.array([-u, 1 + u, u * (1 + u)]) / d
