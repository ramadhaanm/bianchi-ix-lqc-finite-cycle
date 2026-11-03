# Derivations, Part 1: Effective LQC Bianchi IX, from the Hamiltonian to a verified integrator

*2026-10-24. Every equation below is checked by computer: [01_symbolic_checks.py](01_symbolic_checks.py) (16 symbolic checks) and [02_numerical_checks.py](02_numerical_checks.py) (5 numerical checks). Run them with `.venv/Scripts/python derivations/<file>`.*

Conventions follow Blackmore & Wilson-Ewing, CQG 40 (2023) [arXiv:2303.16747]: the "K" loop quantization, with G = ħ = ℓ_Pl = 1, γ = 0.2375 and Δ = 4√3 πγ.

## 1. Geometry

Diagonal Bianchi IX: $ds^2 = -N^2dt^2 + \sum_i a_i^2\,\sigma_i^2$, with $\sigma_i$ the left-invariant 1-forms on $SU(2)\simeq S^3$.
Triad variables: $p_1 = a_2a_3,\ p_2 = a_1a_3,\ p_3 = a_1a_2$, so that $q = \det h = p_1p_2p_3$.

**Spatial curvature.** It was computed from scratch: Christoffel symbols of the metric in Euler coordinates, compared at six random exact points.
$$ {}^3R = -\frac{\sum_i a_i^4 - 2\sum_{i<j}a_i^2a_j^2}{2a_1^2a_2^2a_3^2}, \qquad {}^3R = -\frac{16\pi G\,U_{IX}}{q}, $$
$$ U_{IX} = \frac{1}{32\pi G}\Big[\frac{p_2^2p_3^2}{p_1^2}+\frac{p_1^2p_3^2}{p_2^2}+\frac{p_1^2p_2^2}{p_3^2} - 2(p_1^2+p_2^2+p_3^2)\Big]. $$
In the isotropic limit, ${}^3R = 3/(2a^2)$: a round $S^3$ of radius $2a$.

## 2. Classical Hamiltonian constraint

From $\frac{N\sqrt q}{16\pi G}\big[(K_{ab}K^{ab}-K^2) - {}^3R\big]$, with $\{K_i,p_j\} = 8\pi G\,\delta_{ij}$:
$$ C = \frac{N}{\sqrt q}\Big[-\frac{1}{8\pi G}\big(p_1p_2K_1K_2 + p_2p_3K_2K_3 + p_3p_1K_3K_1\big) + U\Big]. $$
Verified:
- both terms match the ADM expression;
- with N = 1, $K_i = \dot a_i$ and $\dot p_1 = \frac{d}{dt}(a_2a_3)$;
- the Kasner solution solves the Bianchi I case;
- $k_i = p_iK_i/\sum_j p_jK_j$.

## 3. Effective LQC ("K" quantization), lapse $N=\sqrt q$

$$ C_{\rm eff} = -\frac{q}{8\pi G\gamma^2\Delta}\sum_{\rm cyc}\sin\theta_j\sin\theta_k + U, \qquad \theta_i = \gamma\bar\mu_iK_i,\quad \bar\mu_i = \sqrt{\frac{\Delta\,p_i}{p_jp_k}}. $$
Verified:
- $\Delta\to0$ recovers Sec. 2, because $\bar\mu_j\bar\mu_k = \Delta/p_i$;
- the post-bounce "shift" used by Blackmore & Wilson-Ewing, $\theta_i\to\theta_i+\pi$ for all $i$ together, is an **exact canonical symmetry**: $\pi/(\gamma\bar\mu_i) = \partial_{p_i}F$ with $F = 2\pi\sqrt q/(\gamma\sqrt\Delta)$.

## 4. Reformulation used for the numerics (new)

Take $x_i = \ln p_i$ and $\theta_i$ as variables, with the **constant lapse $N = \gamma\sqrt\Delta$**, so proper time is $t = \gamma\sqrt\Delta\,T$:
$$ x_i' = \cos\theta_i(\sin\theta_j + \sin\theta_k), \qquad
   \theta_i' = -S + \tfrac12\sum_{l\neq i}x_l'(\theta_l-\theta_i) + F_i(x), \qquad S = \sum_{\rm cyc}\sin\theta_j\sin\theta_k, $$
$$ \text{constraint: } S = W(x) = \tfrac{\gamma^2\Delta}{4}\Big[\sum_i E_i - 2\sum_i e^{x_i-x_j-x_k}\Big],\quad E_i = e^{x_j+x_k-3x_i},\quad
   F_i = \tfrac{\gamma^2\Delta}{4}\big[-2E_i+2E_j+2E_k-4e^{x_i-x_j-x_k}\big]. $$
*Derivation sketch:* $\theta_i' = \theta_i\,(\ln\bar\mu_i)' + \gamma\bar\mu_iK_i'$, with $(\ln\bar\mu_i)' = \tfrac12(x_i'-x_j'-x_k')$. The explicit terms in $\theta$ combine so that **only differences $\theta_l-\theta_i$ survive**; the $l=i$ term cancels.

Verified symbolically against Hamilton's equations for $C_{\rm eff}$, and numerically at 20 random phase-space points (agreement to 2×10⁻¹²).

Advantages over the canonical form:
- All rates are O(1) instead of ~10⁷.
- $\theta$ enters only through sin, cos and differences, which is what makes the shift symmetry manifest.
- Classical GR is the same system with $\sin\theta\to\theta$, $\cos\theta\to1$.
- Kasner exponents: $k_i = H_i/\sum H$, with $H_i = \tfrac12(x_j'+x_k'-x_i')$. They are geometric and independent of the lapse.

## 5. Exact angle representation (new)

Near a bounce, one direction can sit at $\theta_i\approx-\pi$ while the others are still ≈ 0. A common shift cannot help, and the physical offset from $-\pi$ is lost to rounding. This loss is why Blackmore & Wilson-Ewing needed 128-bit arithmetic.

**Fix:** store $\theta_i = n_i\pi + \delta_i$, with an integer $n_i$ and $|\delta_i|\lesssim\pi/2$. Then use
- $\sin\theta_i = (-1)^{n_i}\sin\delta_i$, and
- $\theta_l-\theta_i = (n_l-n_i)\pi + (\delta_l-\delta_i)$,

both of which are computed without cancellation. $\delta_i$ is re-wrapped whenever $|\delta_i| > \pi/2 + 0.1$ (this is a change of representation, not an approximation). The common integer is removed using the exact symmetry. Verified to give the same right-hand side to 6×10⁻¹⁴.

## 6. Numerical baseline: published results reproduced

| Test | Result |
|---|---|
| A. LQC Bianchi I bounce, $k_i\to\frac23-k_i$, initial data of Fig. 1 of 2303.16747 | $k_{\rm out} = (-0.287831, 0.434273, 0.853558)$ vs. paper $(-0.287832, 0.434273, 0.853558)$ |
| B. Classical Bianchi II wall bounce (BKL rule) | agrees to 10⁻¹⁰: $(-0.31578947, 0.52631579, 0.78947368)$ |
| C. LQC Bianchi IX, 6 bounce/recollapse cycles | runs cleanly; bounce volume grows from cycle to cycle (as also seen in 2501.14891); first 4 extrema converged to < 10⁻⁸ |

## 7. First finding: sensitivity to initial data

Compare two runs that differ only in integration tolerance (rtol 10⁻¹² vs 10⁻¹³). Their difference in $\ln V$ at successive extrema is:

| Extremum | B1 | R1 | B2 | R2 | B3 | R3 | B4 | R4 |
|---|---|---|---|---|---|---|---|---|
| $\lvert\Delta\ln V\rvert$ | 0 | 6e-13 | 2e-10 | 3e-9 | 5e-7 | 3e-5 | 5e-2 | 8e-2 |

(B = bounce, R = recollapse.) The difference grows about **×70 per volume extremum** (geometric mean), i.e. roughly 3–4 decades per cycle.

- **What it means for the numerics:** with double precision, an individual trajectory is reliable for about **3 cycles**. Each extra digit of working precision buys less than one more cycle. Basin outcomes beyond cycle 3 need extended precision (Julia `Double64`/`Float128`, or C++ quad precision); ensemble statistics (RQ2) can rely on shadowing but must be tested for convergence.
- **What it means for the physics (preliminary):** the volumes at bounces and recollapses are **geometric scalars**, so the growth of their separation does not depend on the time coordinate. It is a candidate invariant chaos indicator for RQ1, alongside fractal basin boundaries. *Caveat:* numerical round-off is not a controlled perturbation. The proper test perturbs initial data along the constraint surface (next step).

## 7b. Julia port and extended precision (2026-10-25)

[julia/BianchiLQC.jl](../julia/BianchiLQC.jl) implements the same equations as `WrappedModel`, written generically over the number type (Float64, `Double64`, `Float128`, `BigFloat`). It uses the Vern9 integrator. Initial data are now in closed form: the constraint is linear in $\sin\theta_3$, so $\theta_3 = \arcsin\big((W - s_1s_2)/(s_1+s_2)\big)$. Run the checks with `julia --project=julia julia/test_baseline.jl`; all pass.

- **Float64:** about 50× faster than the Python version (0.3 s for 3 cycles); the extrema match Python.
- **Double64:** about 32 digits; one cycle takes ~5–10 s at reltol 10⁻²⁴–10⁻²⁶. Difference in $\ln V$ between reltol 10⁻²⁴ and 10⁻²⁶ at successive extrema:

| Extremum | B1 | R1 | B2 | R2 | B3 | R3 | B4 | R4 | B5 | R5 | B6 | R6 | B7 | R7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $\lvert\Delta\ln V\rvert$ | 5e-27 | 4e-26 | 1e-23 | 3e-23 | 5e-21 | 3e-19 | 5e-16 | 3e-16 | 8e-13 | 3e-12 | 2e-9 | 4e-8 | 6e-6 | 3e-3 |

  Individual trajectories are resolved through **6 cycles** (Float64: 3). The growth of about 2–3 decades per extremum matches the Float64 estimate, which supports reading it as sensitivity of the dynamics rather than a numerical artifact.

**Implementation caveats found and fixed:**
- In Julia, `k * pi` is evaluated in Float64 even inside Double64 code. The working-precision π is now stored in the model.
- `kasner(m, u)` reads the current offsets `m.n`, so it must be evaluated for the state those offsets belong to.
- The rule $k_i\to\frac23-k_i$ holds only asymptotically. At $\theta\sim10^{-3}$ it is satisfied to ~7×10⁻⁸ (the $O(\theta^2)$ corrections). The tests now say so explicitly; the earlier Python check had silently used a 10⁻⁵ relative tolerance.

## 8. Next derivation steps

1. **Invariant outcome for the basins.** Candidate: at the n-th recollapse (a volume maximum, a geometric event), which $a_i$ is largest, or the sign pattern of $k_i$ just after it. Then show the label is invariant under time reparametrization and under permutations of the axes.
2. **Classical baseline.** Reproduce the Cornish–Levin fractal basins (gr-qc/9605029) with the classical system of Sec. 4.
3. **Controlled perturbations.** Measure the growth of $\delta\ln V_{\rm extremum}$ under perturbations on the constraint surface; compare with the ×70 estimate.
4. **Speed.** Port the right-hand side to compiled code (numba or Julia) before building basin grids (~10⁵–10⁶ trajectories).
