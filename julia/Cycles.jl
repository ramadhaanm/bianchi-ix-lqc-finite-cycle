"""
Phase-based definition of the geometric events of one bounce-recollapse cycle, shared by the
analysis scripts. Needs BianchiLQC to be loaded first (include BianchiLQC.jl, then this file).

Events of cycle n (independent of the integration time parametrization):
  On  : LQC onset, max_i |sin theta_i| first reaches 0.5 while contracting (approached in stages:
        each evolve call restarts the clock, and after a long contraction T ~ 1e8 would limit the
        event-time resolution).
  B   : the first volume minimum after On.
  Off : end of the quantum phase. max |sin theta| falls below 0.5 *and then reaches the classical
        level 0.05 without returning to 0.5*. If it returns above 0.5 first (multiple bounce), the
        search continues; the number of such re-entries is counted.
  R   : the recollapse, the first volume maximum in the classical regime (max |sin theta| < 0.01).
        Volume maxima at higher |sin theta| are part of a multiple bounce; they are skipped and
        counted.
Validity is monitored along every integration segment, including offset look-ahead and
multiple bounces. The bounce-only minimum is retained for diagnostics, not selection.
"""

const ONSET = 0.5
const CLASSICAL_END = 0.05
const CLASSICAL_MAX = 0.01

"""Contracting classical state -> On (staged)."""
function to_onset!(m, u; onset = ONSET, kw...)
    worst = 0.0
    for thr in (1e-3, 1e-2, 0.1, onset)
        lqc_strength(m, u) < thr || continue
        u, e = evolve_until!(m, u, (m, v) -> lqc_strength(m, v) - thr, +1; kw...)
        worst = max(worst, Float64(e))
    end
    return u, worst
end

"""From B: the end of the quantum phase. Returns (u at Off, re-entries, worst constraint error)."""
function to_offset!(m, u; onset = ONSET, kw...)
    reentries = 0; worst = 0.0
    while true
        u, e = evolve_until!(m, u, (m, v) -> lqc_strength(m, v) - onset, -1; kw...)   # candidate Off
        worst = max(worst, Float64(e))
        uoff = copy(u); noff = copy(m.n)
        # positive between CLASSICAL_END and ONSET: stops at whichever level is crossed first
        g(m, v) = (lqc_strength(m, v) - CLASSICAL_END) * (onset - lqc_strength(m, v))
        u, e = evolve_until!(m, u, g, -1; kw...)
        worst = max(worst, Float64(e))
        if lqc_strength(m, u) < (CLASSICAL_END + onset) / 2
            m.n .= noff
            return uoff, reentries, worst          # confirmed: the quantum phase ended at uoff
        end
        reentries += 1                              # back above 0.5: still in the quantum phase
    end
end

"""From Off: the first classical volume maximum. Returns (u at R, skipped quantum maxima, worst error)."""
function to_recollapse!(m, u; kw...)
    skipped = 0; worst = 0.0
    while true
        u, e = evolve!(m, u, -1; kw...)
        worst = max(worst, Float64(e))
        lqc_strength(m, u) < CLASSICAL_MAX && return u, skipped, worst
        skipped += 1
        u, e = evolve!(m, u, +1; kw...)             # pass the following minimum
        worst = max(worst, Float64(e))
    end
end

"""
    run_cycle!(m, u; record, kw...) -> (u, info)

One cycle On -> B -> Off -> R starting from a contracting classical state. `record(name, u)` is
called at each event. info: (reentries, skipped_max, lna_min_B, worst constraint error).
"""
function run_cycle!(m, u; record = (name, u) -> nothing,
                    monitor = ValidityMonitor(u), onset = ONSET, kw...)
    u, e1 = to_onset!(m, u; onset=onset, monitor=monitor, kw...); record(:On, u)
    u, e2 = evolve!(m, u, +1; monitor=monitor, kw...); record(:B, u)
    lna_min = minimum(loga(u))
    u, re, e3 = to_offset!(m, u; onset=onset, monitor=monitor, kw...); record(:Off, u)
    u, sk, e4 = to_recollapse!(m, u; monitor=monitor, kw...); record(:R, u)
    return u, (reentries = re, skipped_max = sk, lna_min_B = Float64(lna_min),
        lna_min_path = Float64(monitor.minimum), cutmask = monitor.cutmask,
        err = max(e1, e2, e3, e4))
end
