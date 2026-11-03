#!/usr/bin/env bash
# Re-run of the v3 pipeline (same jobs as review/run_reanalysis.ps1) after the fix of the
# constraint-residual normalization (julia/BianchiLQC.jl, constraint_error), 2026-10-26.
# The superseded outputs are in results/revised_v3_constraint_false_alarms_2026-10-26/.
# Order: ordinary-precision ensembles for K and A first, extended-precision six-cycle runs last.
set -u
cd "$(dirname "$0")/.."
JULIA="${JULIA:-$LOCALAPPDATA/Microsoft/WindowsApps/julia.exe}"
mkdir -p results/revised
run() {   # run <runner> <job args...>
    local runner=$1; shift
    local log="results/revised/$(IFS=_; echo "$*").log"
    echo "$(date '+%F %T') starting $*"
    "$JULIA" --threads=auto --project=julia "$runner" "$@" > "$log" 2>&1
    echo "$(date '+%F %T') finished $* (exit $?)"
}
for Q in K A; do
    run julia/reanalyze.jl phase 20 $Q 1000
    run julia/reanalyze.jl phase 10 $Q 1000
    run julia/reanalyze.jl phase 30 $Q 1000
done
for Q in K A; do
    run julia/reanalyze.jl labels 20 $Q 400
    run julia/reanalyze.jl phase 20 $Q 200 0.3
    run julia/reanalyze.jl phase 20 $Q 200 0.7
done
run julia/reanalyze.jl basin 20 K 100
for Q in K A; do
    run julia/run_cached.jl long 20 $Q 64
done
echo "$(date '+%F %T') ALL RUNS COMPLETE"
