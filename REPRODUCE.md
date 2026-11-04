# Reproducing the results

All figures and tables of the manuscript are generated from the raw trajectory records in `results/revised/`. Use Julia 1.13.2 with `julia/Project.toml` and `julia/Manifest.toml`, and Python 3.14 with `requirements.txt`. On Windows:

```powershell
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
```

On Linux/macOS use `.venv/bin/python`. The numerical scripts are platform independent.

## Verification

```powershell
.venv/Scripts/python derivations/01_symbolic_checks.py
.venv/Scripts/python derivations/02_numerical_checks.py
.venv/Scripts/python derivations/03_A_quantization.py
julia --project=julia julia/test_baseline.jl
julia --project=julia julia/test_A.jl
julia --project=julia julia/test_validity.jl
.venv/Scripts/python analysis/test_reanalysis.py
```

`test_validity.jl` checks termination at the actual sub-Planckian crossing, rejection of initially invalid data, detection of a genuine constraint violation, and the constraint-residual normalization described below. `test_reanalysis.py` checks zero-count confidence bounds, high/low storage of tiny separations, aggregation of phase gains, absence of later events after an exit, and withholding of an interval for a single pair.

The normalized constraint residual is the absolute residual divided by the sum of `sin^2(theta_i)` plus the positive potential scale (for the A quantization, plus the `|eta_i s_i|` terms). This normalization stays finite when two holonomy sines vanish together, as happens at bounces. A residual of 1e-3 terminates a run; it is not an accuracy certificate for trajectory separations.

## Production runs

Production calculations use `julia/reanalyze.jl`:

```text
julia --threads=8 --project=julia julia/reanalyze.jl MODE L QUANT N [ONSET] [BUDGET]
```

| Mode | Initial scale L | Quantization | N | Numerical controls |
|---|---|---|---:|---|
| phase | 10, 20, 30 | K and A | 1000 each | Three cycles, eps 1e-9 and 1e-10, reltol 1e-12 and control 1e-13 |
| labels | 20 | K and A | 400 each | Three cycles, 17 radii from 1e-2 to 1e-10, control 1e-13 |
| long | 20 | K and A | 64 each | Six cycles, Double64, eps 1e-24 and 1e-25, reltol 1e-28 and control 1e-29 |
| phase | 20 | K and A | 200 each | Repeat at ONSET 0.3 and 0.7 |
| basin | 20 | K | 100 | 100 by 100 grid, two cycles |

The default onset is 0.5 and the per-search budget is 300000 steps. N is the number of base points, except in basin mode where it is the grid side length. Phase mode also records the five-volume-e-fold turnaround event after each recollapse; other modes leave these slots empty by design. Each point has an independent fixed seed, so extending an ensemble does not alter existing initial conditions or perturbation directions.

To run the complete production matrix, either run the PowerShell driver once for each quantization (concurrently, on a machine with enough cores):

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_reanalysis.ps1 K
powershell -NoProfile -ExecutionPolicy Bypass -File tools/run_reanalysis.ps1 A
```

or run everything sequentially with `bash tools/run_all.sh`. Logs and resumable CSV records are written to `results/revised/`; completion is marked `DONE` in each log, and a file with only some rows is not a completed ensemble. The six-cycle runs take several hours. Thread count controls the scheduling of independent jobs; it does not change their seeds or accuracy settings.

For the long runs the drivers use `julia/run_cached.jl` with the same arguments as `reanalyze.jl`. It memoizes the immutable Double64 Vern9 coefficient tables in `julia/VernerDoubleCache.jl` and does not change tolerances, event guards or coefficient values. Each launch verifies all three tables byte for byte against the generic constructors and records the cache source SHA256 in its log. To compare cached and uncached execution on a small sample, generate N=2 K and A files with `reanalyze.jl` in one copy of the pinned project and with `run_cached.jl` in another, then run `python tools/validate_cache.py --reference-dir PATH_TO_UNCACHED_RESULTS --candidate-dir PATH_TO_CACHED_RESULTS --reference-n 2 --candidate-n 2`.

## Raw records

Each raw file records the schema, parameters, source hash, initial points, perturbations, cumulative cutoff masks, event coordinates and stop reason. Resumption refuses mismatched code or settings. No trajectory continues after a physical or numerical exit. The `u,v` columns are Float64 summaries of the base point. For exact extended-precision reproduction, regenerate the point from its integer ID using the runner's `MersenneTwister(110000 + id)` rule and its Double64 arithmetic. Event logarithms, which determine the tiny separations, are saved with both high and low components.

## Analysis

Once all production jobs are complete (the analysis writes tables and figures to `paper/generated/` and `paper/figs/`, so create `paper/` first):

```powershell
mkdir paper
.venv/Scripts/python analysis/reanalyze.py
.venv/Scripts/python analysis/check_data.py
```

The analysis checks completeness before writing `results/revised/summary.json`, `results/revised/RESULTS.md`, the figures and the manuscript tables and macros. `--partial` instead writes `partial_summary.json` for monitoring and cannot update manuscript artifacts.

The complete-cycle statistic is the median of each pair's telescoping total gain, and all phase statistics use the same resolved cohort for that cycle. The phase-marker comparison uses the first 200 point IDs and retains the pairs that pass the controls at all three marker thresholds. Confidence intervals use 2000 point-level bootstrap resamples; label-fraction intervals are Wilson binomial intervals. Quality controls require a signal-to-noise ratio of at least 100, a relative vector discrepancy of at most 10% across the two perturbation sizes, and a separation of at most 0.1. These are conditional statistics and do not describe excluded or unresolved trajectories.

## Validation evidence

`tools/validation/` holds the output of the data audit (52,928 trajectory records and 558,654 saved valid events), the statistical regression tests, the pathwise-validity tests (27 pathwise and 3 residual-normalization checks) and the cache-verifier checks. See `tools/validation/VALIDATION.md`.

`make_release.py` writes `MANIFEST_SHA256.json` and builds a versioned archive of the code, data and scripts.
