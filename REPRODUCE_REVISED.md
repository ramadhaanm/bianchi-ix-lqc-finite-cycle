Reproducing the revised finite-cycle sensitivity study

The maintained manuscript is `paper/main.tex`. All revised LQC figures and tables use **schema 3** records in `results/revised/`, generated with pathwise length guards. Earlier files directly under `results/` are historical exploratory data and must not be mixed with the corrected LQC records. The retained classical reference files are `cl_basins_fig2_200_d64.csv` and `uncertainty_classical_d64.txt`.

Use Julia 1.13.2 with `julia/Project.toml` and `julia/Manifest.toml`, and Python 3.14 with `requirements.txt`. On Windows:

```powershell
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
```

On Linux/macOS use `.venv/bin/python`. The numerical scripts themselves are platform independent.

Verification commands:

```powershell
.venv/Scripts/python derivations/01_symbolic_checks.py
.venv/Scripts/python derivations/02_numerical_checks.py
.venv/Scripts/python derivations/03_A_quantization.py
julia --project=julia julia/test_baseline.jl
julia --project=julia julia/test_A.jl
julia --project=julia julia/test_validity.jl
.venv/Scripts/python analysis/test_reanalysis.py
```

`test_validity.jl` checks the audit counterexample at two tolerances, termination at the actual sub-Planckian crossing, rejection of initially invalid data, detection of a genuine constraint violation, and (since 2026-10-27) the constraint-normalization regression described below.

**Constraint-residual normalization fix (2026-10-27).** The relative residual in `constraint_error` (`julia/BianchiLQC.jl`) was normalized by the sum of |sin_i sin_j|. That sum collapses when two holonomy sines pass through zero together, which typically happens at bounces. Absolute residuals of about 1e-11 then appeared as relative errors of order one, and healthy trajectories were terminated as "constraint" exits: 144/1000 (K, L=20), 412/1000 (K, L=30) and 53/1000 (A, L=30).

The normalization is now the sum of sin_i^2 plus the positive potential scale (plus the |eta_i s_i| terms for the A quantization); the manuscript appendix gives the formula. The supplied false-alarm regression completes a formerly rejected cycle while keeping the corrected residual below 1e-6, and a separate check still detects a genuine constraint violation. The complete production matrix was rerun after this fix. Comparing the old and new six-cycle files gives identical event coordinates, outcomes, and all other non-error columns in all 256 records per model. The normalized residual columns change, so the complete serialized records are not identical. The six-cycle scientific statistics are unchanged. Evidence is supplied in `review/validation/six_cycle_normalization_comparison.log`.

The full production matrix completed on 27 October 2026 using `bash review/rerun_fixed.sh`: the same production ensembles as `run_reanalysis.ps1`, executed sequentially with all available threads. The local superseded outputs are preserved in `results/revised_v3_constraint_false_alarms_2026-10-26/`; that historical directory is excluded from the release and from the manuscript statistics. `test_reanalysis.py` checks zero-count confidence bounds, high/low storage of tiny separations, aggregation of phase gains, absence of later events after an exit, and withholding an interval for a single pair. The exit check uses the supplied small smoke dataset.

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

The default onset is 0.5 and the per-search budget is 300000 steps. N is the number of base points except for basin mode, where it is the grid side length. Phase mode also records the five-volume-e-fold turnaround event after each recollapse; other modes leave P slots empty by design. Each point has an independent fixed seed, so extending an ensemble does not alter existing initial conditions or perturbation directions.

For the complete production matrix, run the supplied PowerShell driver once for each quantization, concurrently on a machine with sufficient cores:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File review/run_reanalysis.ps1 K
powershell -NoProfile -ExecutionPolicy Bypass -File review/run_reanalysis.ps1 A
```

This changes execution policy only for the launched process. Each driver uses eight Julia threads. Individual commands in the table are an equivalent platform-independent route. Logs and resumable CSV records go under `results/revised/`. Completion is marked `DONE` in each log. A file with only some rows is not a completed ensemble.

Thread count controls scheduling of independent point/role jobs; it does not change their seeds or accuracy settings. The current post-fix long-run logs are `results/revised/long_20_K_64.log` and `results/revised/long_20_A_64.log`. Earlier continuation logs describe the superseded production run.

The full six-cycle rerun is part of the declared production matrix and can take several hours. Exploratory and performance pilots with N other than 64 are excluded from the manuscript and release. The previous mixed-precision outputs are not used. New six-cycle statistics require all 256 trajectory jobs per quantization, with consistent high/low storage and two perturbation sizes.

For the long runs, the driver uses `julia/run_cached.jl` with the same arguments as `reanalyze.jl`. This memoizes the immutable Double64 Vern9 coefficient tables in `julia/VernerDoubleCache.jl`; it does not change tolerances, event guards, or coefficient values. Each launch verifies all three tables byte for byte against the original generic constructors and records the cache source SHA256 in its log. Cached and uncached execution use the same CSV simulation-source hash for a given source version; the normalization fix changes that source hash. The uncached runner remains an equivalent, slower reproduction route.

A historical validation on 26 October, before the residual-normalization fix, compared 16 complete controlled trajectory records against uncached reference records and found exact agreement in every saved field. The unchanged cache's historical report is included as `review/validation/cache_historical_2026-10-26.log`; its date and scope are explained in `review/validation/VALIDATION.md`. It is distinct from the old/new normalization comparison. The current N=64 files were rerun with the cached runner and are not labelled as uncached references.

For an optional fresh cache comparison, use two isolated copies of the pinned project and generate N=2 K and A files with `reanalyze.jl` in one copy and `run_cached.jl` in the other. Keep the model sources identical. Compare the result directories with `.venv/Scripts/python review/validate_cache.py --reference-dir PATH_TO_UNCACHED_RESULTS --candidate-dir PATH_TO_CACHED_RESULTS --reference-n 2 --candidate-n 2`. The verifier checks metadata compatibility and every field of the first two points in both models. Pilot files never enter the scientific analysis and are excluded from the release.

Each raw file records schema, parameters, source hash, initial points, perturbations, cumulative cutoff masks, event coordinates, and stop reason. High/low components preserve Double64 separations. Resumption refuses mismatched code or settings. Do not erase metadata to resume altered calculations: preserve the previous files and start a distinct run. No trajectory continues after a physical or numerical exit.

The `u,v` columns are Float64 summaries of the base point. For exact extended-precision reproduction, regenerate the point from its integer ID using the runner's `MersenneTwister(110000 + id)` rule and its Double64 arithmetic; do not round-trip the displayed `u,v` summaries. Event logarithms, which determine the tiny separations, are saved with both high and low components. A confidence interval of a median is not estimated when only one pair is retained.

Once all production jobs are complete:

```powershell
.venv/Scripts/python analysis/reanalyze.py
.venv/Scripts/python analysis/check_revised_data.py
.venv/Scripts/python review/build_pdf.py
```

The analysis checks completeness before writing `results/revised/summary.json`, `results/revised/RESULTS.md`, revised figures, or tables/macros in `paper/generated/`. Manuscript numbers are generated from raw records. `--partial` instead writes `partial_summary.json` for monitoring and cannot update manuscript artifacts. The earlier `analysis/make_figures.py` must not rebuild the revised paper.

The complete-cycle statistic is the median of each pair's telescoping total gain. All phase statistics use the same resolved cohort for that cycle. The phase-marker comparison uses the first 200 point IDs and retains only the intersection of pairs passing controls at all three marker thresholds. Confidence intervals use 2000 point-level bootstrap resamples; label-fraction intervals are Wilson binomial intervals. Quality controls require signal-to-noise at least 100, relative vector discrepancy at most 10% across the two perturbation sizes, and separation at most 0.1. These conditional statistics do not establish the behavior of excluded or unresolved trajectories.

The release includes current audit and regression output under `review/validation/`, with source hashes and provenance in `VALIDATION.md` and `sources.json`. The audited production and smoke datasets contain 52,928 records and 558,654 saved valid events. Five statistical tests and 27 pathwise plus 3 residual-normalization checks passed. The exploratory persistent-band calculation is retained locally as historical work; it is excluded from the revised manuscript and supporting release because it does not implement the main analysis's numerical controls or per-pair gain statistics.

After compiling the final PDF, build the submission package:

```powershell
.venv/Scripts/python make_revised_release.py
```

This creates a separate versioned archive without overwriting `code_release/`. It excludes third-party reference PDFs. A public repository or DOI is not invented by local packaging. Licensing remains with the author; the archive has an explicit rights notice rather than granting redistribution rights on the author's behalf.
