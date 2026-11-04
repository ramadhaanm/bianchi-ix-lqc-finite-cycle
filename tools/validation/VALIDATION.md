# Validation

The logs in this directory were produced from the code and data in this repository.

- `data_audit.log`: `analysis/check_data.py` passed for 52,928 trajectory records and 558,654 saved valid events, including a 16-record smoke dataset. The audit checks source hashes, completeness, irreversible cutoffs and exits, constraint bounds, and the recorded event definitions.
- `statistics.log`: `analysis/test_reanalysis.py` passed all five statistical regression tests.
- `validity.log`: `julia --project=julia julia/test_validity.jl` passed 27 pathwise checks and 3 constraint-normalization checks.
- `cache_verifier_checks.log`: acceptance and rejection checks for `tools/validate_cache.py`. They show that the helper accepts matching data and rejects mismatched source metadata. They are not a cached/uncached simulation experiment.

The SHA256 hashes of all files are listed in `MANIFEST_SHA256.json`. Every raw data file also records the combined hash of the Julia model sources that produced it, and `analysis/reanalyze.py` refuses data whose recorded hash differs from the sources in `julia/`.

The cached Double64 integrator (`julia/run_cached.jl`) verifies, at every launch, that its three coefficient tables match the generic constructors byte for byte; the long-run logs record the cache hash. For a trajectory-level comparison of cached and uncached execution, generate N=2 datasets with `julia/reanalyze.jl` and `julia/run_cached.jl` in separate copies of the pinned project and compare them with `tools/validate_cache.py` (see `REPRODUCE.md`).

The equation checks (`derivations/`, `julia/test_baseline.jl` and `julia/test_A.jl`) are run with the commands in `REPRODUCE.md`.

These checks validate the numerical implementation and the finite-time statistics within their declared conditional populations. They do not establish that the effective Hamiltonian applies to every quantum state, or any generic late-time behavior.
