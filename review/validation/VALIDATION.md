# Validation supplied with the 27 October 2026 readiness revision

The manuscript uses the completed production matrix after the constraint-residual normalization fix. The readiness revision removes the uncontrolled band paragraph and updates documentation; it does not change the simulation equations, integration settings, production records, analysis definitions, or scientific tables/figures.

## Current evidence

- `data_audit.log`: `analysis/check_revised_data.py` passed for 52,928 trajectory records and 558,654 saved valid events, including the 16-record smoke dataset. The audit checks source hashes, completeness, irreversible cutoffs/exits, constraint bounds, and recorded event definitions.
- `statistics.log`: `analysis/test_reanalysis.py` passed all five statistical regression tests.
- `validity.log`: `julia --project=julia julia/test_validity.jl` passed 27 pathwise checks and 3 constraint-normalization checks.
- `pdf_build.log`: the final manuscript was built using `review/build_pdf.py` (pdfLaTeX, BibTeX, then two pdfLaTeX passes). References and citations resolve; no overfull boxes are reported. RevTeX float-placement and underfull-box warnings are checked by inspecting the resulting pages.
- `sources.json`: SHA256 hashes of the model/runner/cache, relevant analysis/tests, pinned environments, and the current production CSVs. The combined simulation hash uses the same file order as the runner and data loader.
- `six_cycle_normalization_comparison.log`: comparison of all 256 records per quantization against the local superseded matrix. All non-error columns, including event coordinates and outcomes, match exactly. Only normalized residual fields differ. Input SHA256 hashes identify the two versions. This compares the normalization versions, not cached versus uncached execution. Superseded input data are historical and are not part of the release.
- `cache_verifier_checks.log`: acceptance and rejection checks for the comparison helper. These establish that the helper accepts matching data and rejects mismatched source metadata; they are not a new cached/uncached simulation experiment.

These checks validate the stated numerical implementation and finite-time statistics within their declared conditional populations. They do not prove applicability of the effective Hamiltonian to every quantum state, generic late-time behavior, or journal acceptance.

## Historical cache evidence

`cache_historical_2026-10-26.log` is copied from the superseded production directory. Before the residual-normalization fix, it recorded exact agreement in every saved field of 16 cached/uncached controlled trajectories (two points, base, tighter-tolerance control, and two perturbations in each model). The original reference's first two points had been generated without the cache. The cache SHA256 in that report matches the current unchanged cache source.

This historical report is not presented as a fresh whole-trajectory comparison of the post-fix model sources. Every current cached launch separately checks all three immutable Double64 coefficient tables byte for byte against the generic constructors. The current long-run logs record the cache hash. Current N=64 production files were generated using the cache and are not labelled as uncached reference datasets.

For a fresh trajectory comparison, generate N=2 K and A datasets in separate pinned project copies using `julia/reanalyze.jl` and `julia/run_cached.jl`, respectively. Then compare the two result directories using the explicit directory arguments of `review/validate_cache.py`. Preserve the generation logs to establish the execution methods. No pilot dataset contributes to the paper's statistics.

## Earlier equation checks

The manuscript describes the symbolic, numerical, A-derivation, A-implementation, and baseline checks supplied in the project. Commands for running them are in `REPRODUCE_REVISED.md`. Their earlier results are not relabelled as fresh test runs in this readiness revision. The current evidence supplied above covers the data audit, statistical definitions, trajectory validity, residual-normalization regression, PDF build, and provenance of the cache evidence.
