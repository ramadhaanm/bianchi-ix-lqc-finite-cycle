# Finite-cycle sensitivity in effective loop quantum Bianchi IX cosmology

Code, raw data and manuscript source for "Finite-cycle sensitivity and phase-resolved amplification in
effective loop quantum Bianchi IX cosmology" (Muhammad Ramadhan, Department of Informatics Engineering,
Universitas Hasanuddin, Gowa, Indonesia).

- `julia/`: effective Bianchi IX integrators (K and A quantizations), event construction, tests
- `analysis/`: statistics, figures and table generation; `derivations/`: symbolic and numerical checks
- `results/revised/`: raw trajectory records and the summary used by the manuscript
- `paper/`: LaTeX source, generated tables and figures
- `REPRODUCE_REVISED.md`: full reproduction instructions (start here)

Code and data are released under the MIT License (see `LICENSE`).

The manuscript uses the IOP Publishing LaTeX class (`iopjournal.cls`), which is copyright IOP Publishing
and is not redistributed here. Download it from the IOP author template page to recompile `paper/main.tex`.
`paper/iopart-num-doi.bst` is a renamed copy of the `iopart-num.bst` BibTeX style with DOI printing added.
