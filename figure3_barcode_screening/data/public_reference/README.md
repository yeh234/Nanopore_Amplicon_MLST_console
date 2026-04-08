# Public ligation reference data

This repository does not redistribute the published ligation-fidelity matrix used by the original exploratory script.

Place the following files in this directory before running `scripts/01_prepare_public_ligation_data.R`:

- `06_matrix.csv`
- `02_overhangs.csv`

You can override these paths with:

- `BARCODE_LIGATION_MATRIX_PATH`
- `BARCODE_OVERHANGS_PATH`

The new scripts expect `06_matrix.csv` to contain one overhang column plus one column per paired 4-bp sequence.

# Public Data 
https://github.com/potapovneb/CP-LigaseFidelity
