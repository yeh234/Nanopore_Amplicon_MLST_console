# OT-2 Protocols

This directory contains the curated Opentrons OT-2 materials that support the automated liquid-handling portion of the study workflow.

Scope note:

- the OT-2 notebooks were originally developed against larger 24-barcode / 108-combination design libraries during assay development
- for this manuscript, the intended study scope is the 10 selected barcode combinations used in the proof-of-concept validation experiment
- the public-facing repository now uses normalized `paper_10_*` layout filenames for the manuscript-specific subset
- OT-2 execution cells still reference `.csv` files copied onto the robot server, whereas this repository preserves the corresponding `.xlsx` layout sources
- historical ST45 in the original design library corresponds to the current pubMLST profile designation ST508 for the first paper sample

Relation to manuscript sections:

- Section 2.2 `Barcode-bearing primer design and PCR amplification`
  - `step01_combination.ipynb`: barcode-bearing primer combination, DNA-to-position mapping, and PCR-oriented plate setup
  - `step02_MagneticPurification_optimazation.ipynb`: post-PCR magnetic bead cleanup optimization
- Section 2.3 `Golden Gate-based concatenated-amplicon construction`
  - `step03_GoldenGateAssemble.ipynb`: layout-driven Golden Gate assembly setup for ordered concatenated amplicons
  - `step04_MagneticPurification_largeSize.ipynb`: post-assembly magnetic cleanup for larger concatenated products
- Section 2.4 `Automated liquid-handling workflow`
  - all four notebooks use the Opentrons API with OT-2 single-channel pipettes, NEST 96-well PCR plates, 24-tube racks, tip racks, and magnetic-module workflows
- Section 2.5 `Measures to minimize cross-contamination during automated handling`
  - the notebooks repeatedly use fresh-tip transfer patterns and fixed deck positions for magnetic bead cleanup
  - additional contamination-control practices described in the manuscript are laboratory SOP items and are not exhaustively encoded in notebook text

Included notebooks:

- `step01_combination.ipynb`
- `step02_MagneticPurification_optimazation.ipynb`
- `step03_GoldenGateAssemble.ipynb`
- `step04_MagneticPurification_largeSize.ipynb`

Paper-facing layout files:

- `datalayout/paper_10_barcode_combination.xlsx`: manuscript-facing 10-barcode subset
- `datalayout/paper_10_sample_order.xlsx`: manuscript-facing 10-sample order and ST mapping
- `datalayout/paper_10_tube_position.xlsx`: normalized OT-2 tube-position layout for the paper subset
- `datalayout/paper_10_endproduct_layout.xlsx`: normalized end-product layout used in the paper-facing OT-2 materials
- `datalayout/paper_10_endproduct_layout_v1.xlsx`: alternate normalized end-product layout retained with manuscript-facing naming
- `datalayout/paper_10_new_endproduct_layout_v1.xlsx`: normalized later end-product layout variant retained with manuscript-facing naming
- `datalayout/paper_10_38_input_layout.xlsx`: normalized 38-position input layout
- `datalayout/paper_10_38_output_layout.xlsx`: normalized 38-position output layout
- `datalayout/paper_10_96plate_layout.xlsx`: normalized 96-well plate layout
- `datalayout/paper_10_96plate_endproduct.xlsx`: normalized 96-well end-product layout

Each paper-facing workbook also has a same-named `.csv` export in `datalayout/` for OT-2 upload to `/root/script/`. These `.csv` files were exported from the first worksheet of the corresponding `paper_10_*.xlsx` file and are intended as the direct notebook input files.

Uploading CSV files to the OT-2 server:

From the repository root on your workstation, set `OT2_HOST` to the robot's
address and copy the paper-facing CSV files to its Jupyter environment.

```bash
OT2_HOST=your-ot2-host
scp OT-2_protocols/datalayout/paper_10_*.csv "root@${OT2_HOST}:/root/script/"
ssh "root@${OT2_HOST}" 'ls -1 /root/script/paper_10_*.csv'
```

For a specific SSH identity, set `OT2_SSH_KEY` to its local file path and add
`-i "$OT2_SSH_KEY"` to both commands. The commented upload examples in the
notebooks use these two workstation settings. The notebooks read the
`paper_10_*.csv` files from `/root/script/` on the robot.

Historical archive files:

- `datalayout/barcode108_combination.xlsx`: historical full design-library workbook from which the paper-relevant subset was selected
- `datalayout/DNAtoPosition.xlsx`
- `datalayout/16tube_position.xlsx`
- `datalayout/16_endproductLayout_set.xlsx`
- `datalayout/16_endproductLayout_setv1.xlsx`
- `datalayout/16new_endproductLayout_v1.xlsx`
- `datalayout/38_Input_Layout.xlsx`
- `datalayout/38_Output_Layout.xlsx`
- `datalayout/96plate_layout.xlsx`
- `datalayout/96plate_endproduct.xlsx`

Repository notes:

- the notebooks preserve historical run logs from OT-2 development
- the notebooks use normalized `paper_10_*` names for the manuscript-facing subset, while keeping historical development workbooks as archive material
- OT-2 execution still expects exported `.csv` files under `/root/script`; the repository stores the corresponding `.xlsx` source layouts
- unrelated layout, test, and reference tables have been removed from the public-facing OT-2 folder
- the manuscript scope is narrower than the historical OT-2 development scope: only the 10 selected barcode combinations belong to the paper-facing workflow

Paper-only release recommendation:

- cite the `paper_10_*` layout files as the manuscript-facing OT-2 materials
- keep `barcode108_combination.xlsx` and the non-`paper_10_*` layout workbooks only as archived development provenance, or move them out of the paper-facing OT-2 folder

If you later need a stricter paper-only release, keep only the notebooks and layout files explicitly cited in the final manuscript and archive the remaining workflow variants separately.
