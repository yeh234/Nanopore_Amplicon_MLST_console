# PubMLST reference layout

This directory documents the expected local layout for the *Staphylococcus
aureus* MLST allele and profile files used by the workflow.

Expected local files:

- `arcC.fas`
- `aroE.fas`
- `glpF.fas`
- `gmk.fas`
- `pta.fas`
- `tpi.fas`
- `yqiL.fas`
- `saureus_Profile.txt`

Workflow default:

- Step 3 uses `db/saureus_Profile.txt`

Users can provide a local snapshot here or regenerate the files from the
official PubMLST sources.

Suggested reconstruction steps:

- go to the official PubMLST MLST download page
- in the `Staphylococcus aureus` row, download:
  - `profiles`
  - `arcC`
  - `aroE`
  - `glpF`
  - `gmk`
  - `pta`
  - `tpi`
  - `yqiL`
- rename or save them locally as:
  - `saureus_Profile.txt`
  - `arcC.fas`
  - `aroE.fas`
  - `glpF.fas`
  - `gmk.fas`
  - `pta.fas`
  - `tpi.fas`
  - `yqiL.fas`

Known direct profile endpoint:

- `https://rest.pubmlst.org/db/pubmlst_saureus_seqdef/schemes/1/profiles_csv`

Official references:

- API overview: <https://pubmlst.org/api>
- API base URL: <https://rest.pubmlst.org/>
- MLST download page: <https://dev.pubmlst.org/data>
- terms and conditions: <https://pubmlst.org/terms-conditions>
