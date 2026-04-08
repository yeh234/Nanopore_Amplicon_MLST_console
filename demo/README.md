# Demo FASTQ

`Samlst_SUP.demo.fastq.gz` is intended for GitHub-friendly testing.

Recommended build basis:

- First rerun Step 1 once in the current environment against the full raw input.
- Then keep only reads from the original `Samlst_SUP.fastq.gz` that that
  current-environment `read_calls` table assigned to a concrete sample.
- This preserves the per-sample FASTQ content used by downstream assembly/MLST,
  as much as possible while removing reads that mainly inflate the raw upload
  size.

Build command:

```bash
python3 code/07_make_demo_fastq_from_read_calls.py \
  --input-fastq=fastq/Samlst_SUP.fastq.gz \
  --read-calls=output/_app_smoke_step1/read_calls/Samlst_SUP.tsv \
  --output-fastq=demo/Samlst_SUP.demo.fastq.gz
```

Validation in this workspace:

- The demo FASTQ is about 30 MB instead of the original 266 MB.
- The final MLST/ST output matches the current reference for all samples:
  `MLST strain type`, `Genotype Result`, and the 7 allele columns are all
  identical in `demo_publication.tsv`.
- `Reads number` is not bitwise identical, because removing non-demo reads still
  causes a few sample counts to shift during Step 1 reruns in the current
  environment.

If you need exact read counts too:

- A small raw FASTQ is not sufficient in this repo's current Step 1 behavior.
- Use Git LFS, split a larger filtered FASTQ for release, or publish the
  Step 1 `merged_by_sample/` outputs as the exact downstream test fixture.
