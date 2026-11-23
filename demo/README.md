# Demo FASTQ

`Samlst_SUP.demo.fastq.gz` contains 11,950 reads from the full study FASTQ for
a lightweight test of Step 1. See Option B in the [repository README](../README.md).

Recommended build basis:

- First rerun Step 1 once in the current environment against the full raw input.
- Then keep only reads from the original `Samlst_SUP.fastq.gz` that the
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

In the R1 verification, Step 1 assigned 11,939 of these demo reads. A downstream
rerun for samples 1 and 4 recovered their expected seven-locus profiles. For
validation against fixed downstream inputs and their reference outputs, use
Option C in the repository README. The full study input for Table 1 is identified
in Data availability.
