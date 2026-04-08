#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path
from typing import Iterable, TextIO


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def is_missing_like(value: str) -> bool:
    return value.strip().lower() in {"", "na", "n/a", "nan", "null", "none"}


def load_keep_ids(read_calls_path: Path, sample_column: str) -> set[str]:
    keep_ids: set[str] = set()
    with open_text(read_calls_path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in read-calls table: {read_calls_path}")
        missing = {"read_id", sample_column} - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Missing required column(s) in {read_calls_path}: {', '.join(sorted(missing))}"
            )
        for row in reader:
            sample_value = (row.get(sample_column) or "").strip()
            if not is_missing_like(sample_value):
                keep_ids.add((row["read_id"] or "").strip())
    if not keep_ids:
        raise ValueError(f"No sample-assigned read IDs found in {read_calls_path}")
    return keep_ids


def iter_fastq_records(handle: Iterable[str]):
    record: list[str] = []
    for line in handle:
        record.append(line)
        if len(record) == 4:
            yield record
            record = []
    if record:
        raise ValueError("Input FASTQ ended mid-record.")


def read_id_from_header(header_line: str) -> str:
    if not header_line.startswith("@"):
        raise ValueError(f"Invalid FASTQ header: {header_line.rstrip()}")
    return header_line[1:].strip().split()[0]


def build_demo_fastq(input_fastq: Path, output_fastq: Path, keep_ids: set[str]) -> tuple[int, int]:
    output_fastq.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    total = 0

    with open_text(input_fastq, "r") as in_handle, open_text(output_fastq, "w") as out_handle:
        for record in iter_fastq_records(in_handle):
            total += 1
            record_id = read_id_from_header(record[0])
            if record_id in keep_ids:
                out_handle.writelines(record)
                kept += 1

    return kept, total


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create a smaller demo FASTQ by keeping only reads that the existing "
            "read_calls table already assigned to a sample."
        )
    )
    parser.add_argument("--input-fastq", required=True, help="Original FASTQ(.gz) path.")
    parser.add_argument(
        "--read-calls",
        required=True,
        help="Path to Step 1 read_calls TSV for the same input FASTQ.",
    )
    parser.add_argument("--output-fastq", required=True, help="Output demo FASTQ(.gz) path.")
    parser.add_argument(
        "--sample-column",
        default="sample",
        help="Column in read_calls used to decide whether a read is sample-assigned.",
    )
    args = parser.parse_args()

    input_fastq = Path(args.input_fastq)
    read_calls = Path(args.read_calls)
    output_fastq = Path(args.output_fastq)

    keep_ids = load_keep_ids(read_calls, args.sample_column)
    kept, total = build_demo_fastq(input_fastq, output_fastq, keep_ids)

    print(f"Loaded {len(keep_ids)} sample-assigned read IDs from {read_calls}")
    print(f"Wrote {kept} of {total} FASTQ records to {output_fastq}")
    if kept != len(keep_ids):
        missing = len(keep_ids) - kept
        raise SystemExit(
            f"Error: {missing} read IDs from {read_calls} were not found in {input_fastq}"
        )


if __name__ == "__main__":
    main()
