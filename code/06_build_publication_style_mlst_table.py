#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


GENE_TO_COLUMN = {
    "arcC": "arc",
    "aroE": "aro",
    "glpF": "glp",
    "gmk": "gmk",
    "pta": "pta",
    "tpi": "tpi",
    "yqiL": "yqi",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a publication-style MLST table from assembly, mlst, and ST guess outputs."
    )
    parser.add_argument(
        "--assembly-summary",
        required=True,
        help="Path to assembly_polish_summary.tsv",
    )
    parser.add_argument(
        "--mlst-summary",
        required=True,
        help="Path to mlst_sample_summary.tsv",
    )
    parser.add_argument(
        "--st-guess-summary",
        required=True,
        help="Path to st_guess_summary.tsv",
    )
    parser.add_argument(
        "--output-path",
        required=True,
        help="Path to output TSV table.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"Missing input file: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def sample_number(sample_name: str) -> int:
    match = re.search(r"(\d+)$", sample_name)
    if not match:
        raise SystemExit(f"Could not parse sample number from: {sample_name}")
    return int(match.group(1))


def format_strain_type(mlst_row: Dict[str, str], guess_row: Dict[str, str]) -> str:
    standard_st = mlst_row.get("ST", "").strip()
    if standard_st and standard_st != "-":
        return f"ST{standard_st}"
    best_guess = guess_row.get("best_guess_ST", "").strip()
    if best_guess:
        return f"ST{best_guess}"
    return "ST?"


def format_genotype_result(mlst_row: Dict[str, str], guess_row: Dict[str, str]) -> str:
    standard_st = mlst_row.get("ST", "").strip()
    if standard_st and standard_st != "-":
        return "typed_ST"

    best_guess = guess_row.get("best_guess_ST", "").strip()
    if best_guess:
        confidence = guess_row.get("guess_confidence", "").strip()
        return f"best_guess_unique:{confidence}" if confidence else "best_guess_unique"

    candidates = guess_row.get("best_guess_candidates", "").strip()
    if candidates:
        return f"ambiguous:{candidates}"

    typing_status = mlst_row.get("typing_status", "").strip()
    return typing_status or "unresolved"


def build_rows(
    assembly_rows: List[Dict[str, str]],
    mlst_rows: List[Dict[str, str]],
    guess_rows: List[Dict[str, str]],
) -> List[Dict[str, object]]:
    assembly_by_sample = {row["sample"]: row for row in assembly_rows}
    guess_by_sample = {row["sample"]: row for row in guess_rows}

    rows: List[Dict[str, object]] = []
    sorted_rows = sorted(mlst_rows, key=lambda row: sample_number(row["sample"]))
    total_rows = len(sorted_rows)

    print(f"[progress] step4 0/{total_rows} start", flush=True)

    for index, mlst_row in enumerate(sorted_rows, start=1):
        sample = mlst_row["sample"]
        print(f"[progress] step4 {index}/{total_rows} {sample}", flush=True)
        assembly_row = assembly_by_sample.get(sample, {})
        guess_row = guess_by_sample.get(sample, {})

        publication_row: Dict[str, object] = {
            "SampleNo.": sample_number(sample),
            "MLST strain type": format_strain_type(mlst_row, guess_row),
            "Reads number": assembly_row.get("read_count", ""),
            "Genotype Result": format_genotype_result(mlst_row, guess_row),
        }

        for gene, column_name in GENE_TO_COLUMN.items():
            publication_row[column_name] = mlst_row.get(f"{gene}_call", "")

        rows.append(publication_row)

    print(f"[progress] step4 {total_rows}/{total_rows} done", flush=True)
    return rows


def main() -> None:
    args = parse_args()

    assembly_rows = read_tsv(Path(args.assembly_summary))
    mlst_rows = read_tsv(Path(args.mlst_summary))
    guess_rows = read_tsv(Path(args.st_guess_summary))

    rows = build_rows(assembly_rows, mlst_rows, guess_rows)
    write_tsv(
        Path(args.output_path),
        rows,
        [
            "SampleNo.",
            "MLST strain type",
            "Reads number",
            "Genotype Result",
            "arc",
            "aro",
            "glp",
            "gmk",
            "pta",
            "tpi",
            "yqi",
        ],
    )

    print(f"Wrote {args.output_path}")


if __name__ == "__main__":
    main()
