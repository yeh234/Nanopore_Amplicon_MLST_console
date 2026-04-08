#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Torsten Seemann's mlst as the standard Step 6 caller and normalize outputs."
    )
    parser.add_argument(
        "--assembly-dir",
        default="output/minimal_assembly_polish_v1",
        help="Directory containing per-sample assembly outputs.",
    )
    parser.add_argument(
        "--output-dir",
        default="output/minimal_mlst_standard_v1",
        help="Directory for raw mlst outputs and normalized summaries.",
    )
    parser.add_argument(
        "--scheme",
        default="saureus",
        help="MLST scheme passed to mlst --scheme.",
    )
    parser.add_argument(
        "--mlst-bin",
        default="mlst",
        help="Path to mlst executable. Defaults to mlst in PATH.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads passed to mlst.",
    )
    parser.add_argument(
        "--samples",
        default="",
        help="Optional comma-separated sample filter, for example sample1,sample2",
    )
    parser.add_argument(
        "--prefer-polished",
        action="store_true",
        default=True,
        help="Prefer *_assembly_polish_consensus.fasta over draft assemblies when both exist.",
    )
    parser.add_argument(
        "--skipcheck",
        action="store_true",
        help="Pass --skipcheck to mlst.",
    )
    return parser.parse_args()


def discover_assemblies(assembly_dir: Path, sample_filter: Sequence[str]) -> Dict[str, Path]:
    preferred: Dict[str, Path] = {}

    for sample_dir in sorted(path for path in assembly_dir.glob("sample*") if path.is_dir()):
        sample_name = sample_dir.name
        if sample_filter and sample_name not in sample_filter:
            continue

        polished = sorted(sample_dir.glob(f"{sample_name}_assembly_polish_consensus.fasta"))
        draft = sorted(sample_dir.glob(f"{sample_name}_assembly.fasta"))
        if polished:
            preferred[sample_name] = polished[0]
        elif draft:
            preferred[sample_name] = draft[0]

    return preferred


def run_command(command: Sequence[str]) -> None:
    subprocess.run(command, check=True)


def parse_mlst_call(value: str) -> Tuple[str, str]:
    raw = (value or "").strip()
    if raw in {"", "-"}:
        return "", "missing"
    if raw.startswith("~") and raw.endswith("?"):
        return raw[1:-1], "novel_partial_match"
    if raw.endswith("?"):
        return raw[:-1], "partial_match"
    if raw.startswith("~"):
        return raw[1:], "novel_full_length"
    return raw, "exact_match"


def normalize_sample_name(filename: str) -> str:
    for suffix in ("_assembly_polish_consensus.fasta", "_assembly.fasta"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return Path(filename).stem


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_normalized_outputs(
    assembly_paths: Dict[str, Path],
    full_rows: List[Dict[str, str]],
    legacy_rows: List[Dict[str, str]],
    output_dir: Path,
) -> None:
    full_by_file = {row["FILE"]: row for row in full_rows}
    legacy_by_file = {row["FILE"]: row for row in legacy_rows}
    legacy_fields = list(legacy_rows[0].keys()) if legacy_rows else []
    gene_names = [field for field in legacy_fields if field not in {"FILE", "SCHEME", "ST"}]

    sample_rows: List[Dict[str, object]] = []
    gene_rows: List[Dict[str, object]] = []

    for sample_name, assembly_path in sorted(assembly_paths.items()):
        file_label = assembly_path.name
        full = full_by_file.get(file_label, {})
        legacy = legacy_by_file.get(file_label, {})

        sample_row: Dict[str, object] = {
            "sample": sample_name,
            "assembly_fasta": str(assembly_path),
            "scheme": full.get("SCHEME", legacy.get("SCHEME", "")),
            "ST": full.get("ST", legacy.get("ST", "")),
            "STATUS": full.get("STATUS", ""),
            "SCORE": full.get("SCORE", ""),
        }
        exact_loci = 0
        nonmissing_loci = 0

        for gene in gene_names:
            raw_call = legacy.get(gene, "")
            allele_id, status = parse_mlst_call(raw_call)
            if status == "exact_match":
                exact_loci += 1
            if status != "missing":
                nonmissing_loci += 1

            sample_row[f"{gene}_call"] = raw_call
            sample_row[f"{gene}_allele"] = allele_id
            sample_row[f"{gene}_status"] = status
            gene_rows.append(
                {
                    "sample": sample_name,
                    "assembly_fasta": str(assembly_path),
                    "scheme": sample_row["scheme"],
                    "gene": gene,
                    "raw_call": raw_call,
                    "allele_id": allele_id,
                    "status": status,
                }
            )

        sample_row["exact_loci"] = exact_loci
        sample_row["nonmissing_loci"] = nonmissing_loci
        st_value = str(sample_row["ST"]).strip()
        if st_value and st_value != "-":
            sample_row["typing_status"] = "typed_st"
        elif exact_loci == len(gene_names):
            sample_row["typing_status"] = "all_loci_exact_but_st_unresolved"
        elif nonmissing_loci == 0:
            sample_row["typing_status"] = "no_loci_detected"
        else:
            sample_row["typing_status"] = "partial_or_novel_loci"

        sample_rows.append(sample_row)

    sample_fields = [
        "sample",
        "assembly_fasta",
        "scheme",
        "ST",
        "STATUS",
        "SCORE",
    ]
    for gene in gene_names:
        sample_fields.extend([f"{gene}_call", f"{gene}_allele", f"{gene}_status"])
    sample_fields.extend(["exact_loci", "nonmissing_loci", "typing_status"])

    gene_fields = [
        "sample",
        "assembly_fasta",
        "scheme",
        "gene",
        "raw_call",
        "allele_id",
        "status",
    ]

    write_tsv(output_dir / "mlst_sample_summary.tsv", sample_rows, sample_fields)
    write_tsv(output_dir / "mlst_gene_calls.tsv", gene_rows, gene_fields)


def main() -> None:
    args = parse_args()
    assembly_dir = Path(args.assembly_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not assembly_dir.exists():
        raise SystemExit(f"Assembly directory not found: {assembly_dir}")

    sample_filter = [item for item in args.samples.split(",") if item]
    assembly_paths = discover_assemblies(assembly_dir, sample_filter)
    if not assembly_paths:
        raise SystemExit(f"No sample assemblies found under {assembly_dir}")

    fofn_path = output_dir / "assemblies_fofn.txt"
    fofn_path.write_text("".join(f"{path}\n" for path in assembly_paths.values()))

    full_path = output_dir / "mlst_full.tsv"
    legacy_path = output_dir / "mlst_legacy.tsv"

    base_command = [args.mlst_bin]
    if args.skipcheck:
        base_command.append("--skipcheck")

    full_command = base_command + [
        "--full",
        "--nopath",
        "--scheme",
        args.scheme,
        "--threads",
        str(args.threads),
        "--outfile",
        str(full_path),
    ] + [str(path) for path in assembly_paths.values()]

    legacy_command = base_command + [
        "--legacy",
        "--nopath",
        "--scheme",
        args.scheme,
        "--threads",
        str(args.threads),
        "--outfile",
        str(legacy_path),
    ] + [str(path) for path in assembly_paths.values()]

    run_command(full_command)
    run_command(legacy_command)

    full_rows = read_tsv(full_path)
    legacy_rows = read_tsv(legacy_path)
    build_normalized_outputs(assembly_paths, full_rows, legacy_rows, output_dir)

    print(f"Wrote {fofn_path}")
    print(f"Wrote {full_path}")
    print(f"Wrote {legacy_path}")
    print(f"Wrote {output_dir / 'mlst_sample_summary.tsv'}")
    print(f"Wrote {output_dir / 'mlst_gene_calls.tsv'}")


if __name__ == "__main__":
    main()
