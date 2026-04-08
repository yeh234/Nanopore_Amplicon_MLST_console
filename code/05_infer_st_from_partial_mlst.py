#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


GENES = ["arcC", "aroE", "glpF", "gmk", "pta", "tpi", "yqiL"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Infer nearest known ST candidates from partial or novel MLST calls."
    )
    parser.add_argument(
        "--mlst-summary",
        required=True,
        help="Path to mlst_sample_summary.tsv produced by 05_minimal_mlst_standard_caller.py",
    )
    parser.add_argument(
        "--profile-table",
        default="db/saureus_Profile.txt",
        help="Path to pubMLST profile table.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for inferred ST outputs.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="How many candidate STs to show in summary candidate lists.",
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


def joined(values: Sequence[str]) -> str:
    return ";".join(values)


def summarize_candidates(candidates: Sequence[str], top_n: int) -> str:
    if len(candidates) <= top_n:
        return ";".join(candidates)
    return ";".join(candidates[:top_n]) + f";...(+{len(candidates) - top_n} more)"


def confidence_label(top_count: int, exact_supported: int, approx_supported: int) -> str:
    if top_count == 1:
        if exact_supported >= 6:
            return "high_unique"
        if exact_supported >= 5:
            return "medium_unique"
        if exact_supported >= 4 and approx_supported >= 1:
            return "medium_unique"
        return "low_unique"
    if exact_supported >= 6 and top_count <= 3:
        return "medium_ambiguous"
    if exact_supported >= 5 and top_count <= 5:
        return "medium_ambiguous"
    return "ambiguous"


def score_profile(sample_row: Dict[str, str], profile_row: Dict[str, str]) -> Dict[str, object]:
    exact_supported_genes: List[str] = []
    approx_supported_genes: List[str] = []
    exact_conflict_genes: List[str] = []
    approx_conflict_genes: List[str] = []
    missing_genes: List[str] = []

    for gene in GENES:
        allele_id = sample_row.get(f"{gene}_allele", "").strip()
        status = sample_row.get(f"{gene}_status", "").strip()
        if not allele_id:
            missing_genes.append(gene)
            continue

        if status == "exact_match":
            if profile_row[gene] == allele_id:
                exact_supported_genes.append(gene)
            else:
                exact_conflict_genes.append(gene)
        else:
            if profile_row[gene] == allele_id:
                approx_supported_genes.append(gene)
            else:
                approx_conflict_genes.append(gene)

    score_tuple = (
        len(exact_supported_genes),
        len(approx_supported_genes),
        -len(exact_conflict_genes),
        -len(approx_conflict_genes),
    )

    return {
        "candidate_ST": profile_row["ST"],
        "score_tuple": score_tuple,
        "exact_supported_genes": exact_supported_genes,
        "approx_supported_genes": approx_supported_genes,
        "exact_conflict_genes": exact_conflict_genes,
        "approx_conflict_genes": approx_conflict_genes,
        "missing_genes": missing_genes,
    }


def infer_candidates(
    mlst_rows: List[Dict[str, str]],
    profile_rows: List[Dict[str, str]],
    top_n: int,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    summary_rows: List[Dict[str, object]] = []
    candidate_rows: List[Dict[str, object]] = []
    total_rows = len(mlst_rows)

    print(f"[progress] step3 0/{total_rows} start", flush=True)

    for index, sample_row in enumerate(mlst_rows, start=1):
        sample = sample_row["sample"]
        print(f"[progress] step3 {index}/{total_rows} {sample}", flush=True)
        scored = [score_profile(sample_row, profile_row) for profile_row in profile_rows]
        scored.sort(key=lambda item: (item["score_tuple"], item["candidate_ST"]), reverse=True)

        top_score = scored[0]["score_tuple"]
        top_candidates = [item for item in scored if item["score_tuple"] == top_score]
        top_candidate_ids = [item["candidate_ST"] for item in top_candidates]
        lead = top_candidates[0]

        summary_rows.append(
            {
                "sample": sample,
                "source_ST": sample_row.get("ST", ""),
                "source_typing_status": sample_row.get("typing_status", ""),
                "exact_loci": sample_row.get("exact_loci", ""),
                "nonmissing_loci": sample_row.get("nonmissing_loci", ""),
                "best_guess_ST": top_candidate_ids[0] if len(top_candidate_ids) == 1 else "",
                "best_guess_candidates": summarize_candidates(top_candidate_ids, top_n),
                "top_rank_candidate_count": len(top_candidate_ids),
                "matched_exact_loci": len(lead["exact_supported_genes"]),
                "matched_approx_loci": len(lead["approx_supported_genes"]),
                "conflicting_exact_loci": len(lead["exact_conflict_genes"]),
                "conflicting_approx_loci": len(lead["approx_conflict_genes"]),
                "matched_exact_genes": joined(lead["exact_supported_genes"]),
                "matched_approx_genes": joined(lead["approx_supported_genes"]),
                "conflicting_exact_genes": joined(lead["exact_conflict_genes"]),
                "conflicting_approx_genes": joined(lead["approx_conflict_genes"]),
                "guess_confidence": confidence_label(
                    len(top_candidate_ids),
                    len(lead["exact_supported_genes"]),
                    len(lead["approx_supported_genes"]),
                ),
            }
        )

        for rank, item in enumerate(scored[:top_n], start=1):
            candidate_rows.append(
                {
                    "sample": sample,
                    "rank": rank,
                    "candidate_ST": item["candidate_ST"],
                    "is_top_rank": "yes" if item["score_tuple"] == top_score else "no",
                    "matched_exact_loci": len(item["exact_supported_genes"]),
                    "matched_approx_loci": len(item["approx_supported_genes"]),
                    "conflicting_exact_loci": len(item["exact_conflict_genes"]),
                    "conflicting_approx_loci": len(item["approx_conflict_genes"]),
                    "matched_exact_genes": joined(item["exact_supported_genes"]),
                    "matched_approx_genes": joined(item["approx_supported_genes"]),
                    "conflicting_exact_genes": joined(item["exact_conflict_genes"]),
                    "conflicting_approx_genes": joined(item["approx_conflict_genes"]),
                }
            )

    print(f"[progress] step3 {total_rows}/{total_rows} done", flush=True)
    return summary_rows, candidate_rows


def main() -> None:
    args = parse_args()
    mlst_rows = read_tsv(Path(args.mlst_summary))
    profile_rows = read_tsv(Path(args.profile_table))
    output_dir = Path(args.output_dir)

    summary_rows, candidate_rows = infer_candidates(mlst_rows, profile_rows, args.top_n)

    write_tsv(
        output_dir / "st_guess_summary.tsv",
        summary_rows,
        [
            "sample",
            "source_ST",
            "source_typing_status",
            "exact_loci",
            "nonmissing_loci",
            "best_guess_ST",
            "best_guess_candidates",
            "top_rank_candidate_count",
            "matched_exact_loci",
            "matched_approx_loci",
            "conflicting_exact_loci",
            "conflicting_approx_loci",
            "matched_exact_genes",
            "matched_approx_genes",
            "conflicting_exact_genes",
            "conflicting_approx_genes",
            "guess_confidence",
        ],
    )
    write_tsv(
        output_dir / "st_guess_top_candidates.tsv",
        candidate_rows,
        [
            "sample",
            "rank",
            "candidate_ST",
            "is_top_rank",
            "matched_exact_loci",
            "matched_approx_loci",
            "conflicting_exact_loci",
            "conflicting_approx_loci",
            "matched_exact_genes",
            "matched_approx_genes",
            "conflicting_exact_genes",
            "conflicting_approx_genes",
        ],
    )

    print(f"Wrote {output_dir / 'st_guess_summary.tsv'}")
    print(f"Wrote {output_dir / 'st_guess_top_candidates.tsv'}")


if __name__ == "__main__":
    main()
