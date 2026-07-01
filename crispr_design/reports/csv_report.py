"""CSV report generation."""

from __future__ import annotations

import csv
from typing import Optional

from ..core.design import DesignResult
from ..core.off_target import OffTargetHit


def generate_csv_report(
    results: list[DesignResult],
    output_path: str,
    include_off_targets: bool = True,
) -> str:
    """Generate a CSV report from design results.

    Args:
        results: List of DesignResult objects.
        output_path: Path to write CSV file.
        include_off_targets: Whether to include off-target details.

    Returns:
        Path to the generated report.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)

        # Header
        header = [
            "rank", "name", "sequence", "pam", "type",
            "gc_content", "mit_cf", "cfd", "doench_2014", "doench_2016",
            "specificity_score", "off_target_count",
        ]
        if include_off_targets:
            header.extend(["off_target_sites", "min_mismatches", "chromosome", "position"])

        writer.writerow(header)

        for r in results:
            row = [
                r.rank,
                r.guide.name,
                r.guide.sequence,
                r.guide.pam_sequence,
                r.guide.gRNA_type.value,
                round(r.guide.gc_content, 4),
                round(r.scores.mit_cf, 4),
                round(r.scores.cfd, 4),
                round(r.scores.doench_2014, 4),
                round(r.scores.doench_2016, 4),
                round(r.specificity_score, 4),
                len(r.off_targets),
            ]

            if include_off_targets:
                ot_sites = "; ".join(
                    f"{h.chromosome}:{h.position}" for h in r.off_targets[:10]
                )
                min_mm = min((h.mismatches for h in r.off_targets), default=0)
                row.extend([
                    ot_sites,
                    min_mm,
                    r.off_targets[0].chromosome if r.off_targets else "",
                    r.off_targets[0].position if r.off_targets else 0,
                ])

            writer.writerow(row)

    return output_path


def generate_off_target_csv(
    hits: list[OffTargetHit],
    output_path: str,
    guide_name: str = "",
) -> str:
    """Generate a CSV of off-target hits.

    Args:
        hits: List of OffTargetHit objects.
        output_path: Path to write CSV file.
        guide_name: Name of the guide for reference.

    Returns:
        Path to the generated CSV.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "guide", "chromosome", "position", "strand",
            "target_sequence", "mismatches", "mismatch_positions",
            "mismatch_details", "gene_annotation",
        ])

        for hit in hits:
            writer.writerow([
                guide_name,
                hit.chromosome,
                hit.position,
                hit.strand,
                hit.target_sequence,
                hit.mismatches,
                ";".join(str(p) for p in hit.mismatch_positions),
                json_dumps(hit.mismatch_details),
                hit.gene_annotation or "",
            ])

    return output_path


def json_dumps(obj) -> str:
    """Simple JSON-like serialization for CSV."""
    import json
    return json.dumps(obj)
