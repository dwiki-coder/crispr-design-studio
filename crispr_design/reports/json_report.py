"""JSON report generation."""

from __future__ import annotations

import json
from datetime import datetime
from typing import Any, Optional

from ..core.design import DesignResult
from ..core.off_target import OffTargetHit
from ..core.scoring import ScoreResult


def generate_json_report(
    results: list[DesignResult],
    output_path: str,
    metadata: Optional[dict[str, Any]] = None,
) -> str:
    """Generate a JSON report from design results.

    Args:
        results: List of DesignResult objects.
        output_path: Path to write JSON file.
        metadata: Optional metadata to include.

    Returns:
        Path to the generated report.
    """
    report = {
        "report_type": "crispr_design",
        "generated": datetime.now().isoformat(),
        "tool_version": "0.1.0",
        "total_guides": len(results),
        "summary": {
            "total_guides": len(results),
            "top_guide": results[0].guide.to_dict() if results else None,
        },
        "metadata": metadata or {},
        "results": [r.to_dict() for r in results],
    }

    with open(output_path, "w") as f:
        json.dump(report, f, indent=2, default=str)

    return output_path


def format_hit_json(hit: OffTargetHit) -> dict:
    """Format a single off-target hit as JSON-serializable dict."""
    return {
        "chromosome": hit.chromosome,
        "position": hit.position,
        "strand": hit.strand,
        "target_sequence": hit.target_sequence,
        "mismatches": hit.mismatches,
        "mismatch_positions": hit.mismatch_positions,
        "mismatch_details": hit.mismatch_details,
        "gene_annotation": hit.gene_annotation,
    }


def format_scores_json(scores: ScoreResult) -> dict:
    """Format scores as JSON-serializable dict."""
    return scores.to_dict()
