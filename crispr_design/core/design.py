"""Guide optimization and design utilities."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

from .grna import Grna, GrnaType, find_pams_in_sequence
from .off_target import OffTargetHit, find_off_targets
from .scoring import ScoreResult, compute_all_scores

logger = logging.getLogger(__name__)


@dataclass
class DesignResult:
    """Result of a gRNA design scan."""

    guide: Grna
    scores: ScoreResult
    off_targets: list[OffTargetHit] = field(default_factory=list)
    specificity_score: float = 0.0
    rank: int = 0

    def to_dict(self) -> dict:
        return {
            "rank": self.rank,
            "guide": self.guide.to_dict(),
            "scores": self.scores.to_dict(),
            "off_target_count": len(self.off_targets),
            "specificity_score": round(self.specificity_score, 4),
            "off_targets": [hit.to_dict() for hit in self.off_targets[:10]],
        }

    @property
    def off_target_count(self) -> int:
        """Number of off-target sites found."""
        return len(self.off_targets)


def design_guides(
    sequence: str,
    gRNA_type: GrnaType = GrnaType.SPCAS9,
    genome: Optional[dict[str, str]] = None,
    max_mismatches: int = 3,
    min_specificity: float = 0.0,
    max_results: int = 100,
    guide_name_prefix: str = "",
    chromosome: str = "",
    position: int = 0,
) -> list[DesignResult]:
    """Design and rank gRNAs from a DNA sequence.

    Scans the sequence for all valid PAM sites, extracts gRNAs,
    scores them, and optionally checks for off-target sites.

    Args:
        sequence: DNA sequence to scan.
        gRNA_type: Cas enzyme type.
        genome: Reference genome for off-target search.
        max_mismatches: Max mismatches for off-target search.
        min_specificity: Minimum specificity score threshold.
        max_results: Maximum number of results to return.
        guide_name_prefix: Prefix for guide names.
        chromosome: Chromosome name (for context).
        position: Genomic position (for context).

    Returns:
        List of DesignResult sorted by specificity.
    """
    guides = find_pams_in_sequence(sequence, gRNA_type)
    results: list[DesignResult] = []

    for idx, guide in enumerate(guides):
        if guide_name_prefix:
            guide.name = f"{guide_name_prefix}_{idx + 1}"
        else:
            guide.name = f"guide_{idx + 1}_{guide.sequence[:8]}"

        # Compute on-target scores
        scores = compute_all_scores(guide.sequence)

        # Find off-targets if genome provided
        off_targets: list[OffTargetHit] = []
        if genome:
            off_targets = find_off_targets(
                guide,
                genome,
                max_mismatches=max_mismatches,
                include_on_target=False,
                max_results=50,
            )

        # Compute specificity score
        specificity = _compute_specificity(scores, off_targets, max_mismatches)

        if specificity < min_specificity:
            continue

        result = DesignResult(
            guide=guide,
            scores=scores,
            off_targets=off_targets,
            specificity_score=specificity,
        )
        results.append(result)

    # Sort by specificity score (higher = better)
    results.sort(key=lambda r: r.specificity_score, reverse=True)

    # Assign ranks
    for i, r in enumerate(results):
        r.rank = i + 1

    return results[:max_results]


def _compute_specificity(scores: ScoreResult, off_targets: list[OffTargetHit], max_mismatches: int) -> float:
    """Compute a composite specificity score.

    Combines on-target activity prediction with off-target risk assessment.

    Args:
        scores: Scoring result.
        off_targets: List of off-target hits.
        max_mismatches: Maximum mismatches considered.

    Returns:
        Specificity score (0-100, higher = more specific).
    """
    # Penalize for off-target sites (weighted by proximity)
    ot_penalty = 0.0
    for hit in off_targets:
        weight = max(0, 1.0 - (hit.mismatches / (max_mismatches + 1)))
        ot_penalty += weight * 10.0

    # Normalize Doench score to 0-100 scale
    # Doench scores typically range from -20 to +40
    doench_norm = max(0, min(100, (scores.doench_2016 + 20) * 100 / 60))

    # Combine: start with activity score, subtract off-target penalty
    specificity = doench_norm - ot_penalty
    return max(0.0, min(100.0, specificity))


def filter_guides(
    results: list[DesignResult],
    min_gc: float = 0.3,
    max_gc: float = 0.7,
    min_doench: float = -10.0,
    max_off_targets: int = 10,
) -> list[DesignResult]:
    """Filter design results based on quality criteria.

    Args:
        results: Design results to filter.
        min_gc: Minimum GC content.
        max_gc: Maximum GC content.
        min_doench: Minimum Doench score.
        max_off_targets: Maximum number of off-targets allowed.

    Returns:
        Filtered list of design results.
    """
    filtered = []
    for r in results:
        gc = r.guide.gc_content
        if gc < min_gc or gc > max_gc:
            continue
        if r.scores.doench_2016 < min_doench:
            continue
        if len(r.off_targets) > max_off_targets:
            continue
        filtered.append(r)
    return filtered
