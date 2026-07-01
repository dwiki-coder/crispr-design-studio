"""Off-target search engine using seed-and-extend algorithm."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

from ..core.grna import Grna, GrnaType

logger = logging.getLogger(__name__)


@dataclass
class OffTargetHit:
    """Represents a single off-target match."""

    chromosome: str
    position: int
    strand: str
    target_sequence: str
    mismatches: int
    mismatch_positions: list[int] = field(default_factory=list)
    mismatch_details: list[dict] = field(default_factory=list)
    genomic_context: str = ""
    gene_annotation: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "chromosome": self.chromosome,
            "position": self.position,
            "strand": self.strand,
            "target_sequence": self.target_sequence,
            "mismatches": self.mismatches,
            "mismatch_positions": self.mismatch_positions,
            "mismatch_details": self.mismatch_details,
            "genomic_context": self.genomic_context,
            "gene_annotation": self.gene_annotation,
        }


def _count_mismatches(seq1: str, seq2: str) -> tuple[int, list[int], list[dict]]:
    """Compare two sequences and return mismatch count, positions, and details."""
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must be equal length: {len(seq1)} vs {len(seq2)}")

    count = 0
    positions = []
    details = []
    for i, (c1, c2) in enumerate(zip(seq1, seq2)):
        if c1 != c2:
            count += 1
            positions.append(i)
            details.append({"position": i, "query": c1, "target": c2})
    return count, positions, details


def _seed_match(query: str, sequence: str, seed_length: int, max_mismatches: int) -> bool:
    """Check if the seed region (3' end) matches within mismatch tolerance."""
    seed = query[-seed_length:]
    target_seed = sequence[-seed_length:]
    mismatches, _, _ = _count_mismatches(seed, target_seed)
    return mismatches <= max(0, max_mismatches - 1)


def find_off_targets(
    guide: Grna,
    genome: dict[str, str],
    max_mismatches: int = 3,
    seed_length: int = 10,
    max_results: int = 1000,
    include_on_target: bool = True,
) -> list[OffTargetHit]:
    """Find off-target sites using seed-and-extend algorithm.

    The algorithm:
    1. Index genome by seed k-mers for fast lookup
    2. For each seed match, extend to full guide and count mismatches
    3. Filter by total mismatch tolerance

    Args:
        guide: The gRNA to search for.
        genome: Dictionary of {chromosome: sequence}.
        max_mismatches: Maximum mismatches allowed.
        seed_length: Length of seed region (3' end) for fast filtering.
        max_results: Maximum number of results to return.
        include_on_target: Whether to include the on-target site.

    Returns:
        List of OffTargetHit sorted by mismatch count.
    """
    guide_seq = guide.sequence.upper()
    rev_guide_seq = guide.sequence.upper()
    pam_seq = guide.pam_sequence.upper()
    search_seq = guide_seq + pam_seq
    search_len = len(search_seq)
    seed_k = min(seed_length, len(search_seq) - len(pam_seq))

    # Build seed index: map seed k-mers -> list of (chrom, pos)
    seed_index: dict[str, list[tuple[str, int]]] = {}
    for chrom, seq in genome.items():
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - seed_k + 1):
            kmer = seq_upper[i : i + seed_k]
            if "N" not in kmer:
                if kmer not in seed_index:
                    seed_index[kmer] = []
                seed_index[kmer].append((chrom, i))

    # Seed from the 3' end of the guide (most critical region)
    guide_seed = search_seq[-(seed_k + len(pam_seq)):]
    if len(guide_seed) < seed_k:
        guide_seed = search_seq
    seed_kmer = guide_seed[:seed_k]

    results: list[OffTargetHit] = []
    on_target_found = False

    # Check forward strand
    if seed_kmer in seed_index:
        for chrom, pos in seed_index[seed_kmer]:
            start = max(0, pos - len(pam_seq))
            end = min(len(genome[chrom]), start + search_len)
            if end - start < search_len:
                continue
            target = genome[chrom][start:end]

            if not _seed_match(search_seq, target, seed_length, max_mismatches):
                continue

            mismatches, positions, details = _count_mismatches(search_seq, target)
            if mismatches <= max_mismatches:
                if not include_on_target and mismatches == 0:
                    on_target_found = True
                    continue
                if mismatches == 0:
                    on_target_found = True

                hit = OffTargetHit(
                    chromosome=chrom,
                    position=start + 1,  # 1-based
                    strand="+",
                    target_sequence=target,
                    mismatches=mismatches,
                    mismatch_positions=positions,
                    mismatch_details=details,
                )
                results.append(hit)

                if len(results) >= max_results:
                    return results

    # Check reverse complement
    rc_seq = guide.reverse_complement() + guide.pam_sequence
    rc_seed = rc_seq[-(seed_k + len(pam_seq)):]
    if len(rc_seed) < seed_k:
        rc_seed = rc_seq
    rc_seed_kmer = rc_seed[:seed_k]

    if rc_seed_kmer in seed_index:
        for chrom, pos in seed_index[rc_seed_kmer]:
            start = max(0, pos - len(pam_seq))
            end = min(len(genome[chrom]), start + search_len)
            if end - start < search_len:
                continue
            target = genome[chrom][start:end]

            if not _seed_match(rc_seq, target, seed_length, max_mismatches):
                continue

            mismatches, positions, details = _count_mismatches(rc_seq, target)
            if mismatches <= max_mismatches:
                if not include_on_target and mismatches == 0:
                    continue
                hit = OffTargetHit(
                    chromosome=chrom,
                    position=start + 1,
                    strand="-",
                    target_sequence=target,
                    mismatches=mismatches,
                    mismatch_positions=positions,
                    mismatch_details=details,
                )
                results.append(hit)

                if len(results) >= max_results:
                    return results

    # Sort by mismatch count, then by chromosome and position
    results.sort(key=lambda h: (h.mismatches, h.chromosome, h.position))
    return results


def find_off_targets_brute_force(
    guide: Grna,
    genome: dict[str, str],
    max_mismatches: int = 3,
    include_on_target: bool = True,
    max_results: int = 1000,
) -> list[OffTargetHit]:
    """Find off-target sites using brute-force sliding window (for small genomes).

    This is slower than seed-and-extend but guarantees finding all matches.

    Args:
        guide: The gRNA to search for.
        genome: Dictionary of {chromosome: sequence}.
        max_mismatches: Maximum mismatches allowed.
        include_on_target: Whether to include the on-target site.
        max_results: Maximum results to return.

    Returns:
        List of OffTargetHit sorted by mismatch count.
    """
    search_seq = guide.sequence.upper() + guide.pam_sequence.upper()
    search_len = len(search_seq)
    results: list[OffTargetHit] = []

    for chrom, seq in genome.items():
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - search_len + 1):
            target = seq_upper[i : i + search_len]
            if "N" in target:
                continue

            mismatches, positions, details = _count_mismatches(search_seq, target)
            if mismatches <= max_mismatches:
                if not include_on_target and mismatches == 0:
                    continue

                hit = OffTargetHit(
                    chromosome=chrom,
                    position=i + 1,
                    strand="+",
                    target_sequence=target,
                    mismatches=mismatches,
                    mismatch_positions=positions,
                    mismatch_details=details,
                )
                results.append(hit)

                if len(results) >= max_results:
                    results.sort(key=lambda h: (h.mismatches, h.chromosome, h.position))
                    return results

    results.sort(key=lambda h: (h.mismatches, h.chromosome, h.position))
    return results
