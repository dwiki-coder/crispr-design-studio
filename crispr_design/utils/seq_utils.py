"""Sequence utility functions."""

from __future__ import annotations

import re
from typing import Optional


def reverse_complement(sequence: str) -> str:
    """Return reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        Reverse complement sequence.
    """
    comp = str.maketrans("ACGT", "TGCA")
    return sequence.upper()[::-1].translate(comp)


def validate_dna_sequence(sequence: str, allow_n: bool = False) -> bool:
    """Validate that a string contains only valid DNA bases.

    Args:
        sequence: String to validate.
        allow_n: Whether to allow 'N' (unknown base).

    Returns:
        True if valid, False otherwise.
    """
    valid = set("ACGT")
    if allow_n:
        valid = valid | set("N")
    return all(b.upper() in valid for b in sequence.strip())


def count_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence.

    Args:
        sequence: DNA sequence.

    Returns:
        GC content as fraction (0-1).
    """
    if not sequence:
        return 0.0
    gc = sum(1 for b in sequence.upper() if b in "GC")
    return gc / len(sequence)


def hamming_distance(seq1: str, seq2: str) -> int:
    """Compute Hamming distance between two equal-length sequences.

    Args:
        seq1: First sequence.
        seq2: Second sequence.

    Returns:
        Number of positions where bases differ.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be equal length")
    return sum(a != b for a, b in zip(seq1, seq2))


def parse_fasta(filepath: str) -> dict[str, str]:
    """Parse a FASTA file into a dictionary.

    Args:
        filepath: Path to FASTA file.

    Returns:
        Dict of {header: sequence}.
    """
    from Bio import SeqIO

    result = {}
    for record in SeqIO.parse(filepath, "fasta"):
        result[str(record.id)] = str(record.seq).upper()
    return result


def parse_vcf(filepath: str) -> list[dict]:
    """Parse a VCF file into a list of variant dictionaries.

    Args:
        filepath: Path to VCF file.

    Returns:
        List of dicts with variant information.
    """
    variants = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            variant = {
                "chrom": parts[0],
                "pos": int(parts[1]),
                "id": parts[2],
                "ref": parts[3],
                "alt": parts[4],
                "qual": parts[5],
                "filter": parts[6],
                "info": parts[7],
            }
            variants.append(variant)
    return variants


def extract_grna_from_fasta(
    filepath: str,
    guide_length: int = 20,
    pam_pattern: str = "NGG",
) -> list[dict]:
    """Extract potential gRNA sequences from a FASTA file.

    Args:
        filepath: Path to FASTA file.
        guide_length: Length of guide sequence.
        pam_pattern: PAM pattern to search for.

    Returns:
        List of dicts with gRNA info.
    """
    sequences = parse_fasta(filepath)
    guides = []

    for name, seq in sequences.items():
        seq_upper = seq.upper()
        pam_regex = pam_pattern.replace("N", ".")
        pattern = f".{{{guide_length}}}({pam_regex})"

        for match in re.finditer(pattern, seq_upper):
            guide_seq = seq_upper[match.start():match.start() + guide_length]
            pam_seq = match.group(1)

            if any(b not in "ACGT" for b in guide_seq):
                continue

            guides.append({
                "source": name,
                "sequence": guide_seq,
                "pam": pam_seq,
                "position": match.start() + 1,
            })

    return guides


def format_sequence(
    sequence: str,
    width: int = 60,
    show_positions: bool = True,
) -> str:
    """Format a DNA sequence for display.

    Args:
        sequence: DNA sequence.
        width: Characters per line.
        show_positions: Whether to show position numbers.

    Returns:
        Formatted string.
    """
    lines = []
    for i in range(0, len(sequence), width):
        chunk = sequence[i : i + width]
        if show_positions:
            lines.append(f"{i + 1:>6}  {chunk}")
        else:
            lines.append(chunk)
    return "\n".join(lines)
