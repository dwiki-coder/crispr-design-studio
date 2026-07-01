"""Reference genome handling with FASTA indexing and sequence retrieval."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import Optional

from ..core.grna import Grna, GrnaType

logger = logging.getLogger(__name__)


class Genome:
    """Reference genome loaded from FASTA file.

    Uses pyfaidx for efficient random access to genome sequences.
    """

    def __init__(self, fasta_path: str | Path, index: bool = True):
        """Load a genome from FASTA file.

        Args:
            fasta_path: Path to FASTA file.
            index: Whether to create an X-index for fast access.
        """
        self.fasta_path = Path(fasta_path)
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        self._fa = None
        self._sequences: dict[str, str] = {}
        self._loaded: bool = False

        # Try pyfaidx first, fall back to manual loading
        pyfaidx_available = False
        try:
            from pyfaidx import Fasta
            self._fa = Fasta(str(self.fasta_path), as_text=True)
            self._loaded = True
            pyfaidx_available = True
        except ImportError:
            logger.info("pyfaidx not installed, using manual FASTA loading")
        except Exception as e:
            logger.warning(f"pyfaidx failed ({e}), loading manually")

        if not pyfaidx_available:
            self._load_fasta_manual()

        self._chromosomes: list[str] = list(self._sequences.keys())

    def _load_fasta_manual(self):
        """Load FASTA without pyfaidx (full load into memory)."""
        from Bio import SeqIO

        self._sequences.clear()
        for record in SeqIO.parse(str(self.fasta_path), "fasta"):
            self._sequences[str(record.id)] = str(record.seq).upper()

    @property
    def chromosomes(self) -> list[str]:
        """List of chromosome/contig names."""
        if self._fa is not None:
            return [str(c) for c in self._fa.keys()]
        return list(self._sequences.keys())

    @property
    def chromosome_sizes(self) -> dict[str, int]:
        """Get sizes of all chromosomes."""
        sizes = {}
        for chrom in self.chromosomes:
            sizes[chrom] = len(self.get_sequence(chrom))
        return sizes

    def get_sequence(
        self,
        chrom: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Get sequence for a chromosome or region.

        Args:
            chrom: Chromosome name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).
            strand: '+' or '-' for strand-specific retrieval.

        Returns:
            DNA sequence string.
        """
        if self._fa is not None:
            if start is not None and end is not None:
                seq = str(self._fa[chrom][start:end].seq)
            else:
                seq = str(self._fa[chrom].seq)
        else:
            if chrom not in self._sequences:
                raise KeyError(f"Chromosome {chrom} not found in genome")
            full_seq = self._sequences[chrom]
            if start is not None and end is not None:
                seq = full_seq[start:end]
            else:
                seq = full_seq

        if strand == "-":
            seq = self._reverse_complement(seq)

        return seq.upper()

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Return reverse complement of a DNA sequence."""
        comp = str.maketrans("ACGT", "TGCA")
        return seq[::-1].translate(comp)

    def to_dict(self) -> dict[str, str]:
        """Convert entire genome to dictionary (memory intensive!).

        Returns:
            Dict of {chromosome: sequence}.
        """
        if self._fa is not None:
            return {str(c): str(self._fa[c].seq).upper() for c in self._fa.keys()}
        return dict(self._sequences)

    def get_region_context(
        self,
        chrom: str,
        position: int,
        flank_size: int = 50,
    ) -> dict:
        """Get genomic context around a position.

        Args:
            chrom: Chromosome name.
            position: 1-based genomic position.
            flank_size: Number of bases on each side.

        Returns:
            Dictionary with sequence context.
        """
        start = max(0, position - flank_size - 1)
        end = position + flank_size

        seq = self.get_sequence(chrom, start, end)
        return {
            "chromosome": chrom,
            "position": position,
            "flank_size": flank_size,
            "sequence": seq,
            "start": start,
            "end": end,
        }

    def close(self):
        """Close genome file handles."""
        if self._fa is not None:
            # pyfaidx Fasta doesn't have explicit close, but we can clear ref
            self._fa = None


def load_reference_genome(
    fasta_path: str | Path,
    chromosomes: Optional[list[str]] = None,
) -> dict[str, str]:
    """Load reference genome as dictionary.

    Args:
        fasta_path: Path to FASTA file.
        chromosomes: Optional list of chromosomes to load (None = all).

    Returns:
        Dict of {chromosome: sequence}.
    """
    genome = Genome(fasta_path)
    result = {}
    chroms = chromosomes or genome.chromosomes
    for chrom in chroms:
        result[chrom] = genome.get_sequence(chrom)
    genome.close()
    return result


def create_test_genome(sequence_length: int = 10000, seed: int = 42) -> dict[str, str]:
    """Create a synthetic genome for testing.

    Args:
        sequence_length: Length of each chromosome.
        seed: Random seed for reproducibility.

    Returns:
        Dict of {chromosome: sequence}.
    """
    import random
    rng = random.Random(seed)

    bases = "ACGT"
    genome = {}
    for i in range(1, 6):
        chrom = f"chr{i}"
        seq = "".join(rng.choices(bases, k=sequence_length))
        # Insert some known PAM sites (replace existing bases, don't change length)
        for j in range(100, len(seq) - 30, 500):
            pam = rng.choice(["GG", "AG"])
            seq = seq[:j] + rng.choice(bases) + rng.choice(bases) + pam + seq[j + 4:]
        genome[chrom] = seq
    return genome
