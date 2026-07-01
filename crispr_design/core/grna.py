"""gRNA class with PAM detection and sequence utilities."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

logger = logging.getLogger(__name__)


class GrnaType(str, Enum):
    """Supported gRNA types based on Cas enzyme."""

    SPCAS9 = "spcas9"       # S. pyogenes Cas9
    SACAS9 = "sacas9"       # S. aureus Cas9
    SPCAS9N = "spcas9n"     # SpCas9 nickase
    HIFICAS9 = "hificas9"   # High-fidelity Cas9
    SPCAS9V = "spcas9v"     # SpCas9-HySt (hypersensitive)
    XMACAS9 = "xmacas9"     # xMACs Cas9


class PAMType(str, Enum):
    """PAM pattern type."""

    NGG = "NGG"             # SpCas9 canonical
    NNGRRT = "NNGRRT"       # SaCas9
    NNG = "NNG"             # SpCas9 relaxed
    NRCH = "NRCH"           # HiFi variants
    NRY = "NRY"             # SaCas9 variant


# PAM definitions per gRNA type
PAM_DEFINITIONS = {
    GrnaType.SPCAS9: {"pam": r"([ACGT])([AC][GT])([AC][GT])", "regex": r"NN(?=[ACGT]{2}$)", "pattern": "NGG", "length": 2, "name": "NGG"},
    GrnaType.SACAS9: {"pam": r"NNNRRT", "regex": r"NNN.{3}", "pattern": "NNGRRT", "length": 6, "name": "NNGRRT"},
    GrnaType.SPCAS9N: {"pam": r"NN(?=[ACGT]{2}$)", "regex": r"NNN", "pattern": "NGG", "length": 3, "name": "NGG/NAG"},
    GrnaType.HIFICAS9: {"pam": r"NN(?=[ACGT]{2}$)", "regex": r"NN(?=GG)", "pattern": "NGG", "length": 2, "name": "NGG"},
    GrnaType.SPCAS9V: {"pam": r"NNN", "regex": r"NNN", "pattern": "NRY", "length": 3, "name": "NRY"},
    GrnaType.XMACAS9: {"pam": r"NN(?=[ACGT]{2}$)", "regex": r"NN(?=[ACGT]{2}$)", "pattern": "NGG", "length": 2, "name": "NGG"},
}

# Default guide lengths per type
GUIDE_LENGTHS = {
    GrnaType.SPCAS9: 20,
    GrnaType.SACAS9: 20,
    GrnaType.SPCAS9N: 20,
    GrnaType.HIFICAS9: 20,
    GrnaType.SPCAS9V: 20,
    GrnaType.XMACAS9: 20,
}


@dataclass
class Grna:
    """CRISPR guide RNA with PAM detection and scoring.

    Attributes:
        sequence: The gRNA target sequence (20bp for SpCas9, 20bp for SaCas9).
        pam_sequence: The PAM sequence immediately downstream.
        gRNA_type: Type of Cas enzyme.
        name: Optional identifier for the guide.
        chromosome: Chromosome name (if from genomic coordinates).
        position: Genomic position (1-based, if known).
        strand: '+' or '-'.
    """

    sequence: str
    pam_sequence: str = ""
    gRNA_type: GrnaType = GrnaType.SPCAS9
    name: str = ""
    chromosome: str = ""
    position: int = 0
    strand: str = "+"

    _validated: bool = field(default=False, repr=False)

    def __post_init__(self):
        self.sequence = self.sequence.upper().strip()
        self.pam_sequence = self.pam_sequence.upper().strip()
        if not self._validated:
            self._validate()

    def _validate(self) -> None:
        """Validate the gRNA sequence and PAM."""
        valid_bases = set("ACGT")
        seq_clean = self.sequence.replace("T", "U")  # Allow U in RNA notation

        if not self.sequence:
            raise ValueError("gRNA sequence cannot be empty")

        guide_len = GUIDE_LENGTHS[self.gRNA_type]
        if len(self.sequence) != guide_len:
            raise ValueError(
                f"gRNA sequence length must be {guide_len} for {self.gRNA_type.value}, "
                f"got {len(self.sequence)}"
            )

        for base in self.sequence:
            if base not in valid_bases:
                raise ValueError(f"Invalid base '{base}' in gRNA sequence. Expected A, C, G, T.")

        if self.pam_sequence:
            for base in self.pam_sequence:
                if base not in valid_bases:
                    raise ValueError(f"Invalid base '{base}' in PAM sequence. Expected A, C, G, T.")

        self._validated = True

    @classmethod
    def from_target(
        cls,
        target_seq: str,
        gRNA_type: GrnaType = GrnaType.SPCAS9,
        name: str = "",
        chromosome: str = "",
        position: int = 0,
        strand: str = "+",
    ) -> "Grna":
        """Create a Grna from a full target sequence (guide + PAM).

        The target_seq should contain the guide sequence followed by the PAM.
        The PAM is automatically detected based on the gRNA_type.

        Args:
            target_seq: Full sequence including guide + PAM region.
            gRNA_type: Type of Cas enzyme.
            name: Optional name for the guide.
            chromosome: Chromosome name.
            position: Genomic position.
            strand: DNA strand.

        Returns:
            Grna instance with auto-detected PAM.

        Raises:
            ValueError: If no valid PAM is found.
        """
        target_seq = target_seq.upper().strip()
        guide_len = GUIDE_LENGTHS[gRNA_type]

        # Try to find PAM from the 3' end
        pam_info = PAM_DEFINITIONS[gRNA_type]
        pam_len = pam_info["length"]

        if len(target_seq) < guide_len + pam_len:
            raise ValueError(
                f"Target sequence too short. Need at least {guide_len + pam_len} bases "
                f"(guide={guide_len}, pam={pam_len}), got {len(target_seq)}."
            )

        # Extract potential guide and PAM
        guide_seq = target_seq[:guide_len]
        pam_seq = target_seq[guide_len:guide_len + pam_len]

        # Validate PAM
        if not cls._validate_pam(pam_seq, gRNA_type):
            # Try sliding window
            for start in range(len(target_seq) - pam_len - guide_len + 1):
                g = target_seq[start:start + guide_len]
                p = target_seq[start + guide_len:start + guide_len + pam_len]
                if cls._validate_pam(p, gRNA_type):
                    return cls(
                        sequence=g,
                        pam_sequence=p,
                        gRNA_type=gRNA_type,
                        name=name,
                        chromosome=chromosome,
                        position=position,
                        strand=strand,
                    )
            raise ValueError(
                f"No valid {pam_info['name']} PAM found in target sequence '{target_seq}'"
            )

        return cls(
            sequence=guide_seq,
            pam_sequence=pam_seq,
            gRNA_type=gRNA_type,
            name=name,
            chromosome=chromosome,
            position=position,
            strand=strand,
        )

    @classmethod
    def _validate_pam(cls, pam: str, gRNA_type: GrnaType) -> bool:
        """Validate PAM sequence for the given gRNA type."""
        if not pam:
            return False

        pam_upper = pam.upper()

        if gRNA_type in (GrnaType.SPCAS9, GrnaType.HIFICAS9, GrnaType.XMACAS9):
            # SpCas9: NGG (where N = any base, GG required)
            return len(pam_upper) >= 2 and pam_upper[-2:] == "GG"

        elif gRNA_type == GrnaType.SPCAS9N:
            # SpCas9 nickase: NGG or NAG
            return len(pam_upper) >= 2 and pam_upper[-2:] in ("GG", "AG")

        elif gRNA_type == GrnaType.SACAS9:
            # SaCas9: NNGRRT (6bp PAM)
            if len(pam_upper) < 6:
                return False
            r = set("AG")
            t = set("AC")  # T or C
            return pam_upper[2] in r and pam_upper[3] in r and pam_upper[4:6] in ("RT", "RA", "CT", "CA")

        elif gRNA_type == GrnaType.SPCAS9V:
            # SpCas9-HySt: NRY (R=A/G, Y=C/T)
            if len(pam_upper) < 3:
                return False
            return pam_upper[-3] in "ACGT" and pam_upper[-2] in "AG" and pam_upper[-1] in "CT"

        return False

    def reverse_complement(self) -> str:
        """Return reverse complement of the guide sequence."""
        comp = str.maketrans("ACGT", "TGCA")
        return self.sequence[::-1].translate(comp)

    @property
    def gc_content(self) -> float:
        """Calculate GC content of the guide sequence."""
        if not self.sequence:
            return 0.0
        gc = sum(1 for b in self.sequence if b in "GC")
        return gc / len(self.sequence)

    @property
    def gc_percentage(self) -> float:
        """GC content as percentage."""
        return self.gc_content * 100

    @property
    def full_sequence(self) -> str:
        """Full gRNA sequence including PAM."""
        return self.sequence + self.pam_sequence

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "name": self.name or f"{self.gRNA_type.value}_{self.sequence[:8]}",
            "sequence": self.sequence,
            "pam_sequence": self.pam_sequence,
            "gRNA_type": self.gRNA_type.value,
            "chromosome": self.chromosome,
            "position": self.position,
            "strand": self.strand,
            "gc_content": round(self.gc_content, 4),
            "reverse_complement": self.reverse_complement(),
            "full_sequence": self.full_sequence,
        }

    def __str__(self) -> str:
        return (
            f"gRNA({self.name or 'unnamed'}, seq={self.sequence}, "
            f"PAM={self.pam_sequence}, type={self.gRNA_type.value})"
        )


def find_pams_in_sequence(
    sequence: str,
    gRNA_type: GrnaType = GrnaType.SPCAS9,
    guide_length: int = 20,
) -> list[Grna]:
    """Find all PAM sites in a DNA sequence and extract potential gRNAs.

    Args:
        sequence: DNA sequence to scan.
        gRNA_type: Cas enzyme type (determines PAM pattern).
        guide_length: Length of guide sequence upstream of PAM.

    Returns:
        List of Grna objects for all valid PAM sites.
    """
    sequence = sequence.upper()
    guides: list[Grna] = []
    pam_len = PAM_DEFINITIONS[gRNA_type]["length"]

    # Scan for PAMs (search from position guide_length to end)
    for i in range(guide_length, len(sequence) - pam_len + 1):
        pam = sequence[i:i + pam_len]
        if not Grna._validate_pam(pam, gRNA_type):
            continue

        guide_seq = sequence[i - guide_length:i]
        if any(b not in "ACGT" for b in guide_seq):
            continue

        guide = Grna(
            sequence=guide_seq,
            pam_sequence=pam,
            gRNA_type=gRNA_type,
            name=f"pos_{i}_{pam}",
        )
        guides.append(guide)

    return guides


def detect_pam_type(sequence: str) -> dict[str, GrnaType]:
    """Detect which PAM types are present in a sequence.

    Args:
        sequence: DNA sequence to check.

    Returns:
        Dictionary of {gRNA_type: count} for matching PAM types.
    """
    sequence = sequence.upper()
    results: dict[str, list] = {}

    for gtype, info in PAM_DEFINITIONS.items():
        pam_len = info["length"]
        count = 0
        for i in range(len(sequence) - pam_len + 1):
            pam = sequence[i:i + pam_len]
            if Grna._validate_pam(pam, gtype):
                count += 1
        if count > 0:
            results[gtype.value] = count

    return results
