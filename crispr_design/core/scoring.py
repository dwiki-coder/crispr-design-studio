"""Scoring algorithms for gRNA on-target and off-target efficiency.

Implements:
- MIT-CF (Cong 2013)
- CFD (Hsu 2013)
- Doench 2014 (Rule Set 1)
- Doench 2016 (Rule Set 2)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

logger = logging.getLogger(__name__)

# =============================================================================
# MIT-CF Scoring Parameters (Cong et al. 2013)
# =============================================================================
# Position weights for SpCas9 (1-based, positions 1-20 from 5' to 3')
MIT_CF_WEIGHTS = [
    0.0,   # pos 1
    -0.15, # pos 2
    -0.27, # pos 3
    -0.13, # pos 4
    -0.42, # pos 5
    -0.29, # pos 6
    0.0,   # pos 7
    -0.09, # pos 8
    0.0,   # pos 9
    -0.35, # pos 10
    -0.39, # pos 11
    -0.24, # pos 12
    -0.11, # pos 13
    0.0,   # pos 14
    0.0,   # pos 15
    0.0,   # pos 16
    0.0,   # pos 17
    0.0,   # pos 18
    0.0,   # pos 19
    -0.11, # pos 20
]

# MIT-CF base-specific penalties at each position
MIT_CF_BASE_PENALTY = {
    1: {"A": -1.11, "T": -1.11, "C": 0.0, "G": 0.0},
    2: {"A": 0.0, "T": -1.02, "C": -0.19, "G": -0.19},
    3: {"A": 0.0, "T": -1.37, "C": -0.43, "G": -0.43},
    4: {"A": 0.0, "T": -1.06, "C": -0.14, "G": -0.14},
    5: {"A": 0.0, "T": -1.50, "C": -0.52, "G": -0.52},
    6: {"A": 0.0, "T": -1.20, "C": -0.29, "G": -0.29},
    7: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    8: {"A": 0.0, "T": -1.00, "C": -0.15, "G": -0.15},
    9: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    10: {"A": 0.0, "T": -1.35, "C": -0.36, "G": -0.36},
    11: {"A": 0.0, "T": -1.39, "C": -0.38, "G": -0.38},
    12: {"A": 0.0, "T": -1.24, "C": -0.22, "G": -0.22},
    13: {"A": 0.0, "T": -1.11, "C": -0.11, "G": -0.11},
    14: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    15: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    16: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    17: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    18: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    19: {"A": 0.0, "T": -1.00, "C": 0.0, "G": 0.0},
    20: {"A": 0.0, "T": -1.11, "C": -0.11, "G": -0.11},
}

# =============================================================================
# CFD Scoring Parameters (Hsu et al. 2013)
# =============================================================================
# Position-dependent base-specific values for mismatches
CFD_PARAMS = {
    1:  {"A": 0.598, "C": 0.604, "G": 0.599, "T": 0.605},
    2:  {"A": 0.805, "C": 0.316, "G": 0.541, "T": 0.339},
    3:  {"A": 0.719, "C": 0.565, "G": 0.832, "T": 0.528},
    4:  {"A": 0.716, "C": 0.511, "G": 0.872, "T": 0.391},
    5:  {"A": 0.456, "C": 0.331, "G": 0.273, "T": 0.514},
    6:  {"A": 0.274, "C": 0.270, "G": 0.451, "T": 0.188},
    7:  {"A": 0.217, "C": 0.209, "G": 0.243, "T": 0.088},
    8:  {"A": 0.184, "C": 0.167, "G": 0.227, "T": 0.160},
    9:  {"A": 0.162, "C": 0.135, "G": 0.120, "T": 0.110},
    10: {"A": 0.102, "C": 0.113, "G": 0.106, "T": 0.100},
    11: {"A": 0.111, "C": 0.126, "G": 0.120, "T": 0.110},
    12: {"A": 0.141, "C": 0.181, "G": 0.174, "T": 0.149},
    13: {"A": 0.209, "C": 0.260, "G": 0.234, "T": 0.236},
    14: {"A": 0.266, "C": 0.274, "G": 0.310, "T": 0.290},
    15: {"A": 0.279, "C": 0.283, "G": 0.306, "T": 0.266},
    16: {"A": 0.271, "C": 0.280, "G": 0.298, "T": 0.265},
    17: {"A": 0.253, "C": 0.266, "G": 0.272, "T": 0.220},
    18: {"A": 0.226, "C": 0.242, "G": 0.260, "T": 0.186},
    19: {"A": 0.202, "C": 0.212, "G": 0.219, "T": 0.160},
    20: {"A": 0.182, "C": 0.187, "G": 0.197, "T": 0.142},
}


# =============================================================================
# Doench 2014 Scoring Parameters (Rule Set 1)
# =============================================================================
DOENCH_2014_PARAMS = {
    # Position -> base value
    1:  {"A": 0.362, "C": 0.380, "G": 0.393, "T": 0.384},
    2:  {"A": 0.811, "C": -0.182, "G": 0.344, "T": -0.170},
    3:  {"A": 0.749, "C": 0.253, "G": 0.977, "T": 0.242},
    4:  {"A": 0.677, "C": 0.183, "G": 1.057, "T": 0.100},
    5:  {"A": 0.487, "C": 0.161, "G": 0.914, "T": 0.159},
    6:  {"A": 0.434, "C": 0.347, "G": 0.779, "T": 0.265},
    7:  {"A": 0.442, "C": 0.344, "G": 0.767, "T": 0.254},
    8:  {"A": 0.500, "C": 0.422, "G": 0.886, "T": 0.344},
    9:  {"A": 0.617, "C": 0.511, "G": 0.989, "T": 0.479},
    10: {"A": 0.631, "C": 0.504, "G": 0.978, "T": 0.482},
    11: {"A": 0.611, "C": 0.480, "G": 0.915, "T": 0.439},
    12: {"A": 0.521, "C": 0.358, "G": 0.806, "T": 0.309},
    13: {"A": 0.358, "C": 0.170, "G": 0.594, "T": 0.117},
    14: {"A": 0.203, "C": 0.086, "G": 0.404, "T": -0.034},
    15: {"A": 0.141, "C": 0.060, "G": 0.295, "T": -0.087},
    16: {"A": 0.132, "C": 0.060, "G": 0.277, "T": -0.100},
    17: {"A": 0.160, "C": 0.095, "G": 0.312, "T": -0.077},
    18: {"A": 0.185, "C": 0.120, "G": 0.330, "T": -0.055},
    19: {"A": 0.208, "C": 0.140, "G": 0.348, "T": -0.033},
    20: {"A": 0.229, "C": 0.157, "G": 0.362, "T": -0.013},
}
DOENCH_2014_OFFSET = -0.495

# =============================================================================
# Doench 2016 Scoring Parameters (Rule Set 2)
# =============================================================================
DOENCH_2016_PARAMS = {
    1:  {"A": 0.365, "C": 0.384, "G": 0.398, "T": 0.386},
    2:  {"A": 0.823, "C": -0.175, "G": 0.352, "T": -0.162},
    3:  {"A": 0.761, "C": 0.261, "G": 0.989, "T": 0.250},
    4:  {"A": 0.689, "C": 0.191, "G": 1.069, "T": 0.108},
    5:  {"A": 0.495, "C": 0.169, "G": 0.926, "T": 0.167},
    6:  {"A": 0.444, "C": 0.356, "G": 0.791, "T": 0.273},
    7:  {"A": 0.450, "C": 0.352, "G": 0.779, "T": 0.262},
    8:  {"A": 0.508, "C": 0.430, "G": 0.898, "T": 0.352},
    9:  {"A": 0.625, "C": 0.519, "G": 1.001, "T": 0.487},
    10: {"A": 0.639, "C": 0.512, "G": 0.990, "T": 0.490},
    11: {"A": 0.619, "C": 0.488, "G": 0.927, "T": 0.447},
    12: {"A": 0.529, "C": 0.366, "G": 0.818, "T": 0.317},
    13: {"A": 0.366, "C": 0.178, "G": 0.606, "T": 0.125},
    14: {"A": 0.211, "C": 0.094, "G": 0.416, "T": -0.026},
    15: {"A": 0.149, "C": 0.068, "G": 0.307, "T": -0.079},
    16: {"A": 0.140, "C": 0.068, "G": 0.289, "T": -0.092},
    17: {"A": 0.168, "C": 0.103, "G": 0.324, "T": -0.069},
    18: {"A": 0.193, "C": 0.128, "G": 0.338, "T": -0.047},
    19: {"A": 0.216, "C": 0.148, "G": 0.356, "T": -0.025},
    20: {"A": 0.237, "C": 0.165, "G": 0.370, "T": -0.005},
}
DOENCH_2016_OFFSET = -0.501


@dataclass
class ScoreResult:
    """Result from scoring a gRNA or off-target site."""

    guide_name: str
    mit_cf: float = 0.0
    cfd: float = 0.0
    doench_2014: float = 0.0
    doench_2016: float = 0.0

    def to_dict(self) -> dict:
        return {
            "guide_name": self.guide_name,
            "mit_cf": round(self.mit_cf, 4),
            "cfd": round(self.cfd, 4),
            "doench_2014": round(self.doench_2014, 4),
            "doench_2016": round(self.doench_2016, 4),
        }


def score_mit_cf(sequence: str) -> float:
    """Compute MIT-CF score (Cong et al. 2013).

    The score is the sum of position-specific penalties for each base.
    Lower scores indicate better specificity.

    Args:
        sequence: 20bp gRNA sequence.

    Returns:
        MIT-CF score (lower = more specific).
    """
    if len(sequence) != 20:
        raise ValueError(f"MIT-CF requires 20bp sequence, got {len(sequence)}")

    score = 0.0
    for i, base in enumerate(sequence):
        pos = i + 1  # 1-based position
        if pos in MIT_CF_BASE_PENALTY:
            score += MIT_CF_BASE_PENALTY[pos].get(base, 0.0)

    return score


def score_cfd(guide_seq: str, target_seq: str) -> float:
    """Compute CFD score (Hsu et al. 2013).

    The CFD score accounts for position-specific mismatch effects.
    Score of 1.0 = perfect match, lower values = more mismatches.

    Args:
        guide_seq: gRNA sequence (20bp).
        target_seq: Target DNA sequence (20bp).

    Returns:
        CFD score between 0 and 1.
    """
    if len(guide_seq) != len(target_seq):
        raise ValueError("Guide and target must be same length")
    if len(guide_seq) > 20:
        guide_seq = guide_seq[:20]
        target_seq = target_seq[:20]

    score = 1.0
    for i in range(min(len(guide_seq), 20)):
        pos = i + 1
        if guide_seq[i] != target_seq[i]:
            # Mismatch penalty: multiply by (1 - parameter value)
            penalty = CFD_PARAMS.get(pos, {}).get(target_seq[i], 0.2)
            score *= (1.0 - penalty)

    return max(0.0, score)


def score_doench_2014(sequence: str) -> float:
    """Compute Doench 2014 rule set 1 score.

    Predicts on-target activity for SpCas9 guides.

    Args:
        sequence: 20bp gRNA sequence.

    Returns:
        Activity score (higher = more active).
    """
    if len(sequence) != 20:
        raise ValueError(f"Doench 2014 requires 20bp sequence, got {len(sequence)}")

    score = DOENCH_2014_OFFSET
    for i, base in enumerate(sequence):
        pos = i + 1
        if pos in DOENCH_2014_PARAMS:
            score += DOENCH_2014_PARAMS[pos].get(base, 0.0)

    return score


def score_doench_2016(sequence: str) -> float:
    """Compute Doench 2016 rule set 2 score.

    Improved prediction of on-target activity for SpCas9 guides.

    Args:
        sequence: 20bp gRNA sequence.

    Returns:
        Activity score (higher = more active).
    """
    if len(sequence) != 20:
        raise ValueError(f"Doench 2016 requires 20bp sequence, got {len(sequence)}")

    score = DOENCH_2016_OFFSET
    for i, base in enumerate(sequence):
        pos = i + 1
        if pos in DOENCH_2016_PARAMS:
            score += DOENCH_2016_PARAMS[pos].get(base, 0.0)

    return score


def compute_all_scores(sequence: str, target_seq: Optional[str] = None) -> ScoreResult:
    """Compute all scoring algorithms for a gRNA.

    Args:
        sequence: 20bp gRNA sequence.
        target_seq: Optional target sequence for CFD (uses guide_seq if None).

    Returns:
        ScoreResult with all scores.
    """
    if target_seq is None:
        target_seq = sequence

    name = f"{sequence[:8]}"
    result = ScoreResult(guide_name=name)

    try:
        result.mit_cf = score_mit_cf(sequence)
    except ValueError:
        pass

    try:
        result.cfd = score_cfd(sequence, target_seq)
    except ValueError:
        pass

    try:
        result.doench_2014 = score_doench_2014(sequence)
    except ValueError:
        pass

    try:
        result.doench_2016 = score_doench_2016(sequence)
    except ValueError:
        pass

    return result
