"""Shared test fixtures."""

import pytest
import random
import tempfile
import os
from pathlib import Path


@pytest.fixture
def sample_sequence():
    """A known DNA sequence with GG PAMs."""
    return "AAGCAGTGGTATCAAGTCAGAGGG"


@pytest.fixture
def sample_guide_seq():
    """A 20bp guide sequence."""
    return "AAGCAGTGGTATCAAGTCAG"


@pytest.fixture
def sample_pam():
    """A valid SpCas9 PAM."""
    return "AG"


@pytest.fixture
def test_genome():
    """Create a small test genome for testing."""
    rng = random.Random(42)
    bases = "ACGT"
    genome = {}
    for i in range(1, 4):
        chrom = f"chr{i}"
        seq = "".join(rng.choices(bases, k=10000))
        genome[chrom] = seq
    return genome


@pytest.fixture
def test_fasta(tmp_path):
    """Create a temporary FASTA file."""
    fasta_path = tmp_path / "test_genome.fa"
    sequences = {
        "chr1": "A" * 5000 + "GG" + "AAGCAGTGGTATCAAGTCAG" + "T" * 5000,
        "chr2": "C" * 5000 + "GG" + "TTCGATCGATCGATCGATCG" + "G" * 5000,
    }
    with open(fasta_path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")
    return fasta_path


@pytest.fixture
def valid_guide():
    """Create a valid SpCas9 gRNA."""
    from crispr_design.core.grna import Grna
    return Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")


@pytest.fixture
def valid_guide_with_scores(valid_guide):
    """Guide with pre-computed scores."""
    from crispr_design.core.scoring import compute_all_scores
    return valid_guide, compute_all_scores(valid_guide.sequence)


@pytest.fixture
def off_target_hit():
    """Create a sample off-target hit."""
    from crispr_design.core.off_target import OffTargetHit
    return OffTargetHit(
        chromosome="chr1",
        position=12345,
        strand="+",
        target_sequence="AAGCAGTGGTATCAAGTCAGGG",
        mismatches=1,
        mismatch_positions=[5],
        mismatch_details=[{"position": 5, "query": "G", "target": "T"}],
    )


@pytest.fixture
def design_result(valid_guide):
    """Create a sample design result."""
    from crispr_design.core.scoring import compute_all_scores
    from crispr_design.core.design import DesignResult
    scores = compute_all_scores(valid_guide.sequence)
    return DesignResult(
        guide=valid_guide,
        scores=scores,
        off_targets=[],
        specificity_score=85.0,
        rank=1,
    )


@pytest.fixture
def guide_list_file(tmp_path):
    """Create a file with multiple guide sequences for batch testing."""
    filepath = tmp_path / "guides.txt"
    guides = [
        "AAGCAGTGGTATCAAGTCAG",
        "TTCGATCGATCGATCGATCG",
        "GGCCAATTCCGGCCAATTCC",
        "AAAAAAAACCCCCCCCCCCT",
        "GGGGGGGGTTTTTTTTTTTA",
    ]
    with open(filepath, "w") as f:
        for g in guides:
            f.write(g + "\n")
    return filepath
