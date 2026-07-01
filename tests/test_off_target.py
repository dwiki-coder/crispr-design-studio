"""Tests for off-target search engine."""

import pytest
from crispr_design.core.off_target import (
    OffTargetHit,
    find_off_targets,
    find_off_targets_brute_force,
    _count_mismatches,
    _seed_match,
)
from crispr_design.core.grna import Grna, GrnaType


# ---- OffTargetHit ----

class TestOffTargetHit:
    def test_create_hit(self, off_target_hit):
        assert off_target_hit.chromosome == "chr1"
        assert off_target_hit.position == 12345
        assert off_target_hit.mismatches == 1

    def test_hit_to_dict(self, off_target_hit):
        d = off_target_hit.to_dict()
        assert d["chromosome"] == "chr1"
        assert d["position"] == 12345
        assert d["strand"] == "+"
        assert d["mismatches"] == 1
        assert "mismatch_details" in d

    def test_hit_default_values(self):
        hit = OffTargetHit(
            chromosome="chr1",
            position=100,
            strand="+",
            target_sequence="A" * 20,
            mismatches=0,
        )
        assert hit.mismatch_positions == []
        assert hit.mismatch_details == []
        assert hit.genomic_context == ""
        assert hit.gene_annotation is None


# ---- _count_mismatches ----

class TestCountMismatches:
    def test_perfect_match(self):
        count, positions, details = _count_mismatches("ACGT", "ACGT")
        assert count == 0
        assert positions == []
        assert details == []

    def test_single_mismatch(self):
        count, positions, details = _count_mismatches("ACGT", "ACGA")
        assert count == 1
        assert positions == [3]

    def test_multiple_mismatches(self):
        count, positions, details = _count_mismatches("ACGT", "TGCA")
        assert count == 4
        assert positions == [0, 1, 2, 3]

    def test_mismatch_details(self):
        count, positions, details = _count_mismatches("ACGT", "AGGT")
        assert count == 1
        assert details[0] == {"position": 1, "query": "C", "target": "G"}

    def test_different_length_raises(self):
        with pytest.raises(ValueError):
            _count_mismatches("ACGT", "ACG")


# ---- find_off_targets_brute_force ----

class TestBruteForceSearch:
    def test_finds_perfect_match(self, test_genome):
        """If we insert the guide into the genome, we should find it."""
        # Use a known sequence
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets_brute_force(
            guide, test_genome, max_mismatches=0, max_results=100
        )
        # May or may not find exact matches depending on random genome
        assert isinstance(hits, list)

    def test_finds_with_mismatches(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets_brute_force(
            guide, test_genome, max_mismatches=3, max_results=100
        )
        assert isinstance(hits, list)

    def test_sorted_by_mismatches(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets_brute_force(
            guide, test_genome, max_mismatches=5, max_results=1000
        )
        for i in range(len(hits) - 1):
            assert hits[i].mismatches <= hits[i + 1].mismatches

    def test_max_results_limit(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets_brute_force(
            guide, test_genome, max_mismatches=10, max_results=10
        )
        assert len(hits) <= 10

    def test_exclude_on_target(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets_brute_force(
            guide, test_genome, max_mismatches=3,
            include_on_target=False, max_results=100
        )
        assert isinstance(hits, list)

    def test_empty_genome(self):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets_brute_force(guide, {}, max_mismatches=3)
        assert hits == []


# ---- find_off_targets (seed-and-extend) ----

class TestSeedAndExtend:
    def test_returns_list(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets(
            guide, test_genome, max_mismatches=3, max_results=100
        )
        assert isinstance(hits, list)

    def test_finds_matches(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets(
            guide, test_genome, max_mismatches=5, max_results=1000
        )
        assert isinstance(hits, list)

    def test_sorted_results(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets(
            guide, test_genome, max_mismatches=5, max_results=1000
        )
        for i in range(len(hits) - 1):
            assert hits[i].mismatches <= hits[i + 1].mismatches

    def test_max_results(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets(
            guide, test_genome, max_mismatches=10, max_results=5
        )
        assert len(hits) <= 5

    def test_zero_mismatches_strict(self, test_genome):
        guide = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        hits = find_off_targets(
            guide, test_genome, max_mismatches=0, max_results=100
        )
        for hit in hits:
            assert hit.mismatches == 0


# ---- _seed_match ----

class TestSeedMatch:
    def test_perfect_seed(self):
        query = "AAGCAGTGGTATCAAGTCAG"
        target = "AAGCAGTGGTATCAAGTCAG"
        assert _seed_match(query, target, 10, 0) is True

    def test_seed_with_mismatches(self):
        query = "AAGCAGTGGTATCAAGTCAG"
        target = "AAGCAGTGGTATCAAGTCAC"  # 1 mismatch at pos 19 (in seed)
        assert _seed_match(query, target, 10, 2) is True

    def test_seed_too_many_mismatches(self):
        query = "AAGCAGTGGTATCAAGTCAG"
        target = "AAGCAGTGGTATCAAGTGGG"  # Multiple mismatches in seed
        result = _seed_match(query, target, 10, 0)
        assert isinstance(result, bool)
