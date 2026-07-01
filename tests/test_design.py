"""Tests for guide design module."""

import pytest
from crispr_design.core.grna import Grna, GrnaType
from crispr_design.core.design import (
    DesignResult,
    design_guides,
    filter_guides,
    _compute_specificity,
)
from crispr_design.core.scoring import ScoreResult, compute_all_scores


class TestDesignResult:
    def test_create_design_result(self, valid_guide):
        scores = compute_all_scores(valid_guide.sequence)
        result = DesignResult(
            guide=valid_guide,
            scores=scores,
            off_targets=[],
            specificity_score=85.0,
            rank=1,
        )
        assert result.rank == 1
        assert result.specificity_score == 85.0

    def test_design_result_to_dict(self, design_result):
        d = design_result.to_dict()
        assert d["rank"] == 1
        assert "guide" in d
        assert "scores" in d
        assert d["off_target_count"] == 0

    def test_design_result_with_off_targets(self, valid_guide, off_target_hit):
        scores = compute_all_scores(valid_guide.sequence)
        result = DesignResult(
            guide=valid_guide,
            scores=scores,
            off_targets=[off_target_hit],
            specificity_score=70.0,
            rank=1,
        )
        assert result.off_target_count == 1
        d = result.to_dict()
        assert len(d["off_targets"]) == 1


class TestComputeSpecificity:
    def test_no_off_targets(self):
        scores = ScoreResult(guide_name="test", doench_2016=20.0)
        specificity = _compute_specificity(scores, [], 3)
        assert 0 <= specificity <= 100

    def test_with_off_targets(self, off_target_hit):
        scores = ScoreResult(guide_name="test", doench_2016=20.0)
        specificity = _compute_specificity(scores, [off_target_hit], 3)
        assert specificity >= 0

    def test_many_off_targets_lower_score(self, off_target_hit):
        scores = ScoreResult(guide_name="test", doench_2016=20.0)
        few_ot = _compute_specificity(scores, [off_target_hit], 3)
        many_ot = _compute_specificity(scores, [off_target_hit] * 5, 3)
        assert many_ot <= few_ot


class TestDesignGuides:
    def test_design_from_sequence(self):
        """Design guides from a known sequence with GG PAMs."""
        seq = "A" * 20 + "GG" + "A" * 50
        results = design_guides(
            sequence=seq,
            gRNA_type=GrnaType.SPCAS9,
            max_results=10,
        )
        assert len(results) >= 1
        assert all(isinstance(r, DesignResult) for r in results)

    def test_design_returns_ranked(self):
        seq = "A" * 100 + "GG" + "A" * 50 + "GG" + "A" * 50
        results = design_guides(
            sequence=seq,
            gRNA_type=GrnaType.SPCAS9,
            max_results=10,
        )
        if len(results) > 1:
            # Results should be ranked
            assert all(r.rank > 0 for r in results)

    def test_design_with_genome(self, test_genome):
        seq = "A" * 20 + "GG" + "A" * 50
        results = design_guides(
            sequence=seq,
            gRNA_type=GrnaType.SPCAS9,
            genome=test_genome,
            max_mismatches=3,
            max_results=5,
        )
        assert len(results) >= 1
        # Check that off-targets were found
        for r in results:
            assert isinstance(r.off_targets, list)

    def test_design_max_results_limit(self):
        seq = "A" * 200 + "GG" + "A" * 200
        results = design_guides(
            sequence=seq,
            gRNA_type=GrnaType.SPCAS9,
            max_results=3,
        )
        assert len(results) <= 3

    def test_design_min_specificity_filter(self):
        seq = "A" * 20 + "GG" + "A" * 50
        results = design_guides(
            sequence=seq,
            gRNA_type=GrnaType.SPCAS9,
            min_specificity=0.0,
            max_results=10,
        )
        assert isinstance(results, list)

    def test_design_no_pam(self):
        seq = "A" * 100  # No GG PAM
        results = design_guides(
            sequence=seq,
            gRNA_type=GrnaType.SPCAS9,
        )
        assert len(results) == 0

    def test_design_with_guide_prefix(self):
        seq = "A" * 20 + "GG" + "A" * 50
        results = design_guides(
            sequence=seq,
            gRNA_type=GrnaType.SPCAS9,
            guide_name_prefix="MY",
            max_results=10,
        )
        if results:
            assert "MY" in results[0].guide.name


class TestFilterGuides:
    def test_filter_gc_content(self):
        results = []
        for i in range(10):
            seq = "G" * (5 + i) + "A" * (15 - i)  # Varying GC
            try:
                guide = Grna(sequence=seq, pam_sequence="GG")
                scores = compute_all_scores(guide.sequence)
                results.append(DesignResult(
                    guide=guide,
                    scores=scores,
                    specificity_score=50.0,
                    rank=i + 1,
                ))
            except ValueError:
                pass

        filtered = filter_guides(results, min_gc=0.3, max_gc=0.7)
        for r in filtered:
            assert 0.3 <= r.guide.gc_content <= 0.7

    def test_filter_max_off_targets(self, off_target_hit):
        results = []
        for i in range(3):
            try:
                seq = "G" * 5 + "A" * 15
                guide = Grna(sequence=seq, pam_sequence="GG")
                scores = compute_all_scores(guide.sequence)
                num_ot = [off_target_hit] * i
                results.append(DesignResult(
                    guide=guide,
                    scores=scores,
                    off_targets=num_ot,
                    specificity_score=50.0,
                    rank=i + 1,
                ))
            except ValueError:
                pass

        filtered = filter_guides(results, max_off_targets=1)
        for r in filtered:
            assert len(r.off_targets) <= 1

    def test_filter_empty_results(self):
        filtered = filter_guides([])
        assert filtered == []
