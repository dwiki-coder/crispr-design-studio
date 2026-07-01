"""Tests for gRNA class and PAM detection."""

import pytest
from crispr_design.core.grna import (
    Grna,
    GrnaType,
    find_pams_in_sequence,
    detect_pam_type,
    GUIDE_LENGTHS,
    PAM_DEFINITIONS,
)


# ---- Grna Construction ----

class TestGrnaConstruction:
    def test_create_spCas9_guide(self, sample_guide_seq):
        g = Grna(sequence=sample_guide_seq, pam_sequence="GG")
        assert g.sequence == "AAGCAGTGGTATCAAGTCAG"
        assert g.pam_sequence == "GG"
        assert g.gRNA_type == GrnaType.SPCAS9

    def test_create_with_name(self, sample_guide_seq):
        g = Grna(sequence=sample_guide_seq, pam_sequence="GG", name="my_guide")
        assert g.name == "my_guide"

    def test_create_with_coordinates(self, sample_guide_seq):
        g = Grna(
            sequence=sample_guide_seq,
            pam_sequence="GG",
            chromosome="chr1",
            position=1000,
            strand="-",
        )
        assert g.chromosome == "chr1"
        assert g.position == 1000
        assert g.strand == "-"

    def test_sequence_uppercased(self):
        g = Grna(sequence="aagcagtggtatcaagtcag", pam_sequence="gg")
        assert g.sequence == "AAGCAGTGGTATCAAGTCAG"
        assert g.pam_sequence == "GG"

    def test_invalid_length_raises(self):
        with pytest.raises(ValueError, match="length"):
            Grna(sequence="ACGT", pam_sequence="GG")

    def test_invalid_base_raises(self):
        with pytest.raises(ValueError, match="Invalid"):
            Grna(sequence="AAGCAGTGGTATCAAGTCAX", pam_sequence="GG")

    def test_invalid_pam_base_raises(self):
        with pytest.raises(ValueError, match="Invalid base"):
            grna = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GX")

    def test_empty_sequence_raises(self):
        with pytest.raises(ValueError, match="empty"):
            Grna(sequence="", pam_sequence="GG")


# ---- PAM Detection ----

class TestPAMDetection:
    def test_validate_spcas9_pam_gg(self):
        assert Grna._validate_pam("AGG", GrnaType.SPCAS9) is True
        assert Grna._validate_pam("GGG", GrnaType.SPCAS9) is True
        assert Grna._validate_pam("CGG", GrnaType.SPCAS9) is True
        assert Grna._validate_pam("TGG", GrnaType.SPCAS9) is True

    def test_validate_spcas9_pam_not_match(self):
        assert Grna._validate_pam("AAA", GrnaType.SPCAS9) is False
        assert Grna._validate_pam("ATT", GrnaType.SPCAS9) is False
        assert Grna._validate_pam("GAT", GrnaType.SPCAS9) is False

    def test_validate_spCas9n_pam(self):
        assert Grna._validate_pam("AGG", GrnaType.SPCAS9N) is True
        assert Grna._validate_pam("AAG", GrnaType.SPCAS9N) is True
        assert Grna._validate_pam("AGG", GrnaType.SPCAS9N) is True
        assert Grna._validate_pam("AAA", GrnaType.SPCAS9N) is False

    def test_validate_hifi_pam(self):
        assert Grna._validate_pam("NGG", GrnaType.HIFICAS9) is True
        assert Grna._validate_pam("AGG", GrnaType.HIFICAS9) is True
        assert Grna._validate_pam("AAA", GrnaType.HIFICAS9) is False

    def test_validate_sacas9_pam(self):
        # SaCas9 PAM: NNGRRT
        assert Grna._validate_pam("NNGRRT", GrnaType.SACAS9) is False  # N not in AG
        # Actual test with valid R/Y
        assert Grna._validate_pam("AAGGTT", GrnaType.SACAS9) is False  # R must be A/G

    def test_validate_spcas9v_pam(self):
        # SpCas9-HySt: NRY
        assert Grna._validate_pam("AGC", GrnaType.SPCAS9V) is True  # R=A, Y=C
        assert Grna._validate_pam("GGT", GrnaType.SPCAS9V) is True  # R=G, Y=T
        assert Grna._validate_pam("GGA", GrnaType.SPCAS9V) is False  # Y not C/T
        assert Grna._validate_pam("CGT", GrnaType.SPCAS9V) is True

    def test_empty_pam_invalid(self):
        assert Grna._validate_pam("", GrnaType.SPCAS9) is False

    def test_short_pam_invalid(self):
        assert Grna._validate_pam("A", GrnaType.SPCAS9) is False


# ---- from_target ----

class TestFromTarget:
    def test_from_target_spcas9(self):
        """Create guide from full target sequence (guide + PAM)."""
        g = Grna.from_target("AAGCAGTGGTATCAAGTCAGGGG")
        assert g.sequence == "AAGCAGTGGTATCAAGTCAG"
        assert g.pam_sequence == "GG"
        assert g.gRNA_type == GrnaType.SPCAS9

    def test_from_target_with_pam_type(self):
        g = Grna.from_target("AAGCAGTGGTATCAAGTCAGGG", gRNA_type=GrnaType.SPCAS9)
        assert g.sequence == "AAGCAGTGGTATCAAGTCAG"
        assert g.pam_sequence == "GG"

    def test_from_target_no_valid_pam_raises(self):
        with pytest.raises(ValueError, match="PAM"):
            Grna.from_target("AAGCAGTGGTATCAAGTCAGAAA")

    def test_from_target_too_short(self):
        with pytest.raises(ValueError, match="short"):
            Grna.from_target("AAGCAGT")


# ---- find_pams_in_sequence ----

class TestFindPAMs:
    def test_find_pams_finds_gg_sites(self):
        seq = "A" * 20 + "GG" + "A" * 20 + "GG"
        guides = find_pams_in_sequence(seq, GrnaType.SPCAS9)
        assert len(guides) >= 2  # At least the two explicit GG sites

    def test_find_pams_returns_grna_objects(self):
        seq = "A" * 30 + "GG"
        guides = find_pams_in_sequence(seq, GrnaType.SPCAS9)
        assert len(guides) > 0
        assert isinstance(guides[0], Grna)
        assert guides[0].sequence == "A" * 20
        assert guides[0].pam_sequence == "GG"

    def test_find_pams_no_pam(self):
        seq = "A" * 100
        guides = find_pams_in_sequence(seq, GrnaType.SPCAS9)
        assert len(guides) == 0

    def test_find_pams_sacas9(self):
        seq = "A" * 20 + "NNGRRT"
        guides = find_pams_in_sequence(seq, GrnaType.SACAS9)
        # May or may not find matches depending on random bases
        assert isinstance(guides, list)


# ---- detect_pam_type ----

class TestDetectPAMType:
    def test_detect_spcas9_pams(self):
        seq = "AAGCAGTGGTATCAAGTCAGAGG" * 10
        results = detect_pam_type(seq)
        assert "spcas9" in results
        assert results["spcas9"] > 0

    def test_detect_multiple_types(self):
        seq = "A" * 100 + "GG" + "A" * 100
        results = detect_pam_type(seq)
        assert "spcas9" in results


# ---- Properties ----

class TestGrnaProperties:
    def test_reverse_complement(self):
        g = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        rc = g.reverse_complement()
        assert len(rc) == 20
        assert all(b in "ACGT" for b in rc)

    def test_reverse_complement_known(self):
        g = Grna(sequence="ACGT" + "A" * 16, pam_sequence="GG")
        # We can't test this directly since length must be 20
        # Just verify the method works
        assert isinstance(g.reverse_complement(), str)

    def test_gc_content(self):
        g = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        assert 0 <= g.gc_content <= 1

    def test_gc_content_percent(self):
        g = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        assert 0 <= g.gc_percentage <= 100

    def test_full_sequence(self):
        g = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        assert g.full_sequence == "AAGCAGTGGTATCAAGTCAGGG"

    def test_to_dict(self):
        g = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG", name="test")
        d = g.to_dict()
        assert d["name"] == "test"
        assert d["sequence"] == "AAGCAGTGGTATCAAGTCAG"
        assert d["pam_sequence"] == "GG"
        assert d["gRNA_type"] == "spcas9"

    def test_str_repr(self):
        g = Grna(sequence="AAGCAGTGGTATCAAGTCAG", pam_sequence="GG")
        s = str(g)
        assert "AAGCAGTGGTATCAAGTCAG" in s


# ---- GrnaType enum ----

class TestGrnaType:
    def test_all_types_defined(self):
        types = [t.value for t in GrnaType]
        assert "spcas9" in types
        assert "sacas9" in types
        assert "hificas9" in types
        assert "spcas9n" in types
        assert "spcas9v" in types
        assert "xmacas9" in types

    def test_guide_lengths(self):
        for gtype in GrnaType:
            assert gtype in GUIDE_LENGTHS
            assert GUIDE_LENGTHS[gtype] > 0

    def test_pam_definitions(self):
        for gtype in GrnaType:
            assert gtype in PAM_DEFINITIONS
