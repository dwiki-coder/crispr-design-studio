"""Tests for sequence utilities."""

import pytest
from crispr_design.utils.seq_utils import (
    reverse_complement,
    validate_dna_sequence,
    count_gc_content,
    hamming_distance,
    parse_fasta,
    extract_grna_from_fasta,
    format_sequence,
)


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ACGT") == "ACGT"

    def test_atrich(self):
        assert reverse_complement("AAAA") == "TTTT"

    def test_gcrich(self):
        assert reverse_complement("GGGG") == "CCCC"

    def test_mixed(self):
        assert reverse_complement("ACGTACGT") == "ACGTACGT"

    def test_lowercase_input(self):
        assert reverse_complement("acgt") == "ACGT"

    def test_empty(self):
        assert reverse_complement("") == ""


class TestValidateDNASequence:
    def test_valid_sequence(self):
        assert validate_dna_sequence("ACGTACGT") is True

    def test_lowercase_valid(self):
        assert validate_dna_sequence("acgt") is True

    def test_invalid_base(self):
        assert validate_dna_sequence("ACGTN") is False

    def test_with_n_allowed(self):
        assert validate_dna_sequence("ACGTN", allow_n=True) is True

    def test_empty_string(self):
        assert validate_dna_sequence("") is True

    def test_whitespace(self):
        assert validate_dna_sequence("  ACGT  ") is True

    def test_number_invalid(self):
        assert validate_dna_sequence("AC1GT") is False


class TestGCContent:
    def test_all_gc(self):
        assert count_gc_content("GCGC") == 1.0

    def test_all_at(self):
        assert count_gc_content("ATAT") == 0.0

    def test_half_gc(self):
        assert count_gc_content("ACGT") == 0.5

    def test_empty(self):
        assert count_gc_content("") == 0.0

    def test_single_g(self):
        assert count_gc_content("G") == 1.0

    def test_mixed(self):
        gc = count_gc_content("AAGCAGTGGTATCAAGTCAG")
        assert 0 <= gc <= 1


class TestHammingDistance:
    def test_same_sequence(self):
        assert hamming_distance("ACGT", "ACGT") == 0

    def test_single_difference(self):
        assert hamming_distance("ACGT", "ACGA") == 1

    def test_all_different(self):
        assert hamming_distance("AAAA", "TTTT") == 4

    def test_different_lengths_raises(self):
        with pytest.raises(ValueError):
            hamming_distance("ACGT", "ACG")


class TestFormatSequence:
    def test_basic_format(self):
        seq = "A" * 120
        result = format_sequence(seq, width=60)
        lines = result.split("\n")
        assert len(lines) == 2

    def test_with_positions(self):
        seq = "A" * 60
        result = format_sequence(seq, width=60, show_positions=True)
        assert "      1" in result or "  1" in result or "1" in result

    def test_no_positions(self):
        seq = "A" * 60
        result = format_sequence(seq, width=60, show_positions=False)
        assert result == "A" * 60

    def test_multiline(self):
        seq = "A" * 100
        result = format_sequence(seq, width=60, show_positions=False)
        lines = result.split("\n")
        assert len(lines) == 2


class TestParseFasta:
    def test_parse_fasta(self, test_fasta):
        sequences = parse_fasta(str(test_fasta))
        assert len(sequences) == 2
        assert "chr1" in sequences
        assert "chr2" in sequences

    def test_fasta_uppercase(self, test_fasta):
        sequences = parse_fasta(str(test_fasta))
        for seq in sequences.values():
            assert seq == seq.upper()


class TestExtractGrnaFromFasta:
    def test_extract_spcas9(self, test_fasta):
        guides = extract_grna_from_fasta(
            str(test_fasta),
            guide_length=20,
            pam_pattern="NGG",
        )
        assert isinstance(guides, list)

    def test_extract_returns_dicts(self, test_fasta):
        guides = extract_grna_from_fasta(
            str(test_fasta),
            guide_length=20,
            pam_pattern="NGG",
        )
        if guides:
            assert "sequence" in guides[0]
            assert "pam" in guides[0]
            assert "position" in guides[0]
