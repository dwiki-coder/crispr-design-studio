"""Tests for genome handling."""

import pytest
import tempfile
from crispr_design.databases.genome import Genome, create_test_genome


class TestCreateTestGenome:
    def test_creates_genome(self):
        genome = create_test_genome(sequence_length=1000)
        assert isinstance(genome, dict)
        assert len(genome) > 0

    def test_chromosome_names(self):
        genome = create_test_genome()
        assert "chr1" in genome
        assert "chr2" in genome

    def test_sequence_length(self):
        genome = create_test_genome(sequence_length=1000)
        for chrom, seq in genome.items():
            assert len(seq) >= 1000

    def test_valid_bases(self):
        genome = create_test_genome()
        for chrom, seq in genome.items():
            assert all(b in "ACGT" for b in seq)

    def test_reproducible_with_seed(self):
        g1 = create_test_genome(sequence_length=1000, seed=42)
        g2 = create_test_genome(sequence_length=1000, seed=42)
        for chrom in g1:
            assert g1[chrom] == g2[chrom]


class TestGenomeFromFasta:
    def test_load_fasta(self, test_fasta):
        genome = Genome(test_fasta)
        assert genome is not None

    def test_get_sequence(self, test_fasta):
        genome = Genome(test_fasta)
        seq = genome.get_sequence("chr1", 0, 10)
        assert len(seq) == 10
        assert all(b in "ACGT" for b in seq)

    def test_get_full_chromosome(self, test_fasta):
        genome = Genome(test_fasta)
        seq = genome.get_sequence("chr1")
        assert len(seq) > 0

    def test_reverse_complement(self, test_fasta):
        genome = Genome(test_fasta)
        fwd = genome.get_sequence("chr1", 0, 10, strand="+")
        rev = genome.get_sequence("chr1", 0, 10, strand="-")
        assert len(fwd) == len(rev)

    def test_chromosomes_property(self, test_fasta):
        genome = Genome(test_fasta)
        chroms = genome.chromosomes
        assert "chr1" in chroms
        assert "chr2" in chroms

    def test_chromosome_sizes(self, test_fasta):
        genome = Genome(test_fasta)
        sizes = genome.chromosome_sizes
        assert len(sizes) > 0
        for chrom, size in sizes.items():
            assert size > 0

    def test_region_context(self, test_fasta):
        genome = Genome(test_fasta)
        context = genome.get_region_context("chr1", 5000, flank_size=10)
        assert "sequence" in context
        assert context["chromosome"] == "chr1"
        assert context["position"] == 5000

    def test_to_dict(self, test_fasta):
        genome = Genome(test_fasta)
        d = genome.to_dict()
        assert len(d) > 0
        for chrom, seq in d.items():
            assert len(seq) > 0


class TestGenomeErrors:
    def test_nonexistent_fasta(self):
        with pytest.raises(FileNotFoundError):
            Genome("/nonexistent/path.fasta")

    def test_unknown_chromosome(self, test_fasta):
        genome = Genome(test_fasta)
        with pytest.raises(KeyError):
            genome.get_sequence("chr999")
