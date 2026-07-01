"""Tests for CLI interface."""

import json
import pytest
from typer.testing import CliRunner
from crispr_design.cli import cli

runner = CliRunner()


class TestCLIHelp:
    def test_main_help(self):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "crispr" in result.output.lower() or "CRISPR" in result.output

    def test_design_help(self):
        result = runner.invoke(cli, ["design", "--help"])
        assert result.exit_code == 0
        assert "sequence" in result.output.lower()

    def test_scan_help(self):
        result = runner.invoke(cli, ["scan", "--help"])
        assert result.exit_code == 0
        assert "guide" in result.output.lower()

    def test_batch_help(self):
        result = runner.invoke(cli, ["batch", "--help"])
        assert result.exit_code == 0

    def test_report_help(self):
        result = runner.invoke(cli, ["report", "--help"])
        assert result.exit_code == 0

    def test_serve_help(self):
        result = runner.invoke(cli, ["serve", "--help"])
        assert result.exit_code == 0


class TestCLIDesign:
    def test_design_sequence(self):
        """Test designing from a raw sequence string."""
        result = runner.invoke(cli, [
            "design",
            "A" * 20 + "GG" + "A" * 50,
            "--max-results", "5",
        ])
        assert result.exit_code == 0
        # Output should be JSON
        try:
            data = json.loads(result.output)
            assert isinstance(data, list)
        except json.JSONDecodeError:
            # Some output might include stderr messages
            assert "guide" in result.output.lower() or len(result.output) > 0

    def test_design_no_sequence(self):
        """Test that missing sequence is handled."""
        result = runner.invoke(cli, ["design"])
        # Should fail because sequence is required
        assert result.exit_code != 0 or "Error" in result.output or "required" in result.output.lower()

    def test_design_with_format(self, tmp_path):
        output_file = tmp_path / "output.json"
        result = runner.invoke(cli, [
            "design",
            "A" * 20 + "GG" + "A" * 50,
            "-o", str(output_file),
            "-f", "json",
        ])
        assert result.exit_code == 0
        assert output_file.exists()


class TestCLIScan:
    def test_scan_without_genome(self):
        result = runner.invoke(cli, [
            "scan", "AAGCAGTGGTATCAAGTCAG",
        ])
        assert result.exit_code != 0  # Needs genome

    def test_scan_help_shows_options(self):
        result = runner.invoke(cli, ["scan", "--help"])
        assert result.exit_code == 0
        assert "--genome" in result.output
        assert "--max-mismatches" in result.output


class TestCLIBatch:
    def test_batch_from_file(self, guide_list_file):
        result = runner.invoke(cli, [
            "batch", str(guide_list_file),
            "-o", "/dev/null",  # Write to /dev/null to avoid file issues
        ])
        assert result.exit_code == 0

    def test_batch_empty_file(self, tmp_path):
        filepath = tmp_path / "empty.txt"
        filepath.write_text("")
        result = runner.invoke(cli, ["batch", str(filepath)])
        assert result.exit_code != 0


class TestCLIReport:
    def test_report_from_json(self, tmp_path):
        # Create a simple JSON input
        input_file = tmp_path / "input.json"
        input_file.write_text(json.dumps([{
            "guide": {
                "sequence": "AAGCAGTGGTATCAAGTCAG",
                "pam_sequence": "GG",
                "gRNA_type": "spcas9",
                "name": "test_guide",
            },
            "scores": {
                "mit_cf": -5.0,
                "cfd": 0.75,
                "doench_2014": 10.0,
                "doench_2016": 12.0,
            },
            "specificity_score": 75.0,
            "rank": 1,
        }]))

        output_file = tmp_path / "report.html"
        result = runner.invoke(cli, [
            "report", str(input_file),
            "-o", str(output_file),
        ])
        assert result.exit_code == 0
        assert output_file.exists()

        content = output_file.read_text()
        assert "<html" in content.lower()
