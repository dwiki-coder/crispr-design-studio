"""Tests for report generation."""

import json
import pytest
from crispr_design.core.design import DesignResult
from crispr_design.core.scoring import ScoreResult, compute_all_scores
from crispr_design.reports.json_report import (
    generate_json_report,
    format_hit_json,
    format_scores_json,
)
from crispr_design.reports.csv_report import (
    generate_csv_report,
    generate_off_target_csv,
)
from crispr_design.reports.html_report import (
    generate_html_report,
    _format_off_targets,
)


class TestJSONReport:
    def test_generate_json(self, design_result, tmp_path):
        output = tmp_path / "report.json"
        path = generate_json_report([design_result], str(output))
        assert path == str(output)

        with open(output) as f:
            data = json.load(f)
        assert data["total_guides"] == 1
        assert len(data["results"]) == 1

    def test_json_with_metadata(self, design_result, tmp_path):
        output = tmp_path / "report.json"
        metadata = {"genome": "hg38", "user": "test"}
        generate_json_report([design_result], str(output), metadata=metadata)

        with open(output) as f:
            data = json.load(f)
        assert data["metadata"]["genome"] == "hg38"

    def test_json_multiple_guides(self, valid_guide, tmp_path):
        scores = compute_all_scores(valid_guide.sequence)
        results = []
        for i in range(3):
            results.append(DesignResult(
                guide=valid_guide,
                scores=scores,
                specificity_score=80.0 - i * 5,
                rank=i + 1,
            ))

        output = tmp_path / "multi.json"
        generate_json_report(results, str(output))

        with open(output) as f:
            data = json.load(f)
        assert data["total_guides"] == 3

    def test_format_hit_json(self, off_target_hit):
        d = format_hit_json(off_target_hit)
        assert d["chromosome"] == "chr1"
        assert d["mismatches"] == 1

    def test_format_scores_json(self):
        scores = ScoreResult(guide_name="test", mit_cf=-5.0, cfd=0.75)
        d = format_scores_json(scores)
        assert d["mit_cf"] == -5.0
        assert d["cfd"] == 0.75


class TestCSVReport:
    def test_generate_csv(self, design_result, tmp_path):
        output = tmp_path / "report.csv"
        path = generate_csv_report([design_result], str(output))
        assert path == str(output)

        with open(output) as f:
            lines = f.readlines()
        assert len(lines) >= 2  # Header + at least 1 data row
        assert "rank" in lines[0]

    def test_csv_header(self, design_result, tmp_path):
        output = tmp_path / "report.csv"
        generate_csv_report([design_result], str(output))

        with open(output) as f:
            header = f.readline().strip()
        assert "sequence" in header
        assert "mit_cf" in header
        assert "doench_2016" in header

    def test_csv_off_target_csv(self, off_target_hit, tmp_path):
        output = tmp_path / "off_targets.csv"
        generate_off_target_csv([off_target_hit], str(output), "guide1")

        with open(output) as f:
            content = f.read()
        assert "chr1" in content
        assert "12345" in content

    def test_csv_multiple_guides(self, valid_guide, tmp_path):
        scores = compute_all_scores(valid_guide.sequence)
        results = []
        for i in range(5):
            results.append(DesignResult(
                guide=valid_guide,
                scores=scores,
                specificity_score=80.0 - i,
                rank=i + 1,
            ))

        output = tmp_path / "multi.csv"
        generate_csv_report(results, str(output))

        with open(output) as f:
            lines = f.readlines()
        assert len(lines) == 6  # Header + 5 data rows


class TestHTMLReport:
    def test_generate_html(self, design_result, tmp_path):
        output = tmp_path / "report.html"
        path = generate_html_report([design_result], str(output))
        assert path == str(output)

        with open(output) as f:
            content = f.read()
        assert "<html" in content.lower()
        assert "CRISPR" in content

    def test_html_with_title(self, design_result, tmp_path):
        output = tmp_path / "report.html"
        generate_html_report(
            [design_result], str(output),
            title="My Custom Report",
        )

        with open(output) as f:
            content = f.read()
        assert "My Custom Report" in content

    def test_html_with_metadata(self, design_result, tmp_path):
        output = tmp_path / "report.html"
        metadata = {"genome": "hg38", "date": "2024-01-01"}
        generate_html_report(
            [design_result], str(output),
            metadata=metadata,
        )

        with open(output) as f:
            content = f.read()
        assert "hg38" in content

    def test_format_off_targets(self, off_target_hit):
        formatted = _format_off_targets([off_target_hit])
        assert len(formatted) == 1
        assert formatted[0]["chromosome"] == "chr1"
        assert formatted[0]["mismatches"] == 1

    def test_format_empty_off_targets(self):
        formatted = _format_off_targets([])
        assert formatted == []

    def test_html_table_structure(self, design_result, tmp_path):
        output = tmp_path / "report.html"
        generate_html_report([design_result], str(output))

        with open(output) as f:
            content = f.read()
        assert "<table" in content.lower()
        assert "<th>" in content.lower()
        assert "<td>" in content.lower()
