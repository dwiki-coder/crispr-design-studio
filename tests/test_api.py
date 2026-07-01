"""Tests for FastAPI server."""

import pytest
from fastapi.testclient import TestClient
from crispr_design.api.server import app


@pytest.fixture
def client():
    return TestClient(app)


class TestHealthEndpoint:
    def test_health_check(self, client):
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "ok"
        assert "crispr-design" in data["service"]

    def test_health_has_version(self, client):
        response = client.get("/health")
        data = response.json()
        assert "version" in data


class TestGenomesEndpoint:
    def test_list_genomes(self, client):
        response = client.get("/genomes")
        assert response.status_code == 200
        data = response.json()
        assert "genomes" in data
        assert "test" in data["genomes"]


class TestDesignEndpoint:
    def test_basic_design(self, client):
        response = client.post("/design", json={
            "sequence": "A" * 20 + "GG" + "A" * 50,
            "gRNA_type": "spcas9",
            "max_mismatches": 3,
            "max_results": 10,
        })
        assert response.status_code == 200
        data = response.json()
        assert "guides" in data
        assert "num_guides" in data
        assert "elapsed_seconds" in data

    def test_design_with_specific_type(self, client):
        response = client.post("/design", json={
            "sequence": "A" * 20 + "GG" + "A" * 50,
            "gRNA_type": "hificas9",
            "max_results": 5,
        })
        assert response.status_code == 200

    def test_design_invalid_type(self, client):
        response = client.post("/design", json={
            "sequence": "A" * 20 + "GG" + "A" * 50,
            "gRNA_type": "invalid_type",
        })
        assert response.status_code == 400

    def test_design_returns_guides(self, client):
        response = client.post("/design", json={
            "sequence": "A" * 20 + "GG" + "A" * 50 + "GG" + "A" * 50,
            "gRNA_type": "spcas9",
            "max_results": 5,
        })
        data = response.json()
        assert data["num_guides"] <= 5
        for guide in data["guides"]:
            assert "guide" in guide
            assert "scores" in guide


class TestScanEndpoint:
    def test_basic_scan(self, client):
        response = client.post("/scan", json={
            "guide_sequence": "AAGCAGTGGTATCAAGTCAG",
            "pam": "GG",
            "gRNA_type": "spcas9",
            "max_mismatches": 3,
            "max_results": 100,
        })
        assert response.status_code == 200
        data = response.json()
        assert "hits" in data
        assert "num_hits" in data
        assert "guide" in data

    def test_scan_short_sequence(self, client):
        response = client.post("/scan", json={
            "guide_sequence": "ACGT",
            "pam": "GG",
        })
        assert response.status_code in (400, 422)

    def test_scan_invalid_type(self, client):
        response = client.post("/scan", json={
            "guide_sequence": "AAGCAGTGGTATCAAGTCAG",
            "pam": "GG",
            "gRNA_type": "invalid",
        })
        assert response.status_code == 400

    def test_scan_returns_hits(self, client):
        response = client.post("/scan", json={
            "guide_sequence": "AAGCAGTGGTATCAAGTCAG",
            "pam": "GG",
            "max_mismatches": 5,
            "max_results": 10,
        })
        data = response.json()
        for hit in data["hits"][:5]:
            assert "chromosome" in hit
            assert "position" in hit
            assert "mismatches" in hit


class TestBatchEndpoint:
    def test_basic_batch(self, client):
        sequences = [
            "AAGCAGTGGTATCAAGTCAG",
            "TTCGATCGATCGATCGATCG",
            "GGCCAATTCCGGCCAATTCC",
        ]
        response = client.post("/batch", json={
            "sequences": sequences,
            "gRNA_type": "spcas9",
            "max_mismatches": 3,
        })
        assert response.status_code == 200
        data = response.json()
        assert "results" in data
        assert data["processed"] == 3

    def test_batch_empty(self, client):
        response = client.post("/batch", json={
            "sequences": [],
        })
        assert response.status_code == 422  # Validation error

    def test_batch_with_scores(self, client):
        response = client.post("/batch", json={
            "sequences": ["AAGCAGTGGTATCAAGTCAG"],
            "include_scores": True,
            "include_off_targets": False,
        })
        data = response.json()
        for r in data["results"]:
            assert "scores" in r


class TestScoreEndpoint:
    def test_basic_score(self, client):
        response = client.post("/score", json={
            "sequence": "AAGCAGTGGTATCAAGTCAG",
        })
        assert response.status_code == 200
        data = response.json()
        assert "mit_cf" in data
        assert "cfd" in data
        assert "doench_2014" in data
        assert "doench_2016" in data

    def test_score_with_target(self, client):
        response = client.post("/score", json={
            "sequence": "AAGCAGTGGTATCAAGTCAG",
            "target_sequence": "AAGCAGTGGTATCAAGTCAC",
        })
        data = response.json()
        assert data["cfd"] < 1.0

    def test_score_wrong_length(self, client):
        response = client.post("/score", json={
            "sequence": "ACGT",
        })
        assert response.status_code == 400
