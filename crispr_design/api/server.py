"""FastAPI server for CRISPR Design Studio."""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException, Query, File, UploadFile, Form, Request
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from pathlib import Path as PathLib

from ..core.grna import Grna, GrnaType, find_pams_in_sequence
from ..core.design import design_guides, DesignResult
from ..core.off_target import (
    OffTargetHit, find_off_targets, find_off_targets_brute_force,
)
from ..core.scoring import (
    ScoreResult, compute_all_scores, score_mit_cf, score_cfd,
    score_doench_2014, score_doench_2016,
)
from ..databases.genome import Genome, create_test_genome

logger = logging.getLogger(__name__)

app = FastAPI(
    title="CRISPR Design Studio API",
    description="gRNA design and off-target prediction API",
    version="0.1.0",
)

# =============================================================================
# Static UI
# =============================================================================
_web_dir = PathLib(__file__).resolve().parent.parent / "web"
if (_web_dir / "static" / "styles.css").is_file():
    app.mount("/static", StaticFiles(directory=str(_web_dir / "static")), name="static")


@app.get("/", response_class=HTMLResponse)
async def serve_ui():
    """Serve the CRISPR Design Studio web interface."""
    index = _web_dir / "index.html"
    if index.is_file():
        return FileResponse(str(index), media_type="text/html")
    return HTMLResponse(content="<h1>UI not found</h1>", status_code=503)


# =============================================================================
# Pydantic models
# =============================================================================
class GrnaInput(BaseModel):
    """Input for a single gRNA design request."""

    sequence: str = Field(..., description="DNA sequence to scan", min_length=21)
    gRNA_type: str = Field("spcas9", description="Cas enzyme type")
    max_mismatches: int = Field(3, ge=0, le=10, description="Max mismatches")
    max_results: int = Field(50, ge=1, le=1000, description="Max guides to return")
    genome_id: Optional[str] = Field(None, description="Genome identifier")


class BatchInput(BaseModel):
    """Input for batch gRNA processing."""

    sequences: list[str] = Field(..., description="List of gRNA sequences (20bp each)", min_length=1)
    gRNA_type: str = Field("spcas9", description="Cas enzyme type")
    max_mismatches: int = Field(3, ge=0, le=10)
    include_scores: bool = Field(True, description="Include scoring")
    include_off_targets: bool = Field(True, description="Include off-target search")


class ScanInput(BaseModel):
    """Input for off-target scanning."""

    guide_sequence: str = Field(..., description="gRNA sequence (20bp)", min_length=20, max_length=25)
    pam: str = Field("NGG", description="PAM sequence")
    gRNA_type: str = Field("spcas9", description="Cas enzyme type")
    max_mismatches: int = Field(3, ge=0, le=10)
    max_results: int = Field(100, ge=1, le=10000)


class ScoreInput(BaseModel):
    """Input for scoring a gRNA."""

    sequence: str = Field(..., description="gRNA sequence", min_length=1)
    target_sequence: Optional[str] = Field(None, description="Target for CFD scoring")


# =============================================================================
# State
# =============================================================================
_genomes: dict[str, Genome] = {}
_test_genome: Optional[dict[str, str]] = None


def _get_test_genome() -> dict[str, str]:
    """Lazy-load test genome."""
    global _test_genome
    if _test_genome is None:
        logger.info("Creating test genome...")
        _test_genome = create_test_genome(sequence_length=50000)
    return _test_genome


def _get_genome_dict(genome_id: Optional[str] = None) -> dict[str, str]:
    """Get genome dictionary by ID."""
    if not genome_id or genome_id == "test":
        return _get_test_genome()

    if genome_id in _genomes:
        return _genomes[genome_id].to_dict()

    raise HTTPException(status_code=404, detail=f"Genome '{genome_id}' not found")


# =============================================================================
# Routes
# =============================================================================
@app.get("/health")
def health():
    """Health check endpoint."""
    return {"status": "ok", "service": "crispr-design-studio", "version": "0.1.0"}


@app.get("/genomes")
def list_genomes():
    """List available genomes."""
    genomes = {"test": "Synthetic test genome (50kb)"}
    for name in _genomes:
        genomes[name] = f"Custom genome ({len(_genomes[name].chromosomes)} sequences)"
    return {"genomes": genomes}


@app.post("/design")
def design_endpoint(req: GrnaInput):
    """Design gRNAs from a DNA sequence.

    Scans the input sequence for all valid PAM sites, extracts gRNAs,
    scores them, and returns ranked results.
    """
    try:
        gRNA_type = GrnaType(req.gRNA_type.lower())
    except ValueError:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid gRNA type: {req.gRNA_type}. "
                   f"Valid: {[t.value for t in GrnaType]}",
        )

    sequence = req.sequence.upper().strip()
    genome = _get_genome_dict(req.genome_id)

    start = time.time()
    results = design_guides(
        sequence=sequence,
        gRNA_type=gRNA_type,
        genome=genome,
        max_mismatches=req.max_mismatches,
        max_results=req.max_results,
    )
    elapsed = time.time() - start

    return {
        "query": {"sequence_length": len(sequence), "gRNA_type": req.gRNA_type},
        "num_guides": len(results),
        "elapsed_seconds": round(elapsed, 3),
        "guides": [r.to_dict() for r in results],
    }


@app.post("/scan")
def scan_endpoint(req: ScanInput):
    """Scan genome for off-target sites of a gRNA.

    Finds all genomic matches within the specified mismatch tolerance.
    """
    try:
        gRNA_type = GrnaType(req.gRNA_type.lower())
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid gRNA type: {req.gRNA_type}")

    if len(req.guide_sequence) < 20:
        raise HTTPException(status_code=400, detail="Guide sequence must be at least 20bp")

    guide = Grna(
        sequence=req.guide_sequence.upper(),
        pam_sequence=req.pam.upper(),
        gRNA_type=gRNA_type,
    )

    genome = _get_test_genome()

    start = time.time()
    hits = find_off_targets_brute_force(
        guide, genome,
        max_mismatches=req.max_mismatches,
        include_on_target=True,
        max_results=req.max_results,
    )
    elapsed = time.time() - start

    return {
        "guide": guide.to_dict(),
        "num_hits": len(hits),
        "elapsed_seconds": round(elapsed, 3),
        "hits": [h.to_dict() for h in hits],
    }


@app.post("/batch")
def batch_endpoint(req: BatchInput):
    """Batch process multiple gRNAs.

    Scores and optionally scans off-targets for each guide.
    """
    try:
        gRNA_type = GrnaType(req.gRNA_type.lower())
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid gRNA type")

    genome = _get_test_genome()
    results = []

    start = time.time()
    for i, seq in enumerate(req.sequences):
        seq = seq.upper().strip()
        if len(seq) < 20:
            continue

        try:
            guide = Grna(sequence=seq[:20], gRNA_type=gRNA_type, name=f"guide_{i + 1}")
            scores = compute_all_scores(guide.sequence)

            off_targets = []
            if req.include_off_targets:
                off_targets = find_off_targets_brute_force(
                    guide, genome,
                    max_mismatches=req.max_mismatches,
                    include_on_target=False,
                    max_results=50,
                )

            result = DesignResult(
                guide=guide,
                scores=scores,
                off_targets=off_targets,
                specificity_score=max(0, 100 - len(off_targets) * 10),
                rank=i + 1,
            )
            results.append(result)
        except ValueError as e:
            logger.warning(f"Skipped guide {i + 1}: {e}")

    elapsed = time.time() - start

    return {
        "total_guides": len(req.sequences),
        "processed": len(results),
        "elapsed_seconds": round(elapsed, 3),
        "results": [r.to_dict() for r in results],
    }


@app.post("/score")
def score_endpoint(req: ScoreInput):
    """Score a gRNA sequence with all available algorithms."""
    seq = req.sequence.upper().strip()

    if len(seq) != 20:
        raise HTTPException(
            status_code=400,
            detail="Scoring requires exactly 20bp sequence",
        )

    scores = compute_all_scores(seq, req.target_sequence)
    return scores.to_dict()


@app.post("/genomes/upload")
async def upload_genome(
    name: str = Form(..., description="Name for this genome"),
    file: UploadFile = File(..., description="FASTA file"),
):
    """Upload a custom genome FASTA file."""
    import tempfile

    content = await file.read()
    with tempfile.NamedTemporaryFile(mode="wb", suffix=".fa", delete=False) as f:
        f.write(content)
        tmp_path = f.name

    try:
        genome = Genome(tmp_path)
        _genomes[name] = genome
        logger.info(f"Loaded genome '{name}' with {len(genome.chromosomes)} chromosomes")
        return {
            "name": name,
            "chromosomes": len(genome.chromosomes),
            "status": "ok",
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to load genome: {e}")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8001)
