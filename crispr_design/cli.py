"""CRISPR Design Studio CLI — main entry point."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import typer

app = typer.Typer(
    name="crispr-design",
    help="CRISPR Design Studio — gRNA design and off-target prediction",
    add_completion=False,
)

# Alias for backward compat
cli = app


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")


# =============================================================================
# design subcommand
# =============================================================================
@cli.command()
def design(
    sequence: str = typer.Argument(
        None, help="DNA sequence to scan for gRNAs (or path to FASTA file)"
    ),
    gRNA_type: str = typer.Option(
        "spcas9", "--type", "-t",
        help="Cas enzyme type (spcas9, sacas9, hificas9, spcas9n, spcas9v, xmacas9)",
    ),
    genome: str = typer.Option(
        None, "--genome", "-g",
        help="Path to reference genome FASTA for off-target search",
    ),
    max_mismatches: int = typer.Option(
        3, "--max-mismatches", "-m", help="Maximum mismatches for off-target search"
    ),
    max_results: int = typer.Option(
        50, "--max-results", "-n", help="Maximum number of guides to return"
    ),
    output: str = typer.Option(
        None, "--output", "-o", help="Output file path (JSON/CSV/HTML). Default: stdout"
    ),
    output_format: str = typer.Option(
        "json", "--format", "-f", help="Output format (json, csv, html)"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """Design gRNAs from a DNA sequence and rank by specificity."""
    from crispr_design.core.grna import GrnaType
    from crispr_design.core.design import design_guides, filter_guides
    from crispr_design.core.off_target import find_off_targets_brute_force
    from crispr_design.databases.genome import create_test_genome
    from crispr_design.reports.json_report import generate_json_report
    from crispr_design.reports.csv_report import generate_csv_report
    from crispr_design.reports.html_report import generate_html_report

    _setup_logging(verbose)

    # Determine gRNA type
    type_map = {
        "spcas9": GrnaType.SPCAS9,
        "sacas9": GrnaType.SACAS9,
        "spcas9n": GrnaType.SPCAS9N,
        "hificas9": GrnaType.HIFICAS9,
        "spcas9v": GrnaType.SPCAS9V,
        "xmacas9": GrnaType.XMACAS9,
    }
    gRNA_type = type_map.get(gRNA_type.lower(), GrnaType.SPCAS9)

    # Load sequence (from string or FASTA file)
    seq_path = Path(sequence) if sequence else None
    if seq_path and seq_path.exists():
        from crispr_design.utils.seq_utils import parse_fasta
        fasta_seqs = parse_fasta(str(seq_path))
        # Use first sequence or concatenate
        seq_parts = []
        for name, seq in fasta_seqs.items():
            seq_parts.append(seq.upper())
            typer.echo(f"Loaded {name} ({len(seq)} bp)", err=True)
        target_sequence = "".join(seq_parts)
    else:
        target_sequence = sequence.upper().strip() if sequence else ""

    if not target_sequence:
        typer.echo("Error: No sequence provided. Pass a DNA sequence or FASTA file path.", err=True)
        raise typer.Exit(1)

    # Load genome if provided
    genome_dict = None
    if genome:
        typer.echo(f"Loading genome from {genome}...", err=True)
        from crispr_design.databases.genome import load_reference_genome
        genome_dict = load_reference_genome(genome)
        typer.echo(f"Loaded {len(genome_dict)} chromosomes", err=True)

    # Design guides
    typer.echo(f"Scanning {len(target_sequence)} bp sequence for {gRNA_type.value} guides...", err=True)
    results = design_guides(
        sequence=target_sequence,
        gRNA_type=gRNA_type,
        genome=genome_dict,
        max_mismatches=max_mismatches,
        max_results=max_results,
    )

    if not results:
        typer.echo("No valid guides found. Check your sequence and PAM pattern.", err=True)
        raise typer.Exit(1)

    typer.echo(f"Found {len(results)} valid guides", err=True)

    # Output results
    if output:
        if output_format == "json":
            generate_json_report(results, output)
        elif output_format == "csv":
            generate_csv_report(results, output)
        elif output_format == "html":
            generate_html_report(results, output)
        else:
            typer.echo(f"Unknown format: {output_format}", err=True)
            raise typer.Exit(1)
        typer.echo(f"Report written to {output}", err=True)
    else:
        # Default: JSON to stdout
        output_data = [r.to_dict() for r in results]
        typer.echo(json.dumps(output_data, indent=2, default=str))


# =============================================================================
# scan subcommand
# =============================================================================
@cli.command()
def scan(
    guide_seq: str = typer.Argument(..., help="gRNA sequence (20bp)"),
    genome: str = typer.Option(
        None, "--genome", "-g",
        help="Path to reference genome FASTA"
    ),
    gRNA_type: str = typer.Option(
        "spcas9", "--type", "-t",
        help="Cas enzyme type"
    ),
    pam: str = typer.Option(
        "NGG", "--pam", "-p",
        help="PAM sequence"
    ),
    max_mismatches: int = typer.Option(
        3, "--max-mismatches", "-m",
        help="Maximum mismatches allowed"
    ),
    output: str = typer.Option(
        None, "--output", "-o",
        help="Output file path"
    ),
    output_format: str = typer.Option(
        "json", "--format", "-f",
        help="Output format (json, csv, html)"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """Scan a genome for off-target sites of a given gRNA."""
    from crispr_design.core.grna import Grna, GrnaType
    from crispr_design.core.off_target import find_off_targets_brute_force
    from crispr_design.reports.json_report import generate_json_report
    from crispr_design.reports.csv_report import generate_off_target_csv

    _setup_logging(verbose)

    type_map = {
        "spcas9": GrnaType.SPCAS9,
        "sacas9": GrnaType.SACAS9,
        "spcas9n": GrnaType.SPCAS9N,
        "hificas9": GrnaType.HIFICAS9,
        "spcas9v": GrnaType.SPCAS9V,
        "xmacas9": GrnaType.XMACAS9,
    }
    gRNA_type = type_map.get(gRNA_type.lower(), GrnaType.SPCAS9)

    # Create guide
    guide = Grna(sequence=guide_seq.upper(), pam_sequence=pam.upper(), gRNA_type=gRNA_type)

    # Load genome
    if genome:
        from crispr_design.databases.genome import load_reference_genome
        genome_dict = load_reference_genome(genome)
    else:
        typer.echo("No genome provided. Use --genome with a FASTA file.", err=True)
        raise typer.Exit(1)

    # Scan
    hits = find_off_targets_brute_force(
        guide, genome_dict, max_mismatches=max_mismatches,
        include_on_target=True, max_results=500
    )

    typer.echo(f"Found {len(hits)} off-target sites (≤{max_mismatches} mismatches)", err=True)

    if output:
        if output_format == "json":
            report = {
                "guide": guide.to_dict(),
                "hits": [h.to_dict() for h in hits],
                "total_hits": len(hits),
            }
            with open(output, "w") as f:
                json.dump(report, f, indent=2, default=str)
        elif output_format == "csv":
            generate_off_target_csv(hits, output, guide.name)
        elif output_format == "html":
            from crispr_design.core.scoring import compute_all_scores
            scores = compute_all_scores(guide.sequence)
            from crispr_design.core.design import DesignResult
            results = [DesignResult(guide=guide, scores=scores, off_targets=hits)]
            from crispr_design.reports.html_report import generate_html_report
            generate_html_report(results, output)
        typer.echo(f"Results written to {output}", err=True)
    else:
        for h in hits:
            typer.echo(
                f"{h.chromosome}:{h.position}({h.strand}) "
                f"{h.target_sequence} "
                f"[{h.mismatches} mismatches at {h.mismatch_positions}]"
            )


# =============================================================================
# batch subcommand
# =============================================================================
@cli.command()
def batch(
    input_file: str = typer.Argument(..., help="File with gRNA sequences (one per line)"),
    genome: str = typer.Option(
        None, "--genome", "-g",
        help="Path to reference genome FASTA"
    ),
    gRNA_type: str = typer.Option(
        "spcas9", "--type", "-t",
        help="Cas enzyme type"
    ),
    max_mismatches: int = typer.Option(
        3, "--max-mismatches", "-m",
        help="Maximum mismatches"
    ),
    output: str = typer.Option(
        "batch_results.json", "--output", "-o",
        help="Output file path"
    ),
    output_format: str = typer.Option(
        "json", "--format", "-f",
        help="Output format (json, csv, html)"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """Batch process multiple gRNAs from a file."""
    from crispr_design.core.grna import Grna, GrnaType
    from crispr_design.core.off_target import find_off_targets_brute_force
    from crispr_design.core.scoring import compute_all_scores
    from crispr_design.core.design import DesignResult
    from crispr_design.reports.html_report import generate_html_report

    _setup_logging(verbose)

    type_map = {
        "spcas9": GrnaType.SPCAS9,
        "sacas9": GrnaType.SACAS9,
        "spcas9n": GrnaType.SPCAS9N,
        "hificas9": GrnaType.HIFICAS9,
        "spcas9v": GrnaType.SPCAS9V,
        "xmacas9": GrnaType.XMACAS9,
    }
    gRNA_type = type_map.get(gRNA_type.lower(), GrnaType.SPCAS9)

    # Read guide sequences
    sequences = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                seq = line.split()[0]  # Handle TSV/CSV
                if len(seq) >= 20:
                    sequences.append(seq[:20])

    if not sequences:
        typer.echo("No valid sequences found in input file.", err=True)
        raise typer.Exit(1)

    typer.echo(f"Processing {len(sequences)} guides...", err=True)

    # Load genome
    genome_dict = None
    if genome:
        from crispr_design.databases.genome import load_reference_genome
        genome_dict = load_reference_genome(genome)

    results = []
    for i, seq in enumerate(sequences):
        typer.echo(f"  [{i + 1}/{len(sequences)}] {seq}", err=True)
        try:
            guide = Grna(sequence=seq, gRNA_type=gRNA_type, name=f"guide_{i + 1}")
            scores = compute_all_scores(guide.sequence)

            off_targets = []
            if genome_dict:
                off_targets = find_off_targets_brute_force(
                    guide, genome_dict, max_mismatches=max_mismatches,
                    include_on_target=False, max_results=50
                )

            results.append(DesignResult(
                guide=guide,
                scores=scores,
                off_targets=off_targets,
                specificity_score=max(0, 100 - len(off_targets) * 10),
                rank=i + 1,
            ))
        except ValueError as e:
            typer.echo(f"  Skipped {seq}: {e}", err=True)

    # Output
    if output_format == "html":
        generate_html_report(results, output, title="CRISPR Batch Analysis")
    elif output_format == "csv":
        from crispr_design.reports.csv_report import generate_csv_report
        generate_csv_report(results, output)
    else:
        report = {
            "total_guides": len(results),
            "guides": [r.to_dict() for r in results],
        }
        with open(output, "w") as f:
            json.dump(report, f, indent=2, default=str)

    typer.echo(f"Batch complete. Results: {output}", err=True)


# =============================================================================
# report subcommand
# =============================================================================
@cli.command()
def report(
    input_file: str = typer.Argument(..., help="JSON results file from design/batch"),
    output: str = typer.Option(
        "report.html", "--output", "-o",
        help="Output HTML report path"
    ),
    title: str = typer.Option(
        "CRISPR Design Studio Report", "--title",
        help="Report title"
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """Generate HTML report from JSON results file."""
    from crispr_design.core.grna import Grna, GrnaType
    from crispr_design.core.scoring import ScoreResult
    from crispr_design.core.design import DesignResult
    from crispr_design.reports.html_report import generate_html_report

    _setup_logging(verbose)

    with open(input_file) as f:
        data = json.load(f)

    # Parse results — handle both raw list and batch format
    if "guides" in data:
        guide_dicts = data["guides"]
    elif "results" in data:
        guide_dicts = data["results"]
    else:
        guide_dicts = data if isinstance(data, list) else [data]

    results = []
    for gd in guide_dicts:
        g = gd.get("guide", gd)
        guide = Grna(
            sequence=g.get("sequence", ""),
            pam_sequence=g.get("pam_sequence", ""),
            gRNA_type=GrnaType(g.get("gRNA_type", "spcas9")),
            name=g.get("name", ""),
        )
        scores_data = gd.get("scores", {})
        scores = ScoreResult(
            guide_name=g.get("name", ""),
            mit_cf=scores_data.get("mit_cf", 0),
            cfd=scores_data.get("cfd", 0),
            doench_2014=scores_data.get("doench_2014", 0),
            doench_2016=scores_data.get("doench_2016", 0),
        )
        r = DesignResult(
            guide=guide,
            scores=scores,
            off_targets=[],
            specificity_score=gd.get("specificity_score", 0),
            rank=gd.get("rank", len(results) + 1),
        )
        # Restore off-targets if present
        for ot in gd.get("off_targets", []):
            from crispr_design.core.off_target import OffTargetHit
            hit = OffTargetHit(
                chromosome=ot.get("chromosome", ""),
                position=ot.get("position", 0),
                strand=ot.get("strand", "+"),
                target_sequence=ot.get("target_sequence", ""),
                mismatches=ot.get("mismatches", 0),
                mismatch_positions=ot.get("mismatch_positions", []),
                mismatch_details=ot.get("mismatch_details", []),
            )
            r.off_targets.append(hit)
        results.append(r)

    generate_html_report(results, output, title=title)
    typer.echo(f"Report written to {output}", err=True)


# =============================================================================
# server subcommand
# =============================================================================
@cli.command()
def serve(
    host: str = typer.Option("127.0.0.1", "--host", help="Host to bind"),
    port: int = typer.Option(8000, "--port", "-p", help="Port to bind"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),
):
    """Start the FastAPI server."""
    import uvicorn
    _setup_logging(verbose)

    typer.echo(f"Starting CRISPR Design Studio API on {host}:{port}", err=True)
    uvicorn.run(
        "crispr_design.api.server:app",
        host=host,
        port=port,
        reload=False,
        log_level="debug" if verbose else "info",
    )


def main():
    """Entry point — calls the typer app to avoid name collision."""
    app = typer.main.get_command(cli)
    app()


if __name__ == "__main__":
    main()
