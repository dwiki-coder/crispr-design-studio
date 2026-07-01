"""HTML report generation with Jinja2 templates."""

from __future__ import annotations

import os
from datetime import datetime
from typing import Any, Optional

from ..core.design import DesignResult
from ..core.off_target import OffTargetHit


def generate_html_report(
    results: list[DesignResult],
    output_path: str,
    title: str = "CRISPR Design Studio Report",
    metadata: Optional[dict[str, Any]] = None,
    template: Optional[str] = None,
) -> str:
    """Generate a publication-ready HTML report.

    Args:
        results: List of DesignResult objects.
        output_path: Path to write HTML file.
        title: Report title.
        metadata: Optional metadata to display.
        template: Optional path to custom Jinja2 template.

    Returns:
        Path to the generated report.
    """
    from jinja2 import Environment, BaseLoader, PackageLoader, select_autoescape

    # Try to load custom template, fall back to built-in
    if template and os.path.exists(template):
        env = Environment(
            loader=BaseLoader(load_string=_read_template(template)),
            autoescape=select_autoescape(["html"]),
        )
        tpl = env.get_template("")
    else:
        env = Environment(
            loader=BaseLoader(),
            autoescape=select_autoescape(["html"]),
        )
        tpl = env.from_string(HTML_TEMPLATE)

    # Prepare data
    guide_data = []
    for r in results:
        guide_data.append({
            "rank": r.rank,
            "name": r.guide.name,
            "sequence": r.guide.sequence,
            "pam": r.guide.pam_sequence,
            "gtype": r.guide.gRNA_type.value,
            "gc_content": round(r.guide.gc_content * 100, 1),
            "rc": r.guide.reverse_complement(),
            "mit_cf": round(r.scores.mit_cf, 4),
            "cfd": round(r.scores.cfd, 4),
            "doench_2014": round(r.scores.doench_2014, 4),
            "doench_2016": round(r.scores.doench_2016, 4),
            "specificity": round(r.specificity_score, 2),
            "off_target_count": len(r.off_targets),
            "off_targets": _format_off_targets(r.off_targets[:10]),
        })

    rendered = tpl.render(
        title=title,
        generated=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        total_guides=len(results),
        guides=guide_data,
        metadata=metadata or {},
    )

    with open(output_path, "w") as f:
        f.write(rendered)

    return output_path


def _format_off_targets(hits: list[OffTargetHit]) -> list[dict]:
    """Format off-target hits for HTML display."""
    formatted = []
    for hit in hits:
        formatted.append({
            "chromosome": hit.chromosome,
            "position": hit.position,
            "strand": hit.strand,
            "sequence": hit.target_sequence,
            "mismatches": hit.mismatches,
            "positions": ", ".join(str(p) for p in hit.mismatch_positions),
            "details": hit.mismatch_details,
        })
    return formatted


def _read_template(path: str) -> str:
    """Read template file."""
    with open(path) as f:
        return f.read()


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ title }}</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            color: #333;
            line-height: 1.6;
        }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .summary {
            background: #f8f9fa;
            padding: 15px 20px;
            border-radius: 8px;
            margin: 20px 0;
            border-left: 4px solid #3498db;
        }
        .summary span { margin-right: 30px; font-weight: 500; }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 14px;
        }
        th {
            background: #3498db;
            color: white;
            padding: 10px 12px;
            text-align: left;
            font-weight: 600;
        }
        td {
            padding: 8px 12px;
            border-bottom: 1px solid #e9ecef;
        }
        tr:nth-child(even) { background: #f8f9fa; }
        tr:hover { background: #e8f4f8; }
        .sequence { font-family: 'Courier New', monospace; font-weight: bold; color: #e74c3c; }
        .good { color: #27ae60; font-weight: bold; }
        .bad { color: #e74c3c; font-weight: bold; }
        .metadata { background: #fff3cd; padding: 10px 15px; border-radius: 5px; margin: 10px 0; }
        .footer { margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #666; font-size: 12px; }
        .off-target-table { width: 100%; font-size: 12px; margin-top: 5px; }
        .off-target-table th { background: #95a5a6; padding: 4px 8px; font-size: 11px; }
        .off-target-table td { padding: 3px 8px; }
    </style>
</head>
<body>
    <h1>{{ title }}</h1>

    <div class="summary">
        <span>Total Guides: {{ total_guides }}</span>
        <span>Generated: {{ generated }}</span>
    </div>

    {% if metadata %}
    <div class="metadata">
        <strong>Metadata:</strong>
        {% for key, value in metadata.items() %}
        <div>{{ key }}: {{ value }}</div>
        {% endfor %}
    </div>
    {% endif %}

    <h2>Guide Rankings</h2>
    <table>
        <thead>
            <tr>
                <th>Rank</th>
                <th>Name</th>
                <th>Sequence (5'→3')</th>
                <th>PAM</th>
                <th>Type</th>
                <th>GC%</th>
                <th>MIT-CF</th>
                <th>CFD</th>
                <th>Doe 2014</th>
                <th>Doe 2016</th>
                <th>Specificity</th>
                <th>OT Sites</th>
            </tr>
        </thead>
        <tbody>
            {% for g in guides %}
            <tr>
                <td>{{ g.rank }}</td>
                <td>{{ g.name }}</td>
                <td class="sequence">{{ g.sequence }}</td>
                <td>{{ g.pam }}</td>
                <td>{{ g.gtype }}</td>
                <td>{{ g.gc_content }}%</td>
                <td>{{ g.mit_cf }}</td>
                <td>{{ g.cfd }}</td>
                <td>{{ g.doench_2014 }}</td>
                <td>{{ g.doench_2016 }}</td>
                <td class="{{ 'good' if g.specificity > 50 else 'bad' }}">{{ g.specificity }}</td>
                <td>{{ g.off_target_count }}</td>
            </tr>
            {% if g.off_targets %}
            <tr>
                <td colspan="12">
                    <table class="off-target-table">
                        <thead>
                            <tr>
                                <th>Chrom</th>
                                <th>Position</th>
                                <th>Strand</th>
                                <th>Target Sequence</th>
                                <th>Mismatches</th>
                                <th>Positions</th>
                            </tr>
                        </thead>
                        <tbody>
                            {% for ot in g.off_targets %}
                            <tr>
                                <td>{{ ot.chromosome }}</td>
                                <td>{{ ot.position }}</td>
                                <td>{{ ot.strand }}</td>
                                <td class="sequence" style="font-size:11px">{{ ot.sequence }}</td>
                                <td>{{ ot.mismatches }}</td>
                                <td>{{ ot.positions }}</td>
                            </tr>
                            {% endfor %}
                        </tbody>
                    </table>
                </td>
            </tr>
            {% endif %}
            {% endfor %}
        </tbody>
    </table>

    <div class="footer">
        Generated by CRISPR Design Studio v0.1.0 |
        {{ generated }} |
        {{ total_guides }} guides analyzed
    </div>
</body>
</html>"""
