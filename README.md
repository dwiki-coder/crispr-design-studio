# CRISPR Design Studio

[![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Tests](https://github.com/dwiki-coder/crispr-design-studio/actions/workflows/test.yml/badge.svg)](https://github.com/dwiki-coder/crispr-design-studio/actions/workflows/test.yml)
[![Coverage](https://img.shields.io/badge/Coverage-92%25-brightgreen)](https://github.com/dwiki-coder/crispr-design-studio)

> CLI + API tool for gRNA design and off-target prediction.
> A scriptable, CLI-first alternative to web-only CRISPR design tools.

CRISPR Design Studio predicts off-target effects for CRISPR guide RNAs with
multiple scoring algorithms, batch processing, and publication-ready reports.
Unlike web-only tools (CRISPOR, CHOPCHOP), this is a **CLI-first, scriptable**
solution that integrates directly into bioinformatics pipelines.

---

## Why This Tool

CRISPOR and CHOPCHOP are the dominant CRISPR design tools, but they are
web-only — you can't script them, batch them, or integrate them into automated
pipelines. CRISPR Design Studio fills this gap with a CLI-first, API-accessible
solution supporting **6 Cas variants** (including HiFi, nickase, HySt, and xMACs)
and **4 independent scoring algorithms**. Researchers can design hundreds of
guides programmatically, scan whole genomes for off-targets, and generate
publication-ready reports — all from the command line.

---

## Metrics

| Metric | Value |
|--------|-------|
| Automated tests | **197** across 9 test files |
| Test-to-source ratio | **62%** (1,653 test LOC / 2,675 source LOC) |
| Cas variants supported | **6** — SpCas9, SaCas9, HiFi Cas9, SpCas9 nickase, SpCas9-HySt, xMACs |
| Scoring algorithms | **4** — MIT-CF, CFD (Hsu 2013), Doench 2014, Doench 2016 |
| Off-target search | Seed-and-extend with configurable mismatch tolerance |
| Batch mode | Hundreds of guides per run |
| Output formats | **3** — JSON, CSV, HTML |
| REST API endpoints | **7** |
| Source code | **4,328 LOC** across 28 Python files |

---

## Who Should Use This

- **Genome editing researchers** designing CRISPR experiments and screening gRNAs
- **Gene therapy companies** building automated, validated editing pipelines
- **Academic core facilities** offering CRISPR design services to multiple labs
- **Bioinformaticians** integrating CRISPR design into broader analysis workflows
- **CRISPR-based drug developers** needing programmatic, reproducible gRNA design

---

## Architecture

```
┌──────────────────────────────────────────────────────────────────┐
│                   User Interface Layer                           │
│  ┌──────────┐  ┌──────────┐  ┌────────────────────────────┐     │
│  │   CLI    │  │  REST    │  │  Python Library (import)    │     │
│  │  (Click) │  │  API     │  │                              │     │
│  │          │  │(FastAPI) │  │  from crispr_design import…  │     │
│  └────┬─────┘  └────┬─────┘  └──────────┬─────────────────┘     │
│       │             │                    │                       │
│       └──────────────┼────────────────────┘                      │
│                      ▼                                            │
│  ┌────────────────────────────────────────────────────────┐      │
│  │             Core Engine (gRNA Design)                   │      │
│  │  ┌──────────┐  ┌──────────┐  ┌──────────┐              │      │
│  │  │ gRNA.py │  │scoring.py│  │off_target│              │      │
│  │  │(PAM     │  │(MIT-CF,  │  │  search  │              │      │
│  │  │detect)  │  │ CFD,     │  │(seed&    │              │      │
│  │  │         │  │Doench)   │  │ extend)  │              │      │
│  │  └──────────┘  └──────────┘  └──────────┘              │      │
│  └─────────────────────────────────────────────────────────┘      │
│                      │                                            │
│  ┌───────────────────┴────────────────────────────────────────┐  │
│  │              Data & Infrastructure Layer                    │  │
│  │  ┌────────────┐  ┌──────────┐  ┌─────┐  ┌─────────────┐  │  │
│  │  │ Biopython  │  │ Polars  │  │Docker│  │ Jinja2       │  │  │
│  │  │ FASTA idx  │  │ DataFrames│  │     │  │ HTML reports │  │  │
│  │  └────────────┘  └──────────┘  └─────┘  └─────────────┘  │  │
│  └───────────────────────────────────────────────────────────┘  │
└──────────────────────────────────────────────────────────────────┘
```

### Pipeline Flow

```
Target Sequence ──► PAM Detection ──► gRNA Extraction
                        │                   │
                        ▼                   ▼
              Cas Variant            Candidate gRNAs
              Matching (6 types)          │
                                          ▼
                                    Scoring (4 algorithms)
                                    ┌────┴───────────────┐
                              MIT-CF  CFD  Doench'14  Doench'16
                                          │
                                          ▼
                                    ┌─────┴─────┐
                                    │ Off-target │
                                    │  Search    │
                                    │(seed&ext)  │
                                    └─────┬──────┘
                                          ▼
                              Ranked Guide List + Report
                              (JSON / CSV / HTML)
```

---

## Features

- **Multi-PAM support**: SpCas9 (NGG), SaCas9 (NNGRRT), HiFi Cas9, SpCas9 nickase, SpCas9-HySt, xMACs
- **Scoring algorithms**: MIT-CF, CFD (Hsu 2013), Doench 2014, Doench 2016
- **Off-target search**: Seed-and-extend algorithm with configurable mismatch tolerance
- **Batch mode**: Process hundreds of guides at once
- **Report generation**: JSON, CSV, and publication-ready HTML
- **FastAPI server**: REST API for programmatic access

---

## Installation

```bash
# Clone and install
git clone https://github.com/dwiki-coder/crispr-design-studio.git
cd crispr-design-studio
pip install -e ".[all]"
```

---

## Quick Start

### Design gRNAs from a sequence

```bash
# Scan a sequence for valid gRNAs
crispr-design design "AAGCAGTGGTATCAAGTCAGAGG" \
    --type spcas9 \
    --max-results 10

# With genome off-target search
crispr-design design "AAGCAGTGGTATCAAGTCAGAGG" \
    --genome hg38.fa \
    --max-mismatches 3 \
    -o results.json -f json
```

### Scan for off-target sites

```bash
crispr-design scan "AAGCAGTGGTATCAAGTCAG" \
    --genome hg38.fa \
    --pam NGG \
    --max-mismatches 3 \
    -o off_targets.csv -f csv
```

### Batch process guides

```bash
# Create input file with one guide per line
echo "AAGCAGTGGTATCAAGTCAG
TTCGATCGATCGATCGATCG
GGCCAATTCCGGCCAATTCC" > guides.txt

crispr-design batch guides.txt \
    --genome hg38.fa \
    --max-mismatches 3 \
    -o batch_results.json
```

### Generate HTML report

```bash
crispr-design report results.json -o report.html
```

### Start API server

```bash
crispr-design serve --host 0.0.0.0 --port 8000
```

Then use the API:

```bash
# Design guides
curl -X POST http://localhost:8000/design \
  -H "Content-Type: application/json" \
  -d '{"sequence": "AAGCAGTGGTATCAAGTCAGAGG", "gRNA_type": "spcas9"}'

# Scan for off-targets
curl -X POST http://localhost:8000/scan \
  -H "Content-Type: application/json" \
  -d '{"guide_sequence": "AAGCAGTGGTATCAAGTCAG", "pam": "GG"}'

# Batch process
curl -X POST http://localhost:8000/batch \
  -H "Content-Type: application/json" \
  -d '{"sequences": ["AAGCAGTGGTATCAAGTCAG", "TTCGATCGATCGATCGATCG"]}'
```

---

## Supported Cas Variants

| Variant | PAM Pattern | Guide Length |
|---------|-------------|-------------|
| SpCas9 | NGG | 20bp |
| SaCas9 | NNGRRT | 20bp |
| HiFi Cas9 | NGG | 20bp |
| SpCas9 nickase | NGG/NAG | 20bp |
| SpCas9-HySt | NRY | 20bp |
| xMACs Cas9 | NGG | 20bp |

---

## Scoring Algorithms

| Algorithm | Reference | Range | Interpretation |
|-----------|-----------|-------|----------------|
| MIT-CF | Cong et al. 2013 | Lower = better | Off-target propensity |
| CFD | Hsu et al. 2013 | 0-1 (1=perfect match) | Mismatch impact |
| Doench 2014 | Doench et al. 2014 | -20 to +40 | On-target activity |
| Doench 2016 | Doench et al. 2016 | -20 to +40 | Improved activity prediction |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Health check |
| `/genomes` | GET | List available genomes |
| `/genomes/upload` | POST | Upload custom genome FASTA |
| `/design` | POST | Design gRNAs from sequence |
| `/scan` | POST | Scan for off-target sites |
| `/batch` | POST | Batch process multiple guides |
| `/score` | POST | Score a single gRNA |

---

## Docker

```bash
# Build
docker build -t crispr-design .

# Run CLI
docker run crispr-design design "AAGCAGTGGTATCAAGTCAGAGG"

# Run API server
docker run -p 8000:8000 crispr-design serve --host 0.0.0.0 --port 8000
```

---

## Development

```bash
# Install dev dependencies
make dev-install

# Run tests
make test

# Run tests with coverage
make test-cov

# Format code
make format

# Lint
make lint
```

---

## Citation

```bibtex
@software{crispr_design_studio,
  author       = {David},
  title        = {{CRISPR Design Studio: gRNA Design and Off-Target Prediction}},
  year         = {2026},
  url          = {https://github.com/dwiki-coder/crispr-design-studio},
  license      = {MIT}
}
```

---

## Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## License

MIT
