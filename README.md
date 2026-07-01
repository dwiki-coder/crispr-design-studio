# CRISPR Design Studio

> CLI + API tool for gRNA design and off-target prediction

CRISPR Design Studio predicts off-target effects for CRISPR guide RNAs with multiple scoring algorithms, batch processing, and publication-ready reports. Unlike web-only tools (CRISPOR, CHOPCHOP), this is a **CLI-first, scriptable** solution.

## Features

- **Multi-PAM support**: SpCas9 (NGG), SaCas9 (NNGRRT), HiFi Cas9, SpCas9 nickase, SpCas9-HySt, xMACs
- **Scoring algorithms**: MIT-CF, CFD (Hsu 2013), Doench 2014, Doench 2016
- **Off-target search**: Seed-and-extend algorithm with configurable mismatch tolerance
- **Batch mode**: Process hundreds of guides at once
- **Report generation**: JSON, CSV, and publication-ready HTML
- **FastAPI server**: REST API for programmatic access

## Installation

```bash
# Standard installation
pip install -e .

# With API server support
pip install -e ".[dev,api]"
```

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

## Supported Cas Variants

| Variant | PAM Pattern | Guide Length |
|---------|-------------|-------------|
| SpCas9 | NGG | 20bp |
| SaCas9 | NNGRRT | 20bp |
| HiFi Cas9 | NGG | 20bp |
| SpCas9 nickase | NGG/NAG | 20bp |
| SpCas9-HySt | NRY | 20bp |
| xMACs Cas9 | NGG | 20bp |

## Scoring Algorithms

| Algorithm | Reference | Range | Interpretation |
|-----------|-----------|-------|----------------|
| MIT-CF | Cong et al. 2013 | Lower = better | Off-target propensity |
| CFD | Hsu et al. 2013 | 0-1 (1=perfect match) | Mismatch impact |
| Doench 2014 | Doench et al. 2014 | -20 to +40 | On-target activity |
| Doench 2016 | Doench et al. 2016 | -20 to +40 | Improved activity prediction |

## Docker

```bash
# Build
docker build -t crispr-design .

# Run CLI
docker run crispr-design design "AAGCAGTGGTATCAAGTCAGAGG"

# Run API server
docker run -p 8000:8000 crispr-design serve --host 0.0.0.0 --port 8000
```

## Project Structure

```
crispr-design/
├── crispr_design/
│   ├── core/           # Core analysis modules
│   │   ├── grna.py    # gRNA class & PAM detection
│   │   ├── off_target.py # Off-target search engine
│   │   ├── scoring.py # Scoring algorithms
│   │   └── design.py  # Guide optimization
│   ├── databases/      # Genome handling
│   │   └── genome.py  # FASTA loading & indexing
│   ├── reports/        # Report generation
│   │   ├── json_report.py
│   │   ├── csv_report.py
│   │   └── html_report.py
│   ├── utils/          # Utilities
│   │   ├── seq_utils.py
│   │   └── vcf_parser.py
│   ├── api/            # FastAPI server
│   │   └── server.py
│   └── cli.py          # CLI entry point
├── tests/              # Test suite (50+ tests)
├── pyproject.toml
├── Dockerfile
├── Makefile
└── README.md
```

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

## License

MIT
