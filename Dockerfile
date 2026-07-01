FROM python:3.11-slim

LABEL maintainer="CRISPR Design Studio Contributors"
LABEL description="CRISPR gRNA design and off-target prediction tool"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy and install dependencies first (layer caching)
COPY pyproject.toml ./
RUN pip install --no-cache-dir .

# Copy application code
COPY . .

# Default command
ENTRYPOINT ["crispr-design"]
CMD ["--help"]

# API server mode
# docker run -p 8000:8000 crispr-design serve --host 0.0.0.0 --port 8000
