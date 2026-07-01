.PHONY: install dev-install test lint format clean docker-build help

SHELL := /bin/bash
PROJECT_NAME := crispr-design
PYTHON := python3
PYTEST := pytest

install: ## Install dependencies
	$(PYTHON) -m pip install .

dev-install: ## Install in development mode with dev dependencies
	$(PYTHON) -m pip install -e ".[dev,api]"

test: ## Run tests
	$(PYTEST) tests/ -v --tb=short

test-cov: ## Run tests with coverage
	$(PYTEST) tests/ -v --tb=short --cov=crispr_design --cov-report=html

lint: ## Run linters
	ruff check crispr_design/

format: ## Format code
	black crispr_design/ tests/
	ruff check crispr_design/ --fix

clean: ## Clean build artifacts
	rm -rf build/ dist/ *.egg-info .pytest_cache .mypy_cache __pycache__
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

docker-build: ## Build Docker image
	docker build -t crispr-design:latest .

demo: ## Run a demo design
	crispr-design design "AAGCAGTGGTATCAAGTCAGAGG" --max-results 5

help: ## Show help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
