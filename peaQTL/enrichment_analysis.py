"""
enrichment_analysis.py

Runs pathway over-representation analysis for cell-type-specific ABC target
genes using Enrichr gene sets through gseapy.
"""

import logging
from pathlib import Path
from typing import Iterable

import gseapy as gp
import pandas as pd

from constants import (
	ABC_DICT_PATH,
	PEAK2GENE_OUTPUT_DIR,
	PEAQTL_RESULTS_DIR,
	SIGNIFICANT_PEAK_CELL_TYPES,
)

logging.basicConfig(
	level=logging.INFO,
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

GENE_SETS = ["KEGG_2021_Human", "GO_Biological_Process_2023"]
TARGET_GENE_COLUMN = "TargetGene"
ADJUSTED_P_VALUE_COLUMN = "Adjusted P-value"
TERM_COLUMN = "Term"

ENRICHMENT_OUTPUT_DIR = PEAQTL_RESULTS_DIR / "pathway_enrichment"
ENRICHMENT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def clean_gene_list(genes: Iterable[object]) -> list[str]:
	"""Return unique non-empty gene symbols while preserving first-seen order."""
	cleaned: list[str] = []
	seen: set[str] = set()

	for gene in genes:
		if pd.isna(gene):
			continue

		gene_symbol = str(gene).strip()
		if not gene_symbol or gene_symbol in seen:
			continue

		seen.add(gene_symbol)
		cleaned.append(gene_symbol)

	return cleaned


def load_target_genes(path: Path) -> list[str]:
	"""Load a CSV file and extract clean unique TargetGene values."""
	if not path.exists():
		raise FileNotFoundError(f"Input file not found: {path}")

	df = pd.read_csv(path, usecols=[TARGET_GENE_COLUMN])
	return clean_gene_list(df[TARGET_GENE_COLUMN])


def load_background_genes(path: Path) -> list[str]:
	logger.info("Loading custom liver background from: %s", path)
	genes = load_target_genes(path)
	if not genes:
		raise ValueError(f"No background genes found in {path}")
	logger.info("Custom liver background contains %d unique genes.", len(genes))
	return genes


def run_enrichment(gene_list: list[str], background: list[str]):
	return gp.enrichr(
		gene_list=gene_list,
		gene_sets=GENE_SETS,
		background=background,
		outdir=None,
	)


def log_top_pathways(results: pd.DataFrame, cell_type: str) -> None:
	top_pathways = results[[TERM_COLUMN, ADJUSTED_P_VALUE_COLUMN]].head(3)
	if top_pathways.empty:
		logger.info("  No significant pathways found for %s.", cell_type)
		return

	logger.info("  Top significant pathways for %s:", cell_type)
	for _, row in top_pathways.iterrows():
		logger.info(
			"    %s | adjusted p-value: %.3g",
			row[TERM_COLUMN],
			row[ADJUSTED_P_VALUE_COLUMN],
		)


def run() -> None:
	background = load_background_genes(ABC_DICT_PATH)

	for cell_type in SIGNIFICANT_PEAK_CELL_TYPES:
		logger.info("Processing cell type: %s", cell_type)
		input_path = PEAK2GENE_OUTPUT_DIR / f"{cell_type}_peak2gene.csv"
		output_path = ENRICHMENT_OUTPUT_DIR / f"{cell_type}_significant_pathways.csv"

		try:
			genes = load_target_genes(input_path)
		except (FileNotFoundError, ValueError) as exc:
			logger.warning("  Skipping %s: %s", cell_type, exc)
			continue

		logger.info("  Unique input genes: %d", len(genes))
		if not genes:
			logger.warning("  Skipping %s because it has no target genes.", cell_type)
			continue

		try:
			enr = run_enrichment(genes, background)
		except Exception as exc:
			logger.exception("  Enrichment failed for %s: %s", cell_type, exc)
			continue

		results = enr.results
		if results is None or results.empty:
			logger.warning("  Enrichr returned no results for %s.", cell_type)
			continue

		significant = (
			results[results[ADJUSTED_P_VALUE_COLUMN] < 0.05]
			.sort_values(ADJUSTED_P_VALUE_COLUMN)
			.reset_index(drop=True)
		)
		significant.to_csv(output_path, index=False)
		logger.info("  Saved %d significant pathways to: %s", len(significant), output_path)
		log_top_pathways(significant, cell_type)


if __name__ == "__main__":
	run()
