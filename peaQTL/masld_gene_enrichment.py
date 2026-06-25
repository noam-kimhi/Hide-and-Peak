"""Hypergeometric test for MASLD gene enrichment among discovered target genes.

For each cell type (and combined across all), this script tests whether the
target genes discovered via the ABC peak-to-gene mapping are enriched for
known MASLD/NASH-associated genes, compared to the liver ABC dictionary
background.

Usage
-----
Default mode (reads from PEAK2GENE_OUTPUT_DIR):

    python -m peaQTL.masld_gene_enrichment

Soft mode (reads from SOFT_PEAK2GENE_OUTPUT_DIR):

    python -m peaQTL.masld_gene_enrichment --soft

The test:
    H0: Our discovered genes are a random draw from the ABC liver background.
    H1: Our discovered genes are enriched for known MASLD-associated genes.

    P-value = P(X >= k) where X ~ Hypergeometric(N, K, n)
        N = total unique genes in ABC liver background
        K = known MASLD genes present in background
        n = number of discovered genes
        k = number of discovered genes that are MASLD-associated
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Final

import pandas as pd
import requests
from scipy.stats import hypergeom, fisher_exact

from constants import (
    ABC_DICT_PATH,
    PEAK2GENE_OUTPUT_DIR,
    PEAK2GENE_19_OUTPUT_DIR,
    PEAQTL_RESULTS_DIR,
    SIGNIFICANT_PEAK_CELL_TYPES,
    SOFT_PEAK2GENE_OUTPUT_DIR,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

OUTPUT_DIR: Final[Path] = PEAQTL_RESULTS_DIR / "masld_enrichment"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Open Targets Platform API configuration
# Disease: MASLD (MONDO_0013209)
# https://platform.opentargets.org/disease/MONDO_0013209/associations
# ---------------------------------------------------------------------------

OPEN_TARGETS_API_URL: Final[str] = "https://api.platform.opentargets.org/api/v4/graphql"
MASLD_DISEASE_ID: Final[str] = "MONDO_0013209"

DISEASE_ASSOCIATIONS_QUERY: Final[str] = """
query DiseaseAssociationsQuery(
  $id: String!
  $index: Int!
  $size: Int!
  $sortBy: String!
  $enableIndirect: Boolean!
  $datasources: [DatasourceSettingsInput!]
) {
  disease(efoId: $id) {
    id
    name
    associatedTargets(
      page: { index: $index, size: $size }
      orderByScore: $sortBy
      enableIndirect: $enableIndirect
      datasources: $datasources
    ) {
      count
      rows {
        target {
          approvedSymbol
        }
        score
      }
    }
  }
}
"""

DATASOURCES: Final[list[dict]] = [
    {"id": "clinical_precedence", "weight": 1, "propagate": True, "required": False},
    {"id": "gwas_credible_sets", "weight": 1, "propagate": True, "required": False},
    {"id": "gene_burden", "weight": 1, "propagate": True, "required": False},
    {"id": "eva", "weight": 1, "propagate": True, "required": False},
    {"id": "genomics_england", "weight": 1, "propagate": True, "required": False},
    {"id": "gene2phenotype", "weight": 1, "propagate": True, "required": False},
    {"id": "uniprot_literature", "weight": 1, "propagate": True, "required": False},
    {"id": "uniprot_variants", "weight": 1, "propagate": True, "required": False},
    {"id": "orphanet", "weight": 1, "propagate": True, "required": False},
    {"id": "clingen", "weight": 1, "propagate": True, "required": False},
    {"id": "cancer_gene_census", "weight": 1, "propagate": True, "required": False},
    {"id": "intogen", "weight": 1, "propagate": True, "required": False},
    {"id": "eva_somatic", "weight": 1, "propagate": True, "required": False},
    {"id": "cancer_biomarkers", "weight": 1, "propagate": True, "required": False},
    {"id": "crispr_screen", "weight": 1, "propagate": True, "required": False},
    {"id": "crispr", "weight": 1, "propagate": True, "required": False},
    {"id": "reactome", "weight": 1, "propagate": True, "required": False},
    {"id": "europepmc", "weight": 0.2, "propagate": True, "required": False},
    {"id": "expression_atlas", "weight": 0.2, "propagate": False, "required": False},
    {"id": "impc", "weight": 0.2, "propagate": True, "required": False},
    {"id": "ot_crispr_validation", "weight": 0.5, "propagate": True, "required": False},
    {"id": "ot_crispr", "weight": 0.5, "propagate": True, "required": False},
    {"id": "encore", "weight": 0.5, "propagate": True, "required": False},
]


def fetch_masld_genes_from_open_targets(
    min_score: float = 0.0,
) -> frozenset[str]:
    """Fetch MASLD-associated genes from the Open Targets Platform API.

    Queries all associated targets for MONDO_0013209 (MASLD) and returns
    their approved gene symbols.

    Parameters
    ----------
    min_score : float
        Minimum overall association score to include a gene (default 0.0 = all).

    Returns
    -------
    frozenset of gene symbols associated with MASLD.
    """
    all_genes: list[str] = []
    page_size = 500
    page_index = 0

    logger.info(
        "Fetching MASLD-associated genes from Open Targets Platform "
        "(disease: %s)...",
        MASLD_DISEASE_ID,
    )

    while True:
        variables = {
            "id": MASLD_DISEASE_ID,
            "index": page_index,
            "size": page_size,
            "sortBy": "score",
            "enableIndirect": True,
            "datasources": DATASOURCES,
        }

        response = requests.post(
            OPEN_TARGETS_API_URL,
            json={"query": DISEASE_ASSOCIATIONS_QUERY, "variables": variables},
            timeout=60,
        )
        response.raise_for_status()
        data = response.json()

        disease_data = data.get("data", {}).get("disease")
        if disease_data is None:
            raise ValueError(
                f"No data returned for disease ID: {MASLD_DISEASE_ID}"
            )

        assoc = disease_data["associatedTargets"]
        total_count = assoc["count"]
        rows = assoc["rows"]

        if not rows:
            break

        genes_before = len(all_genes)
        for row in rows:
            score = row.get("score", 0.0)
            if score >= min_score:
                symbol = row["target"]["approvedSymbol"]
                if symbol:
                    all_genes.append(symbol)
        genes_added = len(all_genes) - genes_before

        logger.info(
            "  Fetched page %d (%d genes so far / %d total, +%d this page)",
            page_index,
            len(all_genes),
            total_count,
            genes_added,
        )

        # Results are sorted by score descending, so if no genes passed
        # the min_score filter on this page, no future page will either.
        if genes_added == 0 and min_score > 0:
            logger.info(
                "  Stopping early: remaining genes score below %.2f",
                min_score,
            )
            break

        page_index += 1
        if page_index * page_size >= total_count:
            break

    genes = frozenset(all_genes)
    logger.info(
        "Retrieved %d MASLD-associated genes from Open Targets (score >= %.2f)",
        len(genes),
        min_score,
    )
    return genes


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Hypergeometric test for MASLD gene enrichment."
    )
    parser.add_argument(
        "--soft",
        action="store_true",
        help="Use soft-filtered peak2gene results.",
    )
    parser.add_argument(
        "--old", action="store_true", help="Use old peak2gene results (before hg38 fix)."
    )
    return parser.parse_args()


def load_background_genes(path: Path) -> set[str]:
    """Load all unique gene symbols from the ABC liver dictionary."""
    if not path.is_file():
        raise FileNotFoundError(f"ABC dictionary not found: {path}")

    df = pd.read_csv(path, usecols=["TargetGene"])
    genes = set(
        df["TargetGene"].dropna().astype(str).str.strip().unique()
    )
    genes.discard("")
    return genes


def load_peak2gene_genes(path: Path) -> set[str]:
    """Load unique target gene symbols from a peak2gene CSV."""
    if not path.is_file():
        return set()

    df = pd.read_csv(path, usecols=["TargetGene"])
    genes = set(
        df["TargetGene"].dropna().astype(str).str.strip().unique()
    )
    genes.discard("")
    return genes


def run_hypergeometric_test(
    discovered_genes: set[str],
    background_genes: set[str],
    masld_genes: set[str],
) -> dict[str, object]:
    """Run one-sided hypergeometric enrichment test.

    Parameters
    ----------
    discovered_genes : genes found by our pipeline
    background_genes : all genes in the ABC liver dictionary
    masld_genes : curated MASLD-associated genes

    Returns
    -------
    Dict with test statistics and gene-level details.
    """
    # Restrict everything to the background universe
    masld_in_background = masld_genes & background_genes
    discovered_in_background = discovered_genes & background_genes
    overlap = discovered_in_background & masld_in_background

    N = len(background_genes)           # population size
    K = len(masld_in_background)        # successes in population
    n = len(discovered_in_background)   # draws
    k = len(overlap)                    # observed successes

    # P(X >= k) using survival function: sf(k-1, N, K, n)
    if n == 0 or K == 0:
        pvalue = 1.0
    else:
        pvalue = float(hypergeom.sf(k - 1, N, K, n))

    # Also compute Fisher's exact test (2x2 contingency table)
    # as a complementary measure
    #
    #                 MASLD    Not-MASLD
    # Discovered:      k        n - k
    # Not discovered:  K - k    N - K - (n - k)
    a = k
    b = n - k
    c = K - k
    d = N - K - b

    if n > 0 and N > 0:
        odds_ratio, fisher_p = fisher_exact(
            [[a, b], [c, d]], alternative="greater"
        )
    else:
        odds_ratio, fisher_p = float("nan"), 1.0

    # Expected overlap under null
    expected_k = (K / N) * n if N > 0 else 0.0
    fold_enrichment = k / expected_k if expected_k > 0 else float("inf") if k > 0 else 0.0

    return {
        "N_background": N,
        "K_masld_in_background": K,
        "n_discovered": n,
        "k_overlap": k,
        "expected_overlap": round(expected_k, 4),
        "fold_enrichment": round(fold_enrichment, 4),
        "hypergeom_pvalue": pvalue,
        "fisher_pvalue": fisher_p,
        "fisher_odds_ratio": odds_ratio,
        "overlap_genes": sorted(overlap),
        "discovered_genes_list": sorted(discovered_in_background),
    }


def print_result(cell_type: str, result: dict[str, object]) -> None:
    """Print a formatted summary for one test."""
    logger.info("")
    logger.info("=" * 60)
    logger.info("Cell type: %s", cell_type)
    logger.info("=" * 60)
    logger.info(
        "  Background (N):           %d genes", result["N_background"]
    )
    logger.info(
        "  MASLD genes in bg (K):    %d genes", result["K_masld_in_background"]
    )
    logger.info(
        "  Discovered genes (n):     %d genes", result["n_discovered"]
    )
    logger.info(
        "  Overlap with MASLD (k):   %d genes", result["k_overlap"]
    )
    logger.info(
        "  Expected overlap:         %.2f genes", result["expected_overlap"]
    )
    logger.info(
        "  Fold enrichment:          %.2fx", result["fold_enrichment"]
    )
    logger.info(
        "  Hypergeometric p-value:   %.4g", result["hypergeom_pvalue"]
    )
    logger.info(
        "  Fisher's exact p-value:   %.4g", result["fisher_pvalue"]
    )
    logger.info(
        "  Fisher's odds ratio:      %.2f", result["fisher_odds_ratio"]
    )

    if result["overlap_genes"]:
        logger.info(
            "  MASLD-associated genes found: %s",
            ", ".join(result["overlap_genes"]),
        )
    else:
        logger.info("  No MASLD-associated genes found.")

    non_masld = sorted(
        set(result["discovered_genes_list"]) - set(result["overlap_genes"])
    )
    if non_masld:
        logger.info("  Other discovered genes: %s", ", ".join(non_masld))


def main() -> None:
    """Run MASLD enrichment test for all cell types."""
    args = parse_arguments()

    if args.soft:
        input_dir = SOFT_PEAK2GENE_OUTPUT_DIR
        mode_label = "soft"
    else:
        input_dir = PEAK2GENE_OUTPUT_DIR
        mode_label = "default"

    if args.old:
        input_dir = PEAK2GENE_19_OUTPUT_DIR
        mode_label = "old"

    logger.info("Mode: %s", mode_label)
    logger.info("Input directory: %s", input_dir)

    # Fetch MASLD-associated genes from Open Targets Platform
    # Filter out genes with association score < 0.05 (minimal evidence threshold)
    masld_genes = fetch_masld_genes_from_open_targets(min_score=0.05)

    # Load background
    background_genes = load_background_genes(ABC_DICT_PATH)
    logger.info("ABC liver background: %d unique genes", len(background_genes))

    # Report MASLD gene coverage in background
    masld_in_bg = masld_genes & background_genes
    logger.info(
        "MASLD gene set (Open Targets): %d genes (%d present in background)",
        len(masld_genes),
        len(masld_in_bg),
    )

    # Per-cell-type tests
    all_discovered: set[str] = set()
    results_rows: list[dict[str, object]] = []

    for cell_type in SIGNIFICANT_PEAK_CELL_TYPES:
        peak2gene_path = input_dir / f"{cell_type}_peak2gene.csv"
        discovered = load_peak2gene_genes(peak2gene_path)
        all_discovered |= discovered

        result = run_hypergeometric_test(
            discovered_genes=discovered,
            background_genes=background_genes,
            masld_genes=masld_genes,
        )

        print_result(cell_type, result)

        results_rows.append({
            "cell_type": cell_type,
            "n_discovered": result["n_discovered"],
            "k_masld_overlap": result["k_overlap"],
            "expected_overlap": result["expected_overlap"],
            "fold_enrichment": result["fold_enrichment"],
            "hypergeom_pvalue": result["hypergeom_pvalue"],
            "fisher_pvalue": result["fisher_pvalue"],
            "odds_ratio": result["fisher_odds_ratio"],
            "masld_genes_found": "; ".join(result["overlap_genes"]),
            "all_genes_found": "; ".join(result["discovered_genes_list"]),
        })

    # Combined test across all cell types
    combined_result = run_hypergeometric_test(
        discovered_genes=all_discovered,
        background_genes=background_genes,
        masld_genes=masld_genes,
    )
    print_result("ALL_COMBINED", combined_result)

    results_rows.append({
        "cell_type": "ALL_COMBINED",
        "n_discovered": combined_result["n_discovered"],
        "k_masld_overlap": combined_result["k_overlap"],
        "expected_overlap": combined_result["expected_overlap"],
        "fold_enrichment": combined_result["fold_enrichment"],
        "hypergeom_pvalue": combined_result["hypergeom_pvalue"],
        "fisher_pvalue": combined_result["fisher_pvalue"],
        "odds_ratio": combined_result["fisher_odds_ratio"],
        "masld_genes_found": "; ".join(combined_result["overlap_genes"]),
        "all_genes_found": "; ".join(combined_result["discovered_genes_list"]),
    })

    # Save results
    output_path = OUTPUT_DIR / f"masld_enrichment_{mode_label}.csv"
    pd.DataFrame(results_rows).to_csv(output_path, index=False)
    logger.info("")
    logger.info("Results saved to: %s", output_path)

    # Final interpretation
    logger.info("")
    logger.info("=" * 60)
    logger.info("INTERPRETATION GUIDE")
    logger.info("=" * 60)
    logger.info(
        "  - Fold enrichment > 1 indicates our genes are enriched for MASLD genes"
    )
    logger.info(
        "  - p-value < 0.05 would indicate statistically significant enrichment"
    )
    logger.info(
        "  - With very few discovered genes (n < 10), the test has low power;"
    )
    logger.info(
        "    non-significance does NOT mean the genes are irrelevant to MASLD"
    )
    logger.info(
        "  - Individual gene-level evidence (literature) is more informative"
    )
    logger.info(
        "    at these sample sizes than formal enrichment statistics"
    )


if __name__ == "__main__":
    main()
