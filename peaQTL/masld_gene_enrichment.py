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
from scipy.stats import hypergeom, fisher_exact

from constants import (
    ABC_DICT_PATH,
    PEAK2GENE_OUTPUT_DIR,
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
# Curated MASLD/NASH-associated gene set
# Sources: GWAS (Anstee et al. 2020, Romeo et al. 2008, Abul-Husn et al. 2018),
#          DisGeNET (C0400966 - NAFLD), KEGG NAFLD pathway (hsa04932),
#          key functional studies in MASLD pathogenesis.
# ---------------------------------------------------------------------------

MASLD_GENES: Final[frozenset[str]] = frozenset({
    # --- GWAS-validated risk loci ---
    "PNPLA3", "TM6SF2", "MBOAT7", "GCKR", "HSD17B13",
    "MARC1", "MTARC1", "GPAM", "APOB", "MTTP",

    # --- Lipid metabolism / de novo lipogenesis ---
    "SREBF1", "SREBF2", "FASN", "SCD", "ACACA", "ACACB",
    "DGAT1", "DGAT2", "PPARA", "PPARG", "PPARGC1A",
    "ACOX1", "CPT1A", "CPT2", "HMGCR", "LDLR",
    "ABCA1", "PCSK9", "ANGPTL3", "ANGPTL4",
    "LPIN1", "LPIN2", "AGPAT2",

    # --- Insulin signaling / glucose metabolism ---
    "INSR", "IRS1", "IRS2", "PIK3CA", "AKT1", "AKT2",
    "FOXO1", "GSK3B", "GCK", "G6PC", "PCK1", "PCK2",
    "ADIPOQ", "ADIPOR1", "ADIPOR2", "LEP", "LEPR",
    "IGFBP1", "IGFBPL1", "IGF1", "IGF1R",

    # --- Inflammation / innate immunity ---
    "TNF", "TNFRSF1A", "IL6", "IL6R", "IL1B", "IL1R1",
    "IL18", "NLRP3", "CASP1", "CCL2", "CCR2",
    "TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4",
    "NFKB1", "RELA", "IKBKB", "TOLLIP", "TLR4", "TLR2",
    "CD14", "CD68", "MARCO", "CLEC4F",

    # --- Oxidative stress / mitochondrial dysfunction ---
    "SOD1", "SOD2", "CAT", "GPX1", "NRF2", "NFE2L2", "KEAP1",
    "OPA1", "MFN1", "MFN2", "DNM1L", "FIS1",
    "CYCS", "BAX", "BCL2", "CASP3", "CASP9",
    "CYP2E1", "CYP4A11",

    # --- Fibrosis / stellate cell activation ---
    "COL1A1", "COL1A2", "COL3A1", "COL4A1",
    "ACTA2", "TIMP1", "TIMP2", "MMP2", "MMP9", "MMP13",
    "LOX", "LOXL2", "CTGF", "CCN2",
    "PDGFA", "PDGFB", "PDGFRA", "PDGFRB",
    "WNT2", "WNT3A", "CTNNB1",
    "PIEZO1", "PIEZO2",

    # --- Notch / Hedgehog signaling (biliary/Kupffer) ---
    "NOTCH1", "NOTCH2", "JAG1", "JAG2", "DLL1", "DLL4",
    "HES1", "HEY1", "RBPJ",
    "SHH", "IHH", "SMO", "GLI1", "GLI2",

    # --- Apoptosis / cell death (hepatocyte ballooning) ---
    "RIPK1", "RIPK3", "MLKL", "FADD", "FAS", "TRAIL",
    "PARP1", "HMGB1",

    # --- Bile acid metabolism ---
    "CYP7A1", "CYP27A1", "NR1H4", "FXR", "ABCB11",
    "SLC10A1", "NTCP",

    # --- Vascular remodeling / angiogenesis ---
    "VEGFA", "VEGFB", "KDR", "FLT1",
    "PECAM1", "VWF", "NOS3", "EDN1",
    "SEMA3C", "SEMA3A", "NRP1", "PLXNA1",

    # --- Iron metabolism (relevant in NASH) ---
    "HFE", "TFR2", "HAMP", "SLC40A1", "FTH1", "FTL",

    # --- Epigenetic regulators implicated in MASLD ---
    "HDAC1", "HDAC3", "SIRT1", "DNMT1", "TET2",

    # --- Platelet / coagulation (NASH-associated) ---
    "GP5", "GP1BA", "THBS1", "F2", "SERPINE1",

    # --- Amino acid / one-carbon metabolism ---
    "MTHFR", "MAT1A", "BHMT", "CBS", "GNMT",

    # --- Autophagy ---
    "BECN1", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1",

    # --- Other validated MASLD-associated genes ---
    "SAMM50", "ERLIN1", "LYPLAL1", "PPP1R3B",
    "SOX9", "KRT19", "KRT7", "EPCAM",
    "CRP", "SAA1", "HP", "ORM1",
    "AOAH", "RND1",
})


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

    logger.info("Mode: %s", mode_label)
    logger.info("Input directory: %s", input_dir)

    # Load background
    background_genes = load_background_genes(ABC_DICT_PATH)
    logger.info("ABC liver background: %d unique genes", len(background_genes))

    # Report MASLD gene coverage in background
    masld_in_bg = MASLD_GENES & background_genes
    logger.info(
        "MASLD curated set: %d genes (%d present in background)",
        len(MASLD_GENES),
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
            masld_genes=MASLD_GENES,
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
        masld_genes=MASLD_GENES,
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
