"""
peak2gene.py

Maps cell-type-specific significant peaks to target genes by intersecting
them with the Stanford ABC Liver reference dictionary via bioframe.overlap.
"""

import logging
import pandas as pd

from constants import (
    ABC_DATA_DIR,
    FILTERED_SIG_PEAKS_DIR,
    FILTERED_SIGNIFICANT_PEAKS_SUFFIX,
    SIGNIFICANT_PEAK_CELL_TYPES,
    PEAQTL_RESULTS_DIR,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ABC_REFERENCE_PATH = ABC_DATA_DIR / "Stanford_ABC_Liver_Dictionary.csv"

PEAK2GENE_OUTPUT_DIR = PEAQTL_RESULTS_DIR / "peak2gene"
PEAK2GENE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_abc_reference(path) -> pd.DataFrame:
    """Load the pre-filtered ABC liver reference dictionary."""
    logger.info("Loading ABC reference dictionary from: %s", path)
    df = pd.read_csv(path, dtype={"chrom": str, "start": "Int64", "end": "Int64"})
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    logger.info("Reference dictionary loaded: %d rows.", len(df))
    return df


def load_cell_type_peaks(cell_type: str) -> pd.DataFrame:
    """Load the filtered significant peaks BED file for a single cell type."""
    path = FILTERED_SIG_PEAKS_DIR / f"{cell_type}{FILTERED_SIGNIFICANT_PEAKS_SUFFIX}"
    logger.info("  Loading peaks from: %s", path)
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        header=None,
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": int, "end": int},
    )
    logger.info("  Peaks loaded: %d rows.", len(df))
    return df


def _normalise_chroms(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure chromosome names carry the 'chr' prefix (e.g. '1' -> 'chr1')."""
    if df.empty:
        return df
    if not df["chrom"].iloc[0].startswith("chr"):
        df = df.copy()
        df["chrom"] = "chr" + df["chrom"].astype(str)
    return df


def map_peaks_to_genes(
    peaks: pd.DataFrame,
    reference: pd.DataFrame,
) -> pd.DataFrame:
    """
    Intersect cell-type peaks with the ABC reference dictionary.

    Uses a pure pandas merge approach to avoid any dependency on bioframe's
    column-suffix behaviour:
      1. Merge left (peaks) and right (reference) on 'chrom' — equality join
         per chromosome.
      2. Filter merged rows to keep only pairs where the intervals overlap:
             peak.start < ref.end  AND  ref.start < peak.end
         (standard half-open interval overlap condition).
      3. Return only the original peak coordinates plus annotation columns.
    """
    peaks = _normalise_chroms(peaks)
    reference = _normalise_chroms(reference)

    logger.info(
        "    Chrom format check — peaks: %s  |  reference: %s",
        peaks["chrom"].iloc[0] if len(peaks) else "empty",
        reference["chrom"].iloc[0] if len(reference) else "empty",
    )

    # Merge on chromosome (reduces search space before interval filtering).
    merged = peaks.merge(
        reference[["chrom", "start", "end", "TargetGene", "ABC.Score"]],
        on="chrom",
        suffixes=("", "_ref"),
    )

    # Keep only rows where the two intervals actually overlap.
    # Half-open interval overlap: [a, b) overlaps [c, d) iff a < d AND c < b
    overlap_mask = (merged["start"] < merged["end_ref"]) & (merged["start_ref"] < merged["end"])
    result = (
        merged.loc[overlap_mask, ["chrom", "start", "end", "TargetGene", "ABC.Score"]]
        .reset_index(drop=True)
    )

    return result


def run() -> None:
    reference = load_abc_reference(ABC_REFERENCE_PATH)

    for cell_type in SIGNIFICANT_PEAK_CELL_TYPES:
        logger.info("Processing cell type: %s", cell_type)

        peaks = load_cell_type_peaks(cell_type)
        n_input_peaks = peaks["chrom"].count()

        mapped = map_peaks_to_genes(peaks, reference)

        n_mapped_peaks = mapped[["chrom", "start", "end"]].drop_duplicates().shape[0]
        logger.info(
            "  Mapped: %d / %d peaks matched at least one gene.",
            n_mapped_peaks,
            n_input_peaks,
        )

        out_path = PEAK2GENE_OUTPUT_DIR / f"{cell_type}_peak2gene.csv"
        mapped.to_csv(out_path, index=False)
        logger.info("  Saved: %s (%d rows)\n", out_path, len(mapped))


if __name__ == "__main__":
    run()
