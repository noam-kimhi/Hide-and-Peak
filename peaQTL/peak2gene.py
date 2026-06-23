"""
peak2gene.py

Maps cell-type-specific significant peaks to target genes by intersecting
them with the Stanford ABC Liver reference dictionary via bioframe.overlap.
"""

import pandas as pd
import bioframe

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

ABC_PEAK2GENE_OUTPUT_DIRABC_DATA_DIR / "Stanford_ABC_Liver_Dictionary.csv"

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
    """Ensure chromosome names carry the 'chr' prefix robustly and drop bad headers."""
    if df.empty:
        return df
    
    df = df.copy()
    # Bulletproof 1: Drop accidental text headers (prevents the 10-row bug)
    df = df[df["chrom"].astype(str).str.lower() != "chrom"]
    
    # Bulletproof 2: Safely add 'chr' prefix only if it's missing (prevents 'chrchr1')
    df["chrom"] = df["chrom"].astype(str).apply(
        lambda x: x if x.startswith("chr") else f"chr{x}"
    )
    return df


def map_peaks_to_genes(
    peaks: pd.DataFrame,
    reference: pd.DataFrame,
) -> pd.DataFrame:
    """
    Intersect cell-type peaks with the ABC reference dictionary via bioframe.
    Contains dynamic column checking to be 100% immune to bioframe version changes.
    """
    peaks = _normalise_chroms(peaks)
    reference = _normalise_chroms(reference)

    logger.info(
        "    Chrom format check — peaks: %s  |  reference: %s",
        peaks["chrom"].iloc[0] if len(peaks) else "empty",
        reference["chrom"].iloc[0] if len(reference) else "empty",
    )

    # Force peaks to only contain the coordinate columns to avoid suffix collisions
    peaks = peaks[["chrom", "start", "end"]].copy()

    # Inner join via interval tree
    overlaps = bioframe.overlap(
        peaks,
        reference,
        cols1=["chrom", "start", "end"],
        cols2=["chrom", "start", "end"],
        how="inner",
        suffixes=("", "_ref"),
        return_index=False
    )

    # Bulletproof 3: Dynamic column resolution 
    # (handles both cases: whether bioframe appended '_ref' or not)
    target_col = "TargetGene_ref" if "TargetGene_ref" in overlaps.columns else "TargetGene"
    abc_col = "ABC.Score_ref" if "ABC.Score_ref" in overlaps.columns else "ABC.Score"

    # Safe extraction
    keep_cols = ["chrom", "start", "end", target_col, abc_col]
    result = overlaps[keep_cols].copy()

    # Final clean rename
    result = result.rename(columns={
        target_col: "TargetGene",
        abc_col: "ABC.Score"
    })

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