
"""Differential chromatin accessibility analysis with PyDESeq2.

For each cell type this script:

1. Loads the pseudobulk count table from data/derived/pseudobulk/.
2. Drops excluded samples (Normal_rep4 and Normal_rep6, which lack
   meaningful data).
3. Pre-filters peaks that are nearly always zero to reduce the multiple-
   testing burden and speed up computation.
4. Runs PyDESeq2 (MASH vs. Normal contrast).
5. Writes per-cell-type results to
   results/peaQTL/differential_peaks/{cell_type}_deseq2_results.csv

Interpretation of output columns
---------------------------------
baseMean     : mean normalised count across all retained samples
log2FoldChange: log2(MASH / Normal) — positive → more open in MASH
lfcSE        : standard error of the log2 fold-change
stat         : Wald test statistic
pvalue       : raw p-value
padj         : Benjamini–Hochberg adjusted p-value
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from constants import (
    CELL_TYPE_STANDARDIZATION,
    PSEUDOBULK_OUTPUT_DIR,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Samples to exclude before any analysis.
EXCLUDED_SAMPLES: list[str] = ["Normal_rep4", "Normal_rep6"]

# Pre-filter: keep a peak only if at least MIN_SAMPLES_NONZERO samples have
# a count strictly greater than zero.  Removes peaks that are almost always
# zero across the retained samples, which would never reach significance and
# inflate the multiple-testing correction.
MIN_SAMPLES_NONZERO: int = 2

OUTPUT_DIR: Path = (
    Path(__file__).resolve().parent.parent
    / "results"
    / "peaQTL"
    / "differential_peaks"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_pseudobulk(cell_type: str) -> pd.DataFrame:
    """
    Load the pseudobulk count table for *cell_type*.

    :return: DataFrame with peaks as rows and samples as columns (raw counts).
    """
    path = PSEUDOBULK_OUTPUT_DIR / f"{cell_type}_pseudobulk.csv"
    return pd.read_csv(path, index_col="peak_id")


def drop_excluded_samples(df: pd.DataFrame) -> pd.DataFrame:
    """Drop EXCLUDED_SAMPLES columns; silently skip any that are absent."""
    to_drop = [c for c in EXCLUDED_SAMPLES if c in df.columns]
    if to_drop:
        print(f"  Dropping excluded samples: {to_drop}")
    return df.drop(columns=to_drop)


def build_metadata(sample_names: list[str]) -> pd.DataFrame:
    """
    Build a sample × condition metadata DataFrame from column names.

    Column names follow the pattern '<Condition>_rep<N>'.  The first
    underscore-delimited token is used as the condition label.

    :param sample_names: Ordered list of retained sample names.
    :return: DataFrame indexed by sample name with a 'condition' column.
    """
    conditions = [name.split("_")[0] for name in sample_names]
    return pd.DataFrame(
        {"condition": conditions},
        index=pd.Index(sample_names, name="sample"),
    )


def filter_low_count_peaks(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Remove peaks where fewer than MIN_SAMPLES_NONZERO samples have count > 0.

    :param counts: Samples × peaks integer DataFrame.
    :return: Filtered DataFrame (same orientation).
    """
    n_nonzero = (counts > 0).sum(axis=0)  # per peak
    keep = n_nonzero >= MIN_SAMPLES_NONZERO
    kept = int(keep.sum())
    print(
        f"  Peak filter (≥{MIN_SAMPLES_NONZERO} non-zero samples): "
        f"{kept:,} / {len(keep):,} peaks retained"
    )
    return counts.loc[:, keep]


def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    cell_type: str,
) -> pd.DataFrame:
    """
    Run PyDESeq2 on *counts* and return the results DataFrame.

    :param counts: Samples × peaks integer DataFrame (already filtered).
    :param metadata: Samples × factors DataFrame with a 'condition' column.
    :param cell_type: Label used only for progress messages.
    :return: Results DataFrame indexed by peak_id.
    """
    print(f"  Fitting DESeq2 model ({counts.shape[0]} samples, "
          f"{counts.shape[1]:,} peaks)…")

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        quiet=False,
    )
    dds.deseq2()

    print("  Computing Wald test statistics (MASH vs. Normal)…")
    stat_res = DeseqStats(
        dds,
        contrast=["condition", "MASH", "Normal"],
        quiet=False,
    )
    stat_res.summary()

    results = stat_res.results_df.copy()
    results.index.name = "peak_id"
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cell_types = sorted(set(CELL_TYPE_STANDARDIZATION.values()))
    summary_records: list[dict] = []

    for cell_type in cell_types:
        print(f"\n{'='*60}")
        print(f"Cell type: {cell_type}")
        print(f"{'='*60}")

        # 1. Load
        csv_path = PSEUDOBULK_OUTPUT_DIR / f"{cell_type}_pseudobulk.csv"
        if not csv_path.exists():
            print(f"  [SKIP] {csv_path.name} not found.")
            continue

        df = load_pseudobulk(cell_type)      # peaks × samples
        df = drop_excluded_samples(df)

        n_samples = df.shape[1]
        print(f"  Retained samples ({n_samples}): {list(df.columns)}")

        # 2. Transpose to samples × peaks (PyDESeq2 convention)
        counts = df.T.astype(int)            # samples × peaks

        # 3. Build metadata
        metadata = build_metadata(list(counts.index))
        condition_counts = metadata["condition"].value_counts().to_dict()
        print(f"  Condition counts: {condition_counts}")

        # 4. Pre-filter
        counts = filter_low_count_peaks(counts)

        # 5. Run DESeq2
        results = run_deseq2(counts, metadata, cell_type)

        # 6. Save
        out_path = OUTPUT_DIR / f"{cell_type}_deseq2_results.csv"
        results.to_csv(out_path)

        n_sig_05 = int((results["padj"] < 0.05).sum())
        n_sig_01 = int((results["padj"] < 0.01).sum())
        n_up     = int(((results["padj"] < 0.05) & (results["log2FoldChange"] > 0)).sum())
        n_down   = int(((results["padj"] < 0.05) & (results["log2FoldChange"] < 0)).sum())

        print(
            f"\n  Results saved → {out_path.name}\n"
            f"  Significant peaks (padj < 0.05) : {n_sig_05:,}  "
            f"(↑MASH: {n_up:,}  ↓MASH: {n_down:,})\n"
            f"  Significant peaks (padj < 0.01) : {n_sig_01:,}"
        )

        summary_records.append({
            "cell_type"        : cell_type,
            "n_samples"        : n_samples,
            "n_peaks_tested"   : len(results),
            "n_sig_padj_05"    : n_sig_05,
            "n_sig_padj_01"    : n_sig_01,
            "n_up_in_MASH"     : n_up,
            "n_down_in_MASH"   : n_down,
        })

    # Cross-cell-type summary
    if summary_records:
        summary_df = pd.DataFrame(summary_records)
        summary_path = OUTPUT_DIR / "differential_peaks_summary.csv"
        summary_df.to_csv(summary_path, index=False)

        print(f"\n{'='*60}")
        print("Summary across cell types")
        print(f"{'='*60}")
        print(summary_df.to_string(index=False))
        print(f"\nSummary saved → {summary_path}")


if __name__ == "__main__":
    main()
