#!/usr/bin/env python3
"""Summarize DESeq2 differential-accessibility results across cell types."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import constants as C


REQUIRED_COLUMNS = {
    "peak_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}

NUMERIC_COLUMNS = [
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
]

# Median of a chi-square distribution with one degree of freedom.
# Used for the descriptive genomic-inflation statistic.
CHI_SQUARE_1DF_MEDIAN = 0.4549364231195727


def threshold_tag(value: float) -> str:
    """Convert a numeric threshold into a CSV-column-safe string."""
    formatted = np.format_float_positional(value, trim="-")
    return formatted.replace(".", "p")


def safe_fraction(numerator: int | float, denominator: int | float) -> float:
    """Return a fraction, or NaN when the denominator is zero."""
    if denominator == 0:
        return float("nan")
    return float(numerator) / float(denominator)


def finite_values(series: pd.Series) -> pd.Series:
    """Return finite numeric values from a Series."""
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric[np.isfinite(numeric)]


def valid_probabilities(series: pd.Series) -> pd.Series:
    """Return numeric probabilities, replacing invalid values with NaN."""
    numeric = pd.to_numeric(series, errors="coerce")
    valid = np.isfinite(numeric) & numeric.between(0.0, 1.0)
    return numeric.where(valid)


def quantile_or_nan(series: pd.Series, quantile: float) -> float:
    """Calculate a quantile, or NaN for an empty Series."""
    values = finite_values(series)

    if values.empty:
        return float("nan")

    return float(values.quantile(quantile))


def min_positive_or_nan(series: pd.Series) -> float:
    """Return the minimum strictly positive value, or NaN."""
    values = finite_values(series)
    values = values[values > 0]

    if values.empty:
        return float("nan")

    return float(values.min())


def load_deseq2_results(path: Path) -> pd.DataFrame:
    """Load and validate one DESeq2 result file."""
    dataframe = pd.read_csv(path)

    missing_columns = REQUIRED_COLUMNS.difference(dataframe.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"{path} is missing required columns: {missing}")

    dataframe = dataframe.copy()

    dataframe["peak_id"] = dataframe["peak_id"].astype("string")

    for column in NUMERIC_COLUMNS:
        dataframe[column] = pd.to_numeric(
            dataframe[column],
            errors="coerce",
        )

    return dataframe


def cell_type_from_path(path: Path) -> str:
    """Extract the cell type from a DESeq2 result filename."""
    if not path.name.endswith(C.DESEQ2_RESULTS_SUFFIX):
        raise ValueError(
            f"Unexpected filename {path.name!r}; expected suffix "
            f"{C.DESEQ2_RESULTS_SUFFIX!r}."
        )

    return path.name.removesuffix(C.DESEQ2_RESULTS_SUFFIX)


def discover_result_files() -> list[Path]:
    """Find all cell-type DESeq2 result files."""
    if not C.DESEQ2_RESULTS_DIR.exists():
        raise FileNotFoundError(
            f"DESeq2 result directory does not exist: "
            f"{C.DESEQ2_RESULTS_DIR}"
        )

    files = sorted(
        path
        for path in C.DESEQ2_RESULTS_DIR.glob(
            f"*{C.DESEQ2_RESULTS_SUFFIX}"
        )
        if path.is_file()
    )

    if not files:
        raise FileNotFoundError(
            f"No files ending with {C.DESEQ2_RESULTS_SUFFIX!r} were found "
            f"inside {C.DESEQ2_RESULTS_DIR}."
        )

    return files


def add_distribution_statistics(
    summary: dict[str, Any],
    prefix: str,
    values: pd.Series,
) -> None:
    """Add robust distribution statistics to a summary dictionary."""
    finite = finite_values(values)

    summary[f"{prefix}_min"] = (
        float(finite.min()) if not finite.empty else float("nan")
    )
    summary[f"{prefix}_q001"] = quantile_or_nan(finite, 0.001)
    summary[f"{prefix}_q01"] = quantile_or_nan(finite, 0.01)
    summary[f"{prefix}_q05"] = quantile_or_nan(finite, 0.05)
    summary[f"{prefix}_q25"] = quantile_or_nan(finite, 0.25)
    summary[f"{prefix}_median"] = quantile_or_nan(finite, 0.50)
    summary[f"{prefix}_q75"] = quantile_or_nan(finite, 0.75)
    summary[f"{prefix}_q95"] = quantile_or_nan(finite, 0.95)
    summary[f"{prefix}_q99"] = quantile_or_nan(finite, 0.99)
    summary[f"{prefix}_max"] = (
        float(finite.max()) if not finite.empty else float("nan")
    )


def summarize_cell_type(
    path: Path,
    dataframe: pd.DataFrame,
) -> dict[str, Any]:
    """Create one wide summary row for one cell type."""
    cell_type = cell_type_from_path(path)

    pvalue = valid_probabilities(dataframe["pvalue"])
    padj = valid_probabilities(dataframe["padj"])

    base_mean = finite_values(dataframe["baseMean"])
    log2fc = pd.to_numeric(
        dataframe["log2FoldChange"],
        errors="coerce",
    )
    lfc_se = pd.to_numeric(dataframe["lfcSE"], errors="coerce")
    statistic = pd.to_numeric(dataframe["stat"], errors="coerce")

    finite_log2fc = log2fc.where(np.isfinite(log2fc))
    finite_lfc_se = lfc_se.where(np.isfinite(lfc_se))
    finite_statistic = statistic.where(np.isfinite(statistic))

    absolute_log2fc = finite_log2fc.abs()

    n_total = len(dataframe)
    n_pvalue_available = int(pvalue.notna().sum())
    n_padj_available = int(padj.notna().sum())
    n_testable_with_effect = int(
        (padj.notna() & finite_log2fc.notna()).sum()
    )

    summary: dict[str, Any] = {
        "cell_type": cell_type,
        "input_file": str(path),
        "n_peaks_total": n_total,
        "n_unique_peak_ids": int(dataframe["peak_id"].nunique(dropna=True)),
        "n_duplicate_peak_ids": int(dataframe["peak_id"].duplicated().sum()),
        "n_rows_complete_numeric": int(
            dataframe[NUMERIC_COLUMNS].notna().all(axis=1).sum()
        ),
        "n_pvalue_available": n_pvalue_available,
        "fraction_pvalue_missing_or_invalid": safe_fraction(
            n_total - n_pvalue_available,
            n_total,
        ),
        "n_padj_available": n_padj_available,
        "fraction_padj_missing_or_invalid": safe_fraction(
            n_total - n_padj_available,
            n_total,
        ),
        "n_padj_and_log2fc_available": n_testable_with_effect,
        "n_zero_pvalues": int((pvalue == 0).sum()),
        "n_zero_padj": int((padj == 0).sum()),
        "min_nonzero_pvalue": min_positive_or_nan(pvalue),
        "min_nonzero_padj": min_positive_or_nan(padj),
        "n_positive_log2fc": int((finite_log2fc > 0).sum()),
        "n_negative_log2fc": int((finite_log2fc < 0).sum()),
        "n_zero_log2fc": int((finite_log2fc == 0).sum()),
        "fraction_positive_log2fc": safe_fraction(
            int((finite_log2fc > 0).sum()),
            int(finite_log2fc.notna().sum()),
        ),
        "fraction_negative_log2fc": safe_fraction(
            int((finite_log2fc < 0).sum()),
            int(finite_log2fc.notna().sum()),
        ),
    }

    add_distribution_statistics(summary, "pvalue", pvalue)
    add_distribution_statistics(summary, "padj", padj)
    add_distribution_statistics(summary, "baseMean", base_mean)
    add_distribution_statistics(summary, "log2FoldChange", finite_log2fc)
    add_distribution_statistics(summary, "abs_log2FoldChange", absolute_log2fc)
    add_distribution_statistics(summary, "lfcSE", finite_lfc_se)
    add_distribution_statistics(summary, "stat", finite_statistic)

    # Raw-P-value summaries and enrichment over the uniform null expectation.
    for threshold in C.DESEQ2_PVALUE_THRESHOLDS:
        tag = threshold_tag(threshold)
        count = int((pvalue < threshold).sum())
        expected_under_uniform = n_pvalue_available * threshold

        summary[f"n_pvalue_lt_{tag}"] = count
        summary[f"fraction_pvalue_lt_{tag}_among_available"] = safe_fraction(
            count,
            n_pvalue_available,
        )
        summary[f"fraction_pvalue_lt_{tag}_among_all_peaks"] = safe_fraction(
            count,
            n_total,
        )
        summary[
            f"pvalue_lt_{tag}_enrichment_over_uniform"
        ] = safe_fraction(
            count,
            expected_under_uniform,
        )

    # Adjusted-P-value/FDR summaries.
    for threshold in C.DESEQ2_PADJ_THRESHOLDS:
        tag = threshold_tag(threshold)
        count = int((padj < threshold).sum())

        summary[f"n_padj_lt_{tag}"] = count
        summary[f"fraction_padj_lt_{tag}_among_available"] = safe_fraction(
            count,
            n_padj_available,
        )
        summary[f"fraction_padj_lt_{tag}_among_all_peaks"] = safe_fraction(
            count,
            n_total,
        )

    # Low-count/baseMean diagnostics.
    for threshold in C.DESEQ2_BASE_MEAN_THRESHOLDS:
        tag = threshold_tag(threshold)
        count = int((dataframe["baseMean"] < threshold).sum())

        summary[f"n_baseMean_lt_{tag}"] = count
        summary[f"fraction_baseMean_lt_{tag}"] = safe_fraction(
            count,
            n_total,
        )

    # How many nominally significant peaks disappear after correction?
    nominal_p_mask = pvalue < 0.05
    fdr_005_mask = padj < 0.05

    summary["n_pvalue_lt_0p05_but_not_padj_lt_0p05"] = int(
        (nominal_p_mask & ~fdr_005_mask).sum()
    )

    # Default display/candidate criterion.
    default_fdr_mask = padj < C.DESEQ2_DEFAULT_PADJ_THRESHOLD
    default_effect_mask = (
        absolute_log2fc >= C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
    )
    default_significant_mask = default_fdr_mask & default_effect_mask

    summary["default_padj_threshold"] = (
        C.DESEQ2_DEFAULT_PADJ_THRESHOLD
    )
    summary["default_abs_log2fc_threshold"] = (
        C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
    )
    summary["default_n_fdr_significant_without_effect_filter"] = int(
        default_fdr_mask.sum()
    )
    summary["default_n_fdr_significant_below_effect_threshold"] = int(
        (default_fdr_mask & ~default_effect_mask).sum()
    )
    summary["default_n_significant"] = int(
        default_significant_mask.sum()
    )
    summary["default_n_significant_positive_log2fc"] = int(
        (default_significant_mask & (finite_log2fc > 0)).sum()
    )
    summary["default_n_significant_negative_log2fc"] = int(
        (default_significant_mask & (finite_log2fc < 0)).sum()
    )
    summary["default_fraction_significant_among_all_peaks"] = safe_fraction(
        int(default_significant_mask.sum()),
        n_total,
    )
    summary[
        "default_fraction_significant_among_testable_peaks"
    ] = safe_fraction(
        int(default_significant_mask.sum()),
        n_testable_with_effect,
    )

    significant_abs_lfc = absolute_log2fc[default_significant_mask]
    significant_base_mean = dataframe.loc[
        default_significant_mask,
        "baseMean",
    ]

    summary["default_significant_median_abs_log2fc"] = (
        float(significant_abs_lfc.median())
        if not significant_abs_lfc.empty
        else float("nan")
    )
    summary["default_significant_median_baseMean"] = (
        float(significant_base_mean.median())
        if not significant_base_mean.empty
        else float("nan")
    )

    # Descriptive lambda-GC based on the DESeq2 Wald statistic.
    # This is diagnostic only: high values can reflect either genuine signal
    # or model/confounding problems.
    valid_statistic = finite_values(finite_statistic)
    if valid_statistic.empty:
        summary["lambda_gc"] = float("nan")
    else:
        summary["lambda_gc"] = float(
            np.median(np.square(valid_statistic))
            / CHI_SQUARE_1DF_MEDIAN
        )

    return summary


def create_threshold_grid(
    path: Path,
    dataframe: pd.DataFrame,
) -> list[dict[str, Any]]:
    """Count significant peaks for every configured threshold combination."""
    cell_type = cell_type_from_path(path)

    padj = valid_probabilities(dataframe["padj"])
    log2fc = pd.to_numeric(
        dataframe["log2FoldChange"],
        errors="coerce",
    ).where(
        lambda values: np.isfinite(values)
    )

    absolute_log2fc = log2fc.abs()

    n_total = len(dataframe)
    testable_mask = padj.notna() & log2fc.notna()
    n_testable = int(testable_mask.sum())

    rows: list[dict[str, Any]] = []

    for padj_threshold in C.DESEQ2_PADJ_THRESHOLDS:
        for abs_log2fc_threshold in C.DESEQ2_ABS_LOG2FC_THRESHOLDS:
            significant_mask = (
                testable_mask
                & (padj < padj_threshold)
                & (absolute_log2fc >= abs_log2fc_threshold)
            )

            n_significant = int(significant_mask.sum())
            n_positive = int(
                (significant_mask & (log2fc > 0)).sum()
            )
            n_negative = int(
                (significant_mask & (log2fc < 0)).sum()
            )

            rows.append(
                {
                    "cell_type": cell_type,
                    "input_file": str(path),
                    "padj_threshold": padj_threshold,
                    "abs_log2fc_threshold": abs_log2fc_threshold,
                    "n_peaks_total": n_total,
                    "n_testable_peaks": n_testable,
                    "n_significant": n_significant,
                    "n_positive_log2fc": n_positive,
                    "n_negative_log2fc": n_negative,
                    "fraction_of_all_peaks": safe_fraction(
                        n_significant,
                        n_total,
                    ),
                    "fraction_of_testable_peaks": safe_fraction(
                        n_significant,
                        n_testable,
                    ),
                }
            )

    return rows


def main() -> None:
    """Run the cross-cell-type analysis."""
    files = discover_result_files()

    summary_rows: list[dict[str, Any]] = []
    threshold_rows: list[dict[str, Any]] = []

    for path in files:
        print(f"Analyzing {path.name}...")
        dataframe = load_deseq2_results(path)

        summary_rows.append(
            summarize_cell_type(path, dataframe)
        )
        threshold_rows.extend(
            create_threshold_grid(path, dataframe)
        )

    summary_dataframe = (
        pd.DataFrame(summary_rows)
        .sort_values("cell_type")
        .reset_index(drop=True)
    )

    threshold_dataframe = (
        pd.DataFrame(threshold_rows)
        .sort_values(
            [
                "cell_type",
                "padj_threshold",
                "abs_log2fc_threshold",
            ],
            ascending=[True, False, True],
        )
        .reset_index(drop=True)
    )

    C.DESEQ2_ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    summary_dataframe.to_csv(
        C.DESEQ2_SUMMARY_CSV,
        index=False,
        float_format="%.10g",
    )
    threshold_dataframe.to_csv(
        C.DESEQ2_THRESHOLD_GRID_CSV,
        index=False,
        float_format="%.10g",
    )

    print(
        f"Saved cell-type summary to:\n"
        f"  {C.DESEQ2_SUMMARY_CSV}"
    )
    print(
        f"Saved threshold sensitivity table to:\n"
        f"  {C.DESEQ2_THRESHOLD_GRID_CSV}"
    )


if __name__ == "__main__":
    main()