"""Analyze candidate hybrid peak-filtering hyperparameters.

This script does not filter or overwrite any matrix. It reads the binarized
sample-by-cell-type matrices, calculates per-peak cell support, evaluates a grid
of candidate hybrid filtering rules, and writes statistics and plots beneath:

    results/preprocessing/peak_filtering/

For a matrix with ``n`` cells, a peak passes a hybrid rule when:

    support >= max(min_cell_support, ceil(min_cell_fraction * n))

A matrix is included in a candidate analysis only when:

    n >= min_group_size

The script evaluates local sample-by-cell-type support only. Replicate-level
support cannot be calibrated correctly until sample-specific peaks have been
projected onto a common consensus-peak coordinate system.
"""

from __future__ import annotations

import gzip
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final, TypeAlias

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from constants import (
    ATAC_SEQ_DIRS,
    BINARIZED_CELL_TYPE_MATRICES_DIR,
    CELL_TYPE_STANDARDIZATION,
    MATRIX_SUFFIX,
    PEAK_FILTERING_RESULTS_DIR,
    PEAK_FILTER_MIN_CELL_FRACTIONS,
    PEAK_FILTER_MIN_CELL_SUPPORTS,
    PEAK_FILTER_MIN_GROUP_SIZES,
    PLOT_DPI,
    PLOTS_MAIN_COLOR,
    PLOTS_SECOND_COLOR,
)


PLOTS_DIR: Final[Path] = PEAK_FILTERING_RESULTS_DIR / "plots"
OVERALL_HEATMAPS_DIR: Final[Path] = PLOTS_DIR / "retained_fraction_overall"
CELL_TYPE_HEATMAPS_DIR: Final[Path] = PLOTS_DIR / "retained_fraction_by_cell_type"
THRESHOLD_CURVES_DIR: Final[Path] = PLOTS_DIR / "threshold_curves"

PlotColor: TypeAlias = str | tuple[float, float, float, float]
SCRIPT_VERSION: Final[str] = "distinct-colors-v2"
HEATMAP_COLORMAP: Final[str] = "viridis"

GROUP_SIZES_PATH: Final[Path] = PEAK_FILTERING_RESULTS_DIR / "group_sizes.csv"
GROUP_SIZE_SENSITIVITY_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "minimum_group_size_sensitivity.csv"
)
PER_MATRIX_SUPPORT_SUMMARY_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "per_matrix_peak_support_summary.csv"
)
SUPPORT_DISTRIBUTIONS_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "peak_support_count_distributions.csv.gz"
)
HYPERPARAMETER_GRID_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "hyperparameter_grid.csv"
)
HYBRID_PER_MATRIX_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "hybrid_threshold_per_matrix.csv.gz"
)
HYBRID_OVERALL_SUMMARY_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "hybrid_threshold_summary_overall.csv"
)
HYBRID_CONDITION_SUMMARY_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "hybrid_threshold_summary_by_condition.csv"
)
HYBRID_CELL_TYPE_SUMMARY_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "hybrid_threshold_summary_by_cell_type.csv"
)
HYBRID_CELL_TYPE_CONDITION_SUMMARY_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR
    / "hybrid_threshold_summary_by_cell_type_and_condition.csv"
)
MANIFEST_PATH: Final[Path] = (
    PEAK_FILTERING_RESULTS_DIR / "analysis_manifest.txt"
)

PREFERRED_CELL_TYPE_ORDER: Final[tuple[str, ...]] = (
    "Hepatocytes",
    "Endothelial",
    "Cholangiocyte",
    "Kupffer",
    "Stellate",
    "T_NK_B",
    "Unknown",
)
CONDITION_ORDER: Final[tuple[str, ...]] = ("MASH", "Normal")
SUPPORT_CATEGORY_ORDER: Final[tuple[str, ...]] = (
    "0",
    "1",
    "2",
    "3-4",
    ">=5",
)


def validate_hyperparameter_values(
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
    minimum_cell_fractions: Sequence[float],
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[float, ...]]:
    """
    Validate and sort all candidate hyperparameter values.

    :param minimum_group_sizes: Candidate minimum numbers of cells required for
        a sample-by-cell-type matrix to be considered usable.
    :param minimum_cell_supports: Candidate minimum absolute peak-support
        counts.
    :param minimum_cell_fractions: Candidate minimum fractions of cells in
        which a peak must be accessible.
    :return: Sorted, unique tuples containing the validated values.
    """
    group_sizes = tuple(sorted({int(value) for value in minimum_group_sizes}))
    supports = tuple(sorted({int(value) for value in minimum_cell_supports}))
    fractions = tuple(
        sorted({float(value) for value in minimum_cell_fractions})
    )

    if not group_sizes or any(value < 1 for value in group_sizes):
        raise ValueError(
            "PEAK_FILTER_MIN_GROUP_SIZES must contain positive integers."
        )

    if not supports or any(value < 1 for value in supports):
        raise ValueError(
            "PEAK_FILTER_MIN_CELL_SUPPORTS must contain positive integers."
        )

    if not fractions or any(value < 0.0 or value > 1.0 for value in fractions):
        raise ValueError(
            "PEAK_FILTER_MIN_CELL_FRACTIONS must contain values in [0, 1]."
        )

    return group_sizes, supports, fractions


def get_cell_type_order() -> tuple[str, ...]:
    """
    Determine and validate the expected standardized cell-type order.

    :return: Ordered standardized cell-type names.
    """
    standardized_cell_types = tuple(
        dict.fromkeys(CELL_TYPE_STANDARDIZATION.values())
    )

    if set(standardized_cell_types) != set(PREFERRED_CELL_TYPE_ORDER):
        raise ValueError(
            "CELL_TYPE_STANDARDIZATION contains an unexpected set of "
            f"standardized cell types: {standardized_cell_types}"
        )

    return PREFERRED_CELL_TYPE_ORDER


def infer_condition(sample_name: str) -> str:
    """
    Infer the biological condition from a sample-directory name.

    :param sample_name: Sample-directory name.
    :return: Either ``MASH`` or ``Normal``.
    """
    lowered_name = sample_name.lower()

    if "mash" in lowered_name:
        return "MASH"

    if "normal" in lowered_name:
        return "Normal"

    raise ValueError(
        f"Could not infer condition from sample name: {sample_name}"
    )


def index_binary_matrices(
    cell_type_order: Sequence[str],
) -> dict[tuple[str, str], Path]:
    """
    Index every available binarized matrix by sample and cell type.

    :param cell_type_order: Expected standardized cell-type names.
    :return: Mapping from ``(sample_directory, cell_type)`` to matrix path.
    """
    known_samples = {sample_directory.name for sample_directory in ATAC_SEQ_DIRS}
    known_cell_types = set(cell_type_order)
    indexed_matrices: dict[tuple[str, str], Path] = {}

    matrix_paths = sorted(
        BINARIZED_CELL_TYPE_MATRICES_DIR.rglob(f"*{MATRIX_SUFFIX}")
    )

    if not matrix_paths:
        raise FileNotFoundError(
            "No binarized matrices were found beneath "
            f"{BINARIZED_CELL_TYPE_MATRICES_DIR}."
        )

    for matrix_path in matrix_paths:
        relative_path = matrix_path.relative_to(
            BINARIZED_CELL_TYPE_MATRICES_DIR
        )

        if len(relative_path.parts) < 3:
            raise ValueError(
                f"Unexpected matrix path structure: {matrix_path}"
            )

        sample_name = relative_path.parts[0]
        cell_type = relative_path.parts[1]

        if sample_name not in known_samples:
            raise ValueError(
                f"Unexpected sample directory beneath the binarized root: "
                f"{sample_name}"
            )

        if cell_type not in known_cell_types:
            raise ValueError(
                f"Unexpected cell type beneath {sample_name}: {cell_type}"
            )

        key = (sample_name, cell_type)

        if key in indexed_matrices:
            raise RuntimeError(
                f"More than one matrix was found for {sample_name} / "
                f"{cell_type}."
            )

        indexed_matrices[key] = matrix_path

    return indexed_matrices


def read_binary_matrix(matrix_path: Path) -> sparse.csr_matrix:
    """
    Read and validate one gzip-compressed binarized Matrix Market file.

    :param matrix_path: Path to the binarized Matrix Market file.
    :return: Validated binary matrix in CSR format.
    """
    with gzip.open(matrix_path, mode="rb") as handle:
        matrix = mmread(handle)

    if not sparse.issparse(matrix):
        raise TypeError(
            f"Expected a sparse matrix in {matrix_path}, "
            f"but found {type(matrix).__name__}."
        )

    matrix = matrix.tocsr()
    raw_nnz = matrix.nnz

    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()

    if matrix.nnz != raw_nnz:
        raise ValueError(
            f"{matrix_path} contained duplicate coordinates or explicit zeros. "
            f"The stored-entry count changed from {raw_nnz:,} to "
            f"{matrix.nnz:,} during cleanup."
        )

    if matrix.nnz > 0:
        if not np.isfinite(matrix.data).all():
            raise ValueError(
                f"{matrix_path} contains non-finite values."
            )

        if not np.all(matrix.data == 1):
            unique_values = np.unique(matrix.data)
            raise ValueError(
                f"{matrix_path} is not binary. Stored values are "
                f"{unique_values.tolist()}."
            )

    return matrix


def calculate_peak_support(
    matrix: sparse.csr_matrix,
    matrix_path: Path,
) -> np.ndarray:
    """
    Calculate the number of accessible cells for every peak row.

    :param matrix: Binarized peak-by-cell matrix.
    :param matrix_path: Matrix path used in validation error messages.
    :return: One integer support count per peak row.
    """
    support = np.asarray(matrix.sum(axis=1)).ravel().astype(
        np.int64,
        copy=False,
    )

    if support.shape[0] != matrix.shape[0]:
        raise RuntimeError(
            f"{matrix_path}: calculated {support.shape[0]:,} support values "
            f"for {matrix.shape[0]:,} matrix rows."
        )

    if support.size > 0:
        if support.min() < 0 or support.max() > matrix.shape[1]:
            raise RuntimeError(
                f"{matrix_path}: peak support is outside the valid range "
                f"[0, {matrix.shape[1]}]."
            )

    if int(support.sum()) != int(matrix.nnz):
        raise RuntimeError(
            f"{matrix_path}: support sum {int(support.sum()):,} does not "
            f"equal matrix nnz {matrix.nnz:,}."
        )

    return support


def count_peaks_at_least(
    support_histogram: np.ndarray,
    minimum_support: int,
) -> int:
    """
    Count peaks whose support is at least a requested threshold.

    :param support_histogram: Histogram in which index ``i`` stores the number
        of peaks with support ``i``.
    :param minimum_support: Minimum support count required.
    :return: Number of peaks meeting the threshold.
    """
    if minimum_support <= 0:
        return int(support_histogram.sum())

    if minimum_support >= support_histogram.size:
        return 0

    return int(support_histogram[minimum_support:].sum())


def support_category_counts(
    support_histogram: np.ndarray,
) -> dict[str, int]:
    """
    Summarize peak counts into interpretable low-support categories.

    :param support_histogram: Histogram of peak support counts.
    :return: Mapping from support category to number of peaks.
    """
    def get_count(index: int) -> int:
        """
        Return a histogram count while handling absent indices.

        :param index: Histogram index.
        :return: Count at the requested index, or zero when absent.
        """
        return (
            int(support_histogram[index])
            if index < support_histogram.size
            else 0
        )

    return {
        "0": get_count(0),
        "1": get_count(1),
        "2": get_count(2),
        "3-4": get_count(3) + get_count(4),
        ">=5": (
            int(support_histogram[5:].sum())
            if support_histogram.size > 5
            else 0
        ),
    }


def calculate_support_summary(
    sample_name: str,
    condition: str,
    cell_type: str,
    matrix_path: Path,
    matrix: sparse.csr_matrix,
    support: np.ndarray,
    support_histogram: np.ndarray,
) -> dict[str, Any]:
    """
    Calculate descriptive peak-support statistics for one matrix.

    :param sample_name: Sample-directory name.
    :param condition: Biological condition.
    :param cell_type: Standardized cell type.
    :param matrix_path: Path to the analyzed matrix.
    :param matrix: Binarized peak-by-cell matrix.
    :param support: Per-peak accessible-cell counts.
    :param support_histogram: Histogram of support counts.
    :return: One per-matrix summary record.
    """
    n_peaks, n_cells = matrix.shape
    categories = support_category_counts(support_histogram)

    if n_peaks == 0:
        raise ValueError(
            f"Matrix has no peak rows: {matrix_path}"
        )

    support_quantiles = np.quantile(
        support,
        [0.25, 0.50, 0.75, 0.90, 0.95, 0.99],
    )

    return {
        "sample_directory": sample_name,
        "condition": condition,
        "cell_type": cell_type,
        "matrix_rows_peaks": n_peaks,
        "matrix_columns_cells": n_cells,
        "nonzero_entries": matrix.nnz,
        "mean_peak_support_cells": float(support.mean()),
        "support_cells_q25": float(support_quantiles[0]),
        "support_cells_median": float(support_quantiles[1]),
        "support_cells_q75": float(support_quantiles[2]),
        "support_cells_q90": float(support_quantiles[3]),
        "support_cells_q95": float(support_quantiles[4]),
        "support_cells_q99": float(support_quantiles[5]),
        "maximum_peak_support_cells": int(support.max()),
        "mean_peak_support_fraction": float(support.mean() / n_cells),
        "support_fraction_q25": float(support_quantiles[0] / n_cells),
        "support_fraction_median": float(support_quantiles[1] / n_cells),
        "support_fraction_q75": float(support_quantiles[2] / n_cells),
        "support_fraction_q90": float(support_quantiles[3] / n_cells),
        "support_fraction_q95": float(support_quantiles[4] / n_cells),
        "support_fraction_q99": float(support_quantiles[5] / n_cells),
        "maximum_peak_support_fraction": float(support.max() / n_cells),
        "n_peaks_support_0": categories["0"],
        "n_peaks_support_1": categories["1"],
        "n_peaks_support_2": categories["2"],
        "n_peaks_support_3_to_4": categories["3-4"],
        "n_peaks_support_at_least_5": categories[">=5"],
        "fraction_peaks_support_0": categories["0"] / n_peaks,
        "fraction_peaks_support_1": categories["1"] / n_peaks,
        "fraction_peaks_support_2": categories["2"] / n_peaks,
        "fraction_peaks_support_3_to_4": categories["3-4"] / n_peaks,
        "fraction_peaks_support_at_least_5": categories[">=5"] / n_peaks,
        "n_unique_support_counts": int(np.count_nonzero(support_histogram)),
        "matrix_path": str(matrix_path),
    }


def build_support_distribution_records(
    sample_name: str,
    condition: str,
    cell_type: str,
    n_cells: int,
    support_histogram: np.ndarray,
) -> list[dict[str, Any]]:
    """
    Build compact support-distribution records for one matrix.

    :param sample_name: Sample-directory name.
    :param condition: Biological condition.
    :param cell_type: Standardized cell type.
    :param n_cells: Number of cells in the matrix.
    :param support_histogram: Histogram of peak support counts.
    :return: Records for support values that occur in at least one peak.
    """
    n_peaks = int(support_histogram.sum())
    records: list[dict[str, Any]] = []

    for support_count in np.flatnonzero(support_histogram):
        peak_count = int(support_histogram[support_count])

        records.append(
            {
                "sample_directory": sample_name,
                "condition": condition,
                "cell_type": cell_type,
                "n_cells": n_cells,
                "support_cells": int(support_count),
                "support_fraction": float(support_count / n_cells),
                "n_peaks": peak_count,
                "fraction_of_matrix_peaks": peak_count / n_peaks,
            }
        )

    return records


def build_hybrid_threshold_records(
    sample_name: str,
    condition: str,
    cell_type: str,
    n_cells: int,
    n_peaks: int,
    support_histogram: np.ndarray,
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
    minimum_cell_fractions: Sequence[float],
) -> list[dict[str, Any]]:
    """
    Evaluate every candidate hybrid threshold for one matrix.

    :param sample_name: Sample-directory name.
    :param condition: Biological condition.
    :param cell_type: Standardized cell type.
    :param n_cells: Number of matrix columns.
    :param n_peaks: Number of matrix rows.
    :param support_histogram: Histogram of peak support counts.
    :param minimum_group_sizes: Candidate minimum usable group sizes.
    :param minimum_cell_supports: Candidate absolute support thresholds.
    :param minimum_cell_fractions: Candidate fractional support thresholds.
    :return: One record for every candidate hyperparameter combination.
    """
    records: list[dict[str, Any]] = []

    for minimum_group_size in minimum_group_sizes:
        eligible = n_cells >= minimum_group_size

        for minimum_cell_support in minimum_cell_supports:
            for minimum_cell_fraction in minimum_cell_fractions:
                fractional_support = math.ceil(
                    minimum_cell_fraction * n_cells
                )
                effective_support = max(
                    minimum_cell_support,
                    fractional_support,
                )

                if eligible:
                    retained_peaks: int | None = count_peaks_at_least(
                        support_histogram=support_histogram,
                        minimum_support=effective_support,
                    )
                    retained_fraction: float | None = (
                        retained_peaks / n_peaks
                    )
                else:
                    retained_peaks = None
                    retained_fraction = None

                records.append(
                    {
                        "sample_directory": sample_name,
                        "condition": condition,
                        "cell_type": cell_type,
                        "n_cells": n_cells,
                        "n_input_peaks": n_peaks,
                        "min_group_size": minimum_group_size,
                        "min_cell_support": minimum_cell_support,
                        "min_cell_fraction": minimum_cell_fraction,
                        "fractional_support_ceiling": fractional_support,
                        "effective_min_support_cells": effective_support,
                        "group_eligible": eligible,
                        "n_retained_peaks": retained_peaks,
                        "fraction_retained_peaks": retained_fraction,
                    }
                )

    return records


def analyze_matrices(
    indexed_matrices: dict[tuple[str, str], Path],
    sample_order: Sequence[str],
    cell_type_order: Sequence[str],
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
    minimum_cell_fractions: Sequence[float],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Analyze all expected sample-by-cell-type groups.

    :param indexed_matrices: Mapping from sample and cell type to matrix path.
    :param sample_order: Expected sample order.
    :param cell_type_order: Expected cell-type order.
    :param minimum_group_sizes: Candidate minimum usable group sizes.
    :param minimum_cell_supports: Candidate absolute support thresholds.
    :param minimum_cell_fractions: Candidate fractional support thresholds.
    :return: Group-size, per-matrix support, support-distribution, and hybrid
        threshold tables.
    """
    group_size_records: list[dict[str, Any]] = []
    support_summary_records: list[dict[str, Any]] = []
    support_distribution_records: list[dict[str, Any]] = []
    hybrid_records: list[dict[str, Any]] = []

    for sample_name in sample_order:
        condition = infer_condition(sample_name)

        for cell_type in cell_type_order:
            matrix_path = indexed_matrices.get((sample_name, cell_type))

            if matrix_path is None:
                group_size_records.append(
                    {
                        "sample_directory": sample_name,
                        "condition": condition,
                        "cell_type": cell_type,
                        "matrix_available": False,
                        "n_cells": 0,
                        "n_peaks": None,
                        "nonzero_entries": None,
                        "matrix_path": None,
                    }
                )
                continue

            print(f"Processing {sample_name} / {cell_type}")

            matrix = read_binary_matrix(matrix_path)
            support = calculate_peak_support(
                matrix=matrix,
                matrix_path=matrix_path,
            )
            support_histogram = np.bincount(
                support,
                minlength=matrix.shape[1] + 1,
            ).astype(np.int64, copy=False)

            if int(support_histogram.sum()) != matrix.shape[0]:
                raise RuntimeError(
                    f"{matrix_path}: support histogram contains "
                    f"{int(support_histogram.sum()):,} peaks, expected "
                    f"{matrix.shape[0]:,}."
                )

            group_size_records.append(
                {
                    "sample_directory": sample_name,
                    "condition": condition,
                    "cell_type": cell_type,
                    "matrix_available": True,
                    "n_cells": matrix.shape[1],
                    "n_peaks": matrix.shape[0],
                    "nonzero_entries": matrix.nnz,
                    "matrix_path": str(matrix_path),
                }
            )

            support_summary_records.append(
                calculate_support_summary(
                    sample_name=sample_name,
                    condition=condition,
                    cell_type=cell_type,
                    matrix_path=matrix_path,
                    matrix=matrix,
                    support=support,
                    support_histogram=support_histogram,
                )
            )

            support_distribution_records.extend(
                build_support_distribution_records(
                    sample_name=sample_name,
                    condition=condition,
                    cell_type=cell_type,
                    n_cells=matrix.shape[1],
                    support_histogram=support_histogram,
                )
            )

            hybrid_records.extend(
                build_hybrid_threshold_records(
                    sample_name=sample_name,
                    condition=condition,
                    cell_type=cell_type,
                    n_cells=matrix.shape[1],
                    n_peaks=matrix.shape[0],
                    support_histogram=support_histogram,
                    minimum_group_sizes=minimum_group_sizes,
                    minimum_cell_supports=minimum_cell_supports,
                    minimum_cell_fractions=minimum_cell_fractions,
                )
            )

    return (
        pd.DataFrame.from_records(group_size_records),
        pd.DataFrame.from_records(support_summary_records),
        pd.DataFrame.from_records(support_distribution_records),
        pd.DataFrame.from_records(hybrid_records),
    )


def build_group_size_sensitivity(
    group_sizes: pd.DataFrame,
    minimum_group_sizes: Sequence[int],
    cell_type_order: Sequence[str],
) -> pd.DataFrame:
    """
    Summarize how candidate group-size thresholds affect matrix availability.

    :param group_sizes: Expected sample-by-cell-type group table.
    :param minimum_group_sizes: Candidate minimum usable group sizes.
    :param cell_type_order: Ordered standardized cell types.
    :return: Group-size threshold sensitivity table.
    """
    records: list[dict[str, Any]] = []

    for minimum_group_size in minimum_group_sizes:
        for cell_type_label in ("ALL", *cell_type_order):
            if cell_type_label == "ALL":
                subset = group_sizes
            else:
                subset = group_sizes.loc[
                    group_sizes["cell_type"] == cell_type_label
                ]

            available = subset["matrix_available"].astype(bool)
            eligible = available & (
                subset["n_cells"] >= minimum_group_size
            )

            record: dict[str, Any] = {
                "min_group_size": minimum_group_size,
                "cell_type": cell_type_label,
                "n_expected_groups": len(subset),
                "n_available_groups": int(available.sum()),
                "n_missing_groups": int((~available).sum()),
                "n_eligible_groups": int(eligible.sum()),
                "n_available_but_too_small": int(
                    (available & ~eligible).sum()
                ),
                "fraction_expected_groups_eligible": float(
                    eligible.mean()
                ),
                "fraction_available_groups_eligible": float(
                    eligible.sum() / available.sum()
                    if available.sum() > 0
                    else 0.0
                ),
            }

            for condition in CONDITION_ORDER:
                condition_mask = subset["condition"] == condition
                condition_available = available & condition_mask
                condition_eligible = eligible & condition_mask

                record[f"n_available_{condition.lower()}"] = int(
                    condition_available.sum()
                )
                record[f"n_eligible_{condition.lower()}"] = int(
                    condition_eligible.sum()
                )

            records.append(record)

    return pd.DataFrame.from_records(records)


def build_hyperparameter_grid(
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
    minimum_cell_fractions: Sequence[float],
) -> pd.DataFrame:
    """
    Build a table containing every evaluated hyperparameter combination.

    :param minimum_group_sizes: Candidate minimum usable group sizes.
    :param minimum_cell_supports: Candidate absolute support thresholds.
    :param minimum_cell_fractions: Candidate fractional support thresholds.
    :return: Complete hyperparameter grid.
    """
    records: list[dict[str, Any]] = []

    for minimum_group_size in minimum_group_sizes:
        for minimum_cell_support in minimum_cell_supports:
            for minimum_cell_fraction in minimum_cell_fractions:
                records.append(
                    {
                        "min_group_size": minimum_group_size,
                        "min_cell_support": minimum_cell_support,
                        "min_cell_fraction": minimum_cell_fraction,
                        "hybrid_rule": (
                            "support >= max("
                            f"{minimum_cell_support}, "
                            f"ceil({minimum_cell_fraction} * n_cells))"
                        ),
                    }
                )

    return pd.DataFrame.from_records(records)


def summarize_hybrid_thresholds(
    hybrid_per_matrix: pd.DataFrame,
    grouping_columns: Sequence[str],
) -> pd.DataFrame:
    """
    Aggregate per-matrix hybrid-threshold results.

    :param hybrid_per_matrix: Per-matrix threshold evaluation table.
    :param grouping_columns: Additional columns defining each summary group.
    :return: Aggregated threshold summary table.
    """
    eligible = hybrid_per_matrix.loc[
        hybrid_per_matrix["group_eligible"]
    ].copy()

    base_grouping_columns = [
        "min_group_size",
        "min_cell_support",
        "min_cell_fraction",
    ]
    all_grouping_columns = [
        *grouping_columns,
        *base_grouping_columns,
    ]

    records: list[dict[str, Any]] = []

    for group_key, group in eligible.groupby(
        all_grouping_columns,
        sort=True,
        dropna=False,
    ):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)

        key_values = dict(zip(all_grouping_columns, group_key))
        retained_fractions = group[
            "fraction_retained_peaks"
        ].astype(float)

        total_input_peaks = int(
            group["n_input_peaks"].sum()
        )
        total_retained_peaks = int(
            group["n_retained_peaks"].sum()
        )

        records.append(
            {
                **key_values,
                "n_eligible_matrices": len(group),
                "n_mash_matrices": int(
                    (group["condition"] == "MASH").sum()
                ),
                "n_normal_matrices": int(
                    (group["condition"] == "Normal").sum()
                ),
                "total_input_peaks": total_input_peaks,
                "total_retained_peaks": total_retained_peaks,
                "pooled_fraction_retained": (
                    total_retained_peaks / total_input_peaks
                ),
                "mean_fraction_retained": float(
                    retained_fractions.mean()
                ),
                "median_fraction_retained": float(
                    retained_fractions.median()
                ),
                "fraction_retained_q25": float(
                    retained_fractions.quantile(0.25)
                ),
                "fraction_retained_q75": float(
                    retained_fractions.quantile(0.75)
                ),
                "minimum_fraction_retained": float(
                    retained_fractions.min()
                ),
                "maximum_fraction_retained": float(
                    retained_fractions.max()
                ),
                "median_effective_min_support_cells": float(
                    group["effective_min_support_cells"].median()
                ),
                "minimum_effective_min_support_cells": int(
                    group["effective_min_support_cells"].min()
                ),
                "maximum_effective_min_support_cells": int(
                    group["effective_min_support_cells"].max()
                ),
            }
        )

    return pd.DataFrame.from_records(records)


def get_plot_colors(n_colors: int) -> list[PlotColor]:
    """
    Return a plot palette with clearly distinguishable colors.

    A single-series plot uses ``PLOTS_MAIN_COLOR``. A two-series plot uses
    ``PLOTS_MAIN_COLOR`` and ``PLOTS_SECOND_COLOR``. Plots with more than two
    series use a qualitative Matplotlib palette so adjacent series remain easy
    to distinguish.

    :param n_colors: Number of colors required by the plot.
    :return: Ordered list of Matplotlib-compatible colors.
    """
    if n_colors < 1:
        raise ValueError("n_colors must be at least 1.")

    if n_colors == 1:
        return [PLOTS_MAIN_COLOR]

    if n_colors == 2:
        return [
            PLOTS_MAIN_COLOR,
            PLOTS_SECOND_COLOR,
        ]

    colormap_name = (
        "tab10"
        if n_colors <= 10
        else "tab20"
    )
    colormap = matplotlib.colormaps[
        colormap_name
    ].resampled(n_colors)

    return [
        colormap(index)
        for index in range(n_colors)
    ]


def get_heatmap_text_color(
    value: float,
    minimum_value: float,
    maximum_value: float,
) -> str:
    """
    Choose readable annotation text for a continuous heatmap cell.

    :param value: Numeric value represented by the heatmap cell.
    :param minimum_value: Lower bound of the heatmap color scale.
    :param maximum_value: Upper bound of the heatmap color scale.
    :return: Either white or black, depending on the normalized cell value.
    """
    if not np.isfinite(value):
        return "black"

    if maximum_value <= minimum_value:
        normalized_value = 0.0
    else:
        normalized_value = (
            value - minimum_value
        ) / (
            maximum_value - minimum_value
        )

    return (
        "white"
        if normalized_value < 0.55
        else "black"
    )



def save_figure(
    figure: plt.Figure,
    output_path: Path,
) -> None:
    """
    Save and close one Matplotlib figure.

    :param figure: Figure to save.
    :param output_path: Destination PNG path.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        output_path,
        dpi=PLOT_DPI,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_group_size_heatmap(
    group_sizes: pd.DataFrame,
    sample_order: Sequence[str],
    cell_type_order: Sequence[str],
) -> None:
    """
    Plot the number of cells in every expected sample-by-cell-type group.

    :param group_sizes: Expected group-size table.
    :param sample_order: Ordered sample names.
    :param cell_type_order: Ordered cell types.
    """
    pivot = group_sizes.pivot(
        index="sample_directory",
        columns="cell_type",
        values="n_cells",
    ).reindex(
        index=sample_order,
        columns=cell_type_order,
    )

    values = pivot.to_numpy(dtype=float)
    transformed_values = np.log10(values + 1.0)

    figure, axis = plt.subplots(
        figsize=(13, 8),
    )

    image = axis.imshow(
        transformed_values,
        aspect="auto",
        cmap=HEATMAP_COLORMAP,
    )

    axis.set_xticks(
        np.arange(len(cell_type_order)),
        labels=cell_type_order,
        rotation=35,
        ha="right",
    )
    axis.set_yticks(
        np.arange(len(sample_order)),
        labels=sample_order,
    )
    axis.set_title(
        "Cells per sample and cell type"
    )
    axis.set_xlabel("Cell type")
    axis.set_ylabel("Sample")

    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = int(values[row_index, column_index])
            label = f"{value:,}" if value > 0 else "—"
            axis.text(
                column_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=7,
                color=get_heatmap_text_color(
                    value=transformed_values[
                        row_index,
                        column_index,
                    ],
                    minimum_value=float(
                        np.nanmin(transformed_values)
                    ),
                    maximum_value=float(
                        np.nanmax(transformed_values)
                    ),
                ),
            )

    colorbar = figure.colorbar(
        image,
        ax=axis,
        fraction=0.025,
        pad=0.02,
    )
    colorbar.set_label("log10(number of cells + 1)")

    figure.tight_layout()

    save_figure(
        figure=figure,
        output_path=PLOTS_DIR / "group_size_heatmap.png",
    )


def plot_usable_groups_overall(
    group_size_sensitivity: pd.DataFrame,
) -> None:
    """
    Plot the number of eligible groups across minimum group-size values.

    :param group_size_sensitivity: Group-size threshold sensitivity table.
    """
    overall = group_size_sensitivity.loc[
        group_size_sensitivity["cell_type"] == "ALL"
    ].sort_values("min_group_size")

    figure, axis = plt.subplots(
        figsize=(9, 5.5),
    )

    colors = get_plot_colors(3)

    series = (
        ("Total", "n_eligible_groups", colors[0]),
        ("MASH", "n_eligible_mash", colors[1]),
        ("Normal", "n_eligible_normal", colors[2]),
    )

    for label, column, color in series:
        axis.plot(
            overall["min_group_size"],
            overall[column],
            marker="o",
            linewidth=2,
            label=label,
            color=color,
        )

    axis.set_title(
        "Usable sample-by-cell-type groups versus minimum cell count"
    )
    axis.set_xlabel("Minimum cells required in a group")
    axis.set_ylabel("Number of eligible groups")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()

    save_figure(
        figure=figure,
        output_path=PLOTS_DIR / "usable_groups_overall.png",
    )


def plot_usable_groups_by_cell_type(
    group_size_sensitivity: pd.DataFrame,
    cell_type_order: Sequence[str],
) -> None:
    """
    Plot eligible-group counts by cell type across group-size thresholds.

    :param group_size_sensitivity: Group-size threshold sensitivity table.
    :param cell_type_order: Ordered standardized cell types.
    """
    colors = get_plot_colors(
        len(cell_type_order)
    )

    figure, axis = plt.subplots(
        figsize=(10, 6),
    )

    for cell_type, color in zip(cell_type_order, colors):
        subset = group_size_sensitivity.loc[
            group_size_sensitivity["cell_type"] == cell_type
        ].sort_values("min_group_size")

        axis.plot(
            subset["min_group_size"],
            subset["n_eligible_groups"],
            marker="o",
            linewidth=1.8,
            label=cell_type,
            color=color,
        )

    axis.set_title(
        "Usable biological samples by cell type"
    )
    axis.set_xlabel("Minimum cells required in a group")
    axis.set_ylabel("Number of eligible samples")
    axis.set_ylim(bottom=0)
    axis.grid(axis="y", alpha=0.25)
    axis.legend(
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
    )
    figure.tight_layout()

    save_figure(
        figure=figure,
        output_path=PLOTS_DIR / "usable_groups_by_cell_type.png",
    )


def plot_support_categories_by_cell_type(
    per_matrix_support: pd.DataFrame,
    cell_type_order: Sequence[str],
) -> None:
    """
    Plot pooled low-support peak categories for every cell type.

    :param per_matrix_support: Per-matrix support summary table.
    :param cell_type_order: Ordered standardized cell types.
    """
    category_columns = {
        "0": "n_peaks_support_0",
        "1": "n_peaks_support_1",
        "2": "n_peaks_support_2",
        "3-4": "n_peaks_support_3_to_4",
        ">=5": "n_peaks_support_at_least_5",
    }

    aggregated = (
        per_matrix_support.groupby("cell_type", sort=False)[
            list(category_columns.values())
        ]
        .sum()
        .reindex(cell_type_order)
    )

    proportions = aggregated.div(
        aggregated.sum(axis=1),
        axis=0,
    )

    colors = get_plot_colors(
        len(SUPPORT_CATEGORY_ORDER)
    )

    figure, axis = plt.subplots(
        figsize=(11, 6),
    )

    bottom = np.zeros(len(cell_type_order), dtype=float)

    for category, color in zip(SUPPORT_CATEGORY_ORDER, colors):
        values = proportions[
            category_columns[category]
        ].to_numpy(dtype=float)

        axis.bar(
            cell_type_order,
            values,
            bottom=bottom,
            label=category,
            color=color,
        )
        bottom += values

    axis.set_title(
        "Peak-support composition by cell type"
    )
    axis.set_xlabel("Cell type")
    axis.set_ylabel("Fraction of sample-specific peak rows")
    axis.set_ylim(0.0, 1.0)
    axis.tick_params(
        axis="x",
        rotation=35,
    )
    axis.legend(
        title="Accessible cells",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
    )
    axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()

    save_figure(
        figure=figure,
        output_path=PLOTS_DIR / "support_categories_by_cell_type.png",
    )


def format_fraction_label(value: float) -> str:
    """
    Format a fractional hyperparameter as a concise percentage label.

    :param value: Fractional value in the interval [0, 1].
    :return: Human-readable percentage label.
    """
    percentage = value * 100.0

    if percentage.is_integer():
        return f"{percentage:.0f}%"

    return f"{percentage:g}%"


def plot_overall_retention_heatmaps(
    overall_summary: pd.DataFrame,
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
    minimum_cell_fractions: Sequence[float],
) -> None:
    """
    Plot overall pooled retained-peak fractions for every threshold grid.

    :param overall_summary: Overall hybrid-threshold summary.
    :param minimum_group_sizes: Candidate minimum usable group sizes.
    :param minimum_cell_supports: Candidate absolute support thresholds.
    :param minimum_cell_fractions: Candidate fractional support thresholds.
    """
    for minimum_group_size in minimum_group_sizes:
        subset = overall_summary.loc[
            overall_summary["min_group_size"] == minimum_group_size
        ]

        pivot = subset.pivot(
            index="min_cell_support",
            columns="min_cell_fraction",
            values="pooled_fraction_retained",
        ).reindex(
            index=minimum_cell_supports,
            columns=minimum_cell_fractions,
        )

        values = pivot.to_numpy(dtype=float)

        figure, axis = plt.subplots(
            figsize=(8.5, 6),
        )

        image = axis.imshow(
            values,
            vmin=0.0,
            vmax=1.0,
            aspect="auto",
            cmap=HEATMAP_COLORMAP,
        )

        axis.set_xticks(
            np.arange(len(minimum_cell_fractions)),
            labels=[
                format_fraction_label(value)
                for value in minimum_cell_fractions
            ],
        )
        axis.set_yticks(
            np.arange(len(minimum_cell_supports)),
            labels=minimum_cell_supports,
        )
        axis.set_title(
            "Pooled retained peak fraction\n"
            f"Minimum group size: {minimum_group_size} cells"
        )
        axis.set_xlabel("Minimum accessible-cell fraction")
        axis.set_ylabel("Minimum accessible cells")

        for row_index in range(values.shape[0]):
            for column_index in range(values.shape[1]):
                value = values[row_index, column_index]

                if np.isnan(value):
                    label = "NA"
                else:
                    label = f"{value:.1%}"

                axis.text(
                    column_index,
                    row_index,
                    label,
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=get_heatmap_text_color(
                        value=value,
                        minimum_value=0.0,
                        maximum_value=1.0,
                    ),
                )

        colorbar = figure.colorbar(
            image,
            ax=axis,
            fraction=0.045,
            pad=0.04,
        )
        colorbar.set_label("Pooled fraction of peaks retained")
        figure.tight_layout()

        save_figure(
            figure=figure,
            output_path=(
                OVERALL_HEATMAPS_DIR
                / f"min_group_size_{minimum_group_size}.png"
            ),
        )


def plot_cell_type_retention_heatmaps(
    cell_type_summary: pd.DataFrame,
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
    minimum_cell_fractions: Sequence[float],
    cell_type_order: Sequence[str],
) -> None:
    """
    Plot retained-peak fractions by cell type for every group-size threshold.

    :param cell_type_summary: Cell-type-specific hybrid-threshold summary.
    :param minimum_group_sizes: Candidate minimum usable group sizes.
    :param minimum_cell_supports: Candidate absolute support thresholds.
    :param minimum_cell_fractions: Candidate fractional support thresholds.
    :param cell_type_order: Ordered standardized cell types.
    """
    threshold_pairs = [
        (minimum_support, minimum_fraction)
        for minimum_support in minimum_cell_supports
        for minimum_fraction in minimum_cell_fractions
    ]
    threshold_labels = [
        (
            f"k≥{minimum_support}\n"
            f"q≥{format_fraction_label(minimum_fraction)}"
        )
        for minimum_support, minimum_fraction in threshold_pairs
    ]

    for minimum_group_size in minimum_group_sizes:
        subset = cell_type_summary.loc[
            cell_type_summary["min_group_size"] == minimum_group_size
        ].copy()

        lookup = {
            (
                row.cell_type,
                int(row.min_cell_support),
                float(row.min_cell_fraction),
            ): float(row.pooled_fraction_retained)
            for row in subset.itertuples(index=False)
        }

        values = np.full(
            (
                len(cell_type_order),
                len(threshold_pairs),
            ),
            np.nan,
            dtype=float,
        )

        for row_index, cell_type in enumerate(cell_type_order):
            for column_index, (
                minimum_support,
                minimum_fraction,
            ) in enumerate(threshold_pairs):
                values[row_index, column_index] = lookup.get(
                    (
                        cell_type,
                        minimum_support,
                        minimum_fraction,
                    ),
                    np.nan,
                )

        figure, axis = plt.subplots(
            figsize=(22, 6.5),
        )

        image = axis.imshow(
            values,
            vmin=0.0,
            vmax=1.0,
            aspect="auto",
            cmap=HEATMAP_COLORMAP,
        )

        axis.set_xticks(
            np.arange(len(threshold_labels)),
            labels=threshold_labels,
            rotation=55,
            ha="right",
            fontsize=8,
        )
        axis.set_yticks(
            np.arange(len(cell_type_order)),
            labels=cell_type_order,
        )
        axis.set_title(
            "Pooled retained peak fraction by cell type\n"
            f"Minimum group size: {minimum_group_size} cells"
        )
        axis.set_xlabel("Hybrid threshold")
        axis.set_ylabel("Cell type")

        colorbar = figure.colorbar(
            image,
            ax=axis,
            fraction=0.018,
            pad=0.02,
        )
        colorbar.set_label("Pooled fraction of peaks retained")
        figure.tight_layout()

        save_figure(
            figure=figure,
            output_path=(
                CELL_TYPE_HEATMAPS_DIR
                / f"min_group_size_{minimum_group_size}.png"
            ),
        )


def plot_threshold_curves(
    overall_summary: pd.DataFrame,
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
) -> None:
    """
    Plot pooled retained-peak fractions as fractional thresholds increase.

    :param overall_summary: Overall hybrid-threshold summary.
    :param minimum_group_sizes: Candidate minimum usable group sizes.
    :param minimum_cell_supports: Candidate absolute support thresholds.
    """
    colors = get_plot_colors(
        len(minimum_cell_supports)
    )

    for minimum_group_size in minimum_group_sizes:
        subset = overall_summary.loc[
            overall_summary["min_group_size"] == minimum_group_size
        ]

        figure, axis = plt.subplots(
            figsize=(9, 5.5),
        )

        for minimum_support, color in zip(
            minimum_cell_supports,
            colors,
        ):
            support_subset = subset.loc[
                subset["min_cell_support"] == minimum_support
            ].sort_values("min_cell_fraction")

            axis.plot(
                support_subset["min_cell_fraction"] * 100.0,
                support_subset["pooled_fraction_retained"],
                marker="o",
                linewidth=2,
                label=f"k ≥ {minimum_support}",
                color=color,
            )

        axis.set_title(
            "Hybrid peak-filtering sensitivity\n"
            f"Minimum group size: {minimum_group_size} cells"
        )
        axis.set_xlabel("Minimum accessible-cell fraction (%)")
        axis.set_ylabel("Pooled fraction of peaks retained")
        axis.set_ylim(0.0, 1.0)
        axis.grid(alpha=0.25)
        axis.legend(
            title="Absolute support",
            bbox_to_anchor=(1.02, 1.0),
            loc="upper left",
        )
        figure.tight_layout()

        save_figure(
            figure=figure,
            output_path=(
                THRESHOLD_CURVES_DIR
                / f"min_group_size_{minimum_group_size}.png"
            ),
        )


def write_manifest(
    minimum_group_sizes: Sequence[int],
    minimum_cell_supports: Sequence[int],
    minimum_cell_fractions: Sequence[float],
) -> None:
    """
    Write a human-readable description of outputs and definitions.

    :param minimum_group_sizes: Evaluated minimum usable group sizes.
    :param minimum_cell_supports: Evaluated absolute support thresholds.
    :param minimum_cell_fractions: Evaluated fractional support thresholds.
    """
    manifest = f"""Peak-filtering hyperparameter analysis
========================================

Script version
--------------
{SCRIPT_VERSION}

This directory contains sensitivity analyses only. No peak or matrix was
filtered or overwritten.

Definitions
-----------
For one sample-by-cell-type binary matrix:

    support(p) = number of cells with value 1 for peak p

A group is eligible when:

    n_cells >= min_group_size

A peak passes the hybrid rule when:

    support(p) >= max(
        min_cell_support,
        ceil(min_cell_fraction * n_cells)
    )

Evaluated values
----------------
min_group_size:
    {list(minimum_group_sizes)}

min_cell_support:
    {list(minimum_cell_supports)}

min_cell_fraction:
    {list(minimum_cell_fractions)}

Tables
------
group_sizes.csv
    Availability and dimensions of every expected sample-by-cell-type group.

minimum_group_size_sensitivity.csv
    Number of eligible groups under each minimum group-size value.

per_matrix_peak_support_summary.csv
    Descriptive support statistics for every available matrix.

peak_support_count_distributions.csv.gz
    Compact support-frequency table for every matrix.

hyperparameter_grid.csv
    Every candidate hyperparameter combination and its rule.

hybrid_threshold_per_matrix.csv.gz
    Retained peak counts and fractions for every matrix and every candidate
    combination.

hybrid_threshold_summary_overall.csv
    Aggregate sensitivity across all eligible matrices.

hybrid_threshold_summary_by_condition.csv
    Aggregate sensitivity separately for MASH and Normal.

hybrid_threshold_summary_by_cell_type.csv
    Aggregate sensitivity separately for each cell type.

hybrid_threshold_summary_by_cell_type_and_condition.csv
    Aggregate sensitivity for every cell-type-by-condition combination.

Plots
-----
plots/group_size_heatmap.png
    Cell counts for all expected sample-by-cell-type groups.

plots/usable_groups_overall.png
    Total, MASH, and Normal eligible groups versus min_group_size.

plots/usable_groups_by_cell_type.png
    Eligible biological samples per cell type versus min_group_size.

plots/support_categories_by_cell_type.png
    Pooled fractions of peaks with support 0, 1, 2, 3-4, or >=5 cells.

plots/retained_fraction_overall/
    Overall k-by-q retained-fraction heatmaps for every min_group_size.

plots/retained_fraction_by_cell_type/
    Cell-type-specific retained-fraction heatmaps for every min_group_size.

plots/threshold_curves/
    Retained-fraction curves across q values, with one line for each k.

Important limitation
--------------------
This analysis evaluates local sample-specific peaks. The replicate-support
hyperparameter cannot yet be calibrated because rows are not aligned across
samples. It should be analyzed after projection onto the common consensus-peak
set.
"""

    MANIFEST_PATH.write_text(
        manifest,
        encoding="utf-8",
    )


def write_tables(
    group_sizes: pd.DataFrame,
    group_size_sensitivity: pd.DataFrame,
    per_matrix_support: pd.DataFrame,
    support_distributions: pd.DataFrame,
    hyperparameter_grid: pd.DataFrame,
    hybrid_per_matrix: pd.DataFrame,
    hybrid_overall: pd.DataFrame,
    hybrid_by_condition: pd.DataFrame,
    hybrid_by_cell_type: pd.DataFrame,
    hybrid_by_cell_type_condition: pd.DataFrame,
) -> None:
    """
    Write all analysis tables to the configured results directory.

    :param group_sizes: Expected group-size table.
    :param group_size_sensitivity: Minimum group-size sensitivity table.
    :param per_matrix_support: Per-matrix support summary.
    :param support_distributions: Compact support-frequency table.
    :param hyperparameter_grid: Evaluated hyperparameter combinations.
    :param hybrid_per_matrix: Per-matrix threshold results.
    :param hybrid_overall: Overall threshold summary.
    :param hybrid_by_condition: Condition-specific threshold summary.
    :param hybrid_by_cell_type: Cell-type-specific threshold summary.
    :param hybrid_by_cell_type_condition: Cell-type-by-condition summary.
    """
    PEAK_FILTERING_RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    group_sizes.to_csv(
        GROUP_SIZES_PATH,
        index=False,
    )
    group_size_sensitivity.to_csv(
        GROUP_SIZE_SENSITIVITY_PATH,
        index=False,
    )
    per_matrix_support.to_csv(
        PER_MATRIX_SUPPORT_SUMMARY_PATH,
        index=False,
    )
    support_distributions.to_csv(
        SUPPORT_DISTRIBUTIONS_PATH,
        index=False,
        compression="gzip",
    )
    hyperparameter_grid.to_csv(
        HYPERPARAMETER_GRID_PATH,
        index=False,
    )
    hybrid_per_matrix.to_csv(
        HYBRID_PER_MATRIX_PATH,
        index=False,
        compression="gzip",
    )
    hybrid_overall.to_csv(
        HYBRID_OVERALL_SUMMARY_PATH,
        index=False,
    )
    hybrid_by_condition.to_csv(
        HYBRID_CONDITION_SUMMARY_PATH,
        index=False,
    )
    hybrid_by_cell_type.to_csv(
        HYBRID_CELL_TYPE_SUMMARY_PATH,
        index=False,
    )
    hybrid_by_cell_type_condition.to_csv(
        HYBRID_CELL_TYPE_CONDITION_SUMMARY_PATH,
        index=False,
    )


def main() -> None:
    """
    Run the complete peak-filtering hyperparameter sensitivity analysis.
    """
    (
        minimum_group_sizes,
        minimum_cell_supports,
        minimum_cell_fractions,
    ) = validate_hyperparameter_values(
        minimum_group_sizes=PEAK_FILTER_MIN_GROUP_SIZES,
        minimum_cell_supports=PEAK_FILTER_MIN_CELL_SUPPORTS,
        minimum_cell_fractions=PEAK_FILTER_MIN_CELL_FRACTIONS,
    )

    cell_type_order = get_cell_type_order()
    sample_order = tuple(
        sample_directory.name
        for sample_directory in ATAC_SEQ_DIRS
    )

    PEAK_FILTERING_RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    PLOTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    indexed_matrices = index_binary_matrices(
        cell_type_order=cell_type_order,
    )

    (
        group_sizes,
        per_matrix_support,
        support_distributions,
        hybrid_per_matrix,
    ) = analyze_matrices(
        indexed_matrices=indexed_matrices,
        sample_order=sample_order,
        cell_type_order=cell_type_order,
        minimum_group_sizes=minimum_group_sizes,
        minimum_cell_supports=minimum_cell_supports,
        minimum_cell_fractions=minimum_cell_fractions,
    )

    group_size_sensitivity = build_group_size_sensitivity(
        group_sizes=group_sizes,
        minimum_group_sizes=minimum_group_sizes,
        cell_type_order=cell_type_order,
    )

    hyperparameter_grid = build_hyperparameter_grid(
        minimum_group_sizes=minimum_group_sizes,
        minimum_cell_supports=minimum_cell_supports,
        minimum_cell_fractions=minimum_cell_fractions,
    )

    hybrid_overall = summarize_hybrid_thresholds(
        hybrid_per_matrix=hybrid_per_matrix,
        grouping_columns=(),
    )
    hybrid_by_condition = summarize_hybrid_thresholds(
        hybrid_per_matrix=hybrid_per_matrix,
        grouping_columns=("condition",),
    )
    hybrid_by_cell_type = summarize_hybrid_thresholds(
        hybrid_per_matrix=hybrid_per_matrix,
        grouping_columns=("cell_type",),
    )
    hybrid_by_cell_type_condition = summarize_hybrid_thresholds(
        hybrid_per_matrix=hybrid_per_matrix,
        grouping_columns=("cell_type", "condition"),
    )

    write_tables(
        group_sizes=group_sizes,
        group_size_sensitivity=group_size_sensitivity,
        per_matrix_support=per_matrix_support,
        support_distributions=support_distributions,
        hyperparameter_grid=hyperparameter_grid,
        hybrid_per_matrix=hybrid_per_matrix,
        hybrid_overall=hybrid_overall,
        hybrid_by_condition=hybrid_by_condition,
        hybrid_by_cell_type=hybrid_by_cell_type,
        hybrid_by_cell_type_condition=hybrid_by_cell_type_condition,
    )

    write_manifest(
        minimum_group_sizes=minimum_group_sizes,
        minimum_cell_supports=minimum_cell_supports,
        minimum_cell_fractions=minimum_cell_fractions,
    )

    plot_group_size_heatmap(
        group_sizes=group_sizes,
        sample_order=sample_order,
        cell_type_order=cell_type_order,
    )
    plot_usable_groups_overall(
        group_size_sensitivity=group_size_sensitivity,
    )
    plot_usable_groups_by_cell_type(
        group_size_sensitivity=group_size_sensitivity,
        cell_type_order=cell_type_order,
    )
    plot_support_categories_by_cell_type(
        per_matrix_support=per_matrix_support,
        cell_type_order=cell_type_order,
    )
    plot_overall_retention_heatmaps(
        overall_summary=hybrid_overall,
        minimum_group_sizes=minimum_group_sizes,
        minimum_cell_supports=minimum_cell_supports,
        minimum_cell_fractions=minimum_cell_fractions,
    )
    plot_cell_type_retention_heatmaps(
        cell_type_summary=hybrid_by_cell_type,
        minimum_group_sizes=minimum_group_sizes,
        minimum_cell_supports=minimum_cell_supports,
        minimum_cell_fractions=minimum_cell_fractions,
        cell_type_order=cell_type_order,
    )
    plot_threshold_curves(
        overall_summary=hybrid_overall,
        minimum_group_sizes=minimum_group_sizes,
        minimum_cell_supports=minimum_cell_supports,
    )

    print("\nPeak-filtering hyperparameter analysis completed.")
    print(
        f"Available matrices analyzed: {len(per_matrix_support):,}"
    )
    print(
        f"Expected sample-by-cell-type groups represented: "
        f"{len(group_sizes):,}"
    )
    print(
        f"Hyperparameter combinations evaluated: "
        f"{len(hyperparameter_grid):,}"
    )
    print(
        f"Results written to: {PEAK_FILTERING_RESULTS_DIR}"
    )


if __name__ == "__main__":
    main()
