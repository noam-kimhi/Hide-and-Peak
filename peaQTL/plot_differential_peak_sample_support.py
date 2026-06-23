#!/usr/bin/env python3
"""Visualize sample-level support for DESeq2 differential ATAC peaks.

The DESeq2 result table remains the source of statistical significance. This
script returns to the per-sample pseudobulk counts and asks whether each
selected peak is supported consistently across Healthy and MASLD samples or is
potentially driven by one or a few samples.

Expected input for each cell type
---------------------------------
1. ``data/derived/pseudobulk/<cell_type>_pseudobulk.csv``
   Rows are peak IDs and columns are the 12 sample pseudobulk counts.
2. ``results/peaQTL/differential_peaks/<cell_type>_deseq2_results.csv``
   Required columns: peak_id, baseMean, log2FoldChange, lfcSE, stat, pvalue,
   padj.

Healthy samples 4 and 6 are excluded from every calculation and plot.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import constants as C


REQUIRED_RESULT_COLUMNS = {
    "peak_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}

NUMERIC_RESULT_COLUMNS = (
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)

HEALTHY_CONDITION = "Healthy"
MASLD_CONDITION = "MASLD"

HEALTHY_COLOR = "#377eb8"
MASLD_COLOR = "#e6550d"
STRONG_COLOR = "#238b45"
MIXED_COLOR = "#fec44f"
WEAK_COLOR = "#d7301f"

STRONG_SUPPORT_LABEL = "Strong descriptive support"
MIXED_SUPPORT_LABEL = "Mixed descriptive support"
WEAK_SUPPORT_LABEL = "Weak/outlier-sensitive"

PEAK_ID_PATTERN = re.compile(r"^chr[^:]+:\d+-\d+$")


def configure_plotting() -> None:
    """Set consistent plotting defaults."""
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.frameon": False,
            "font.size": 9,
            "savefig.facecolor": "white",
        }
    )


def save_figure(figure: Figure, path: Path) -> None:
    """Save and close a figure."""
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=C.PLOT_DPI, bbox_inches="tight")
    plt.close(figure)


def validate_columns(
    dataframe: pd.DataFrame,
    required: set[str],
    path: Path,
) -> None:
    """Raise a useful error if a table is missing required columns."""
    missing = required.difference(dataframe.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: "
            f"{', '.join(sorted(missing))}"
        )


def slugify(value: str) -> str:
    """Create a filesystem-safe value."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return slug.strip("_") or "value"


def cell_type_from_result_path(path: Path) -> str:
    """Extract a cell-type name from a DESeq2 result filename."""
    if not path.name.endswith(C.DESEQ2_RESULTS_SUFFIX):
        raise ValueError(
            f"Unexpected DESeq2 filename {path.name!r}; expected suffix "
            f"{C.DESEQ2_RESULTS_SUFFIX!r}."
        )
    return path.name.removesuffix(C.DESEQ2_RESULTS_SUFFIX)


def discover_result_files() -> list[Path]:
    """Find all per-cell-type DESeq2 result files."""
    files = sorted(
        path
        for path in C.DESEQ2_RESULTS_DIR.glob(
            f"*{C.DESEQ2_RESULTS_SUFFIX}"
        )
        if path.is_file()
    )
    if not files:
        raise FileNotFoundError(
            f"No files ending with {C.DESEQ2_RESULTS_SUFFIX!r} were "
            f"found in {C.DESEQ2_RESULTS_DIR}."
        )
    return files


def load_deseq2_results(path: Path) -> pd.DataFrame:
    """Load and validate one DESeq2 result table."""
    dataframe = pd.read_csv(path)
    validate_columns(dataframe, REQUIRED_RESULT_COLUMNS, path)

    dataframe = dataframe.copy()
    dataframe["peak_id"] = dataframe["peak_id"].astype(str)

    for column in NUMERIC_RESULT_COLUMNS:
        dataframe[column] = pd.to_numeric(
            dataframe[column],
            errors="coerce",
        )

    for column in ("pvalue", "padj"):
        valid = (
            np.isfinite(dataframe[column])
            & dataframe[column].between(0.0, 1.0)
        )
        dataframe.loc[~valid, column] = np.nan

    if dataframe["peak_id"].duplicated().any():
        duplicates = dataframe.loc[
            dataframe["peak_id"].duplicated(), "peak_id"
        ].head(5)
        raise ValueError(
            f"{path} contains duplicate peak IDs, for example: "
            f"{', '.join(duplicates)}"
        )

    return dataframe


def detect_peak_column(dataframe: pd.DataFrame, path: Path) -> str:
    """Identify the peak-ID column in a pseudobulk CSV."""
    for candidate in ("peak_id", "peak", "region"):
        if candidate in dataframe.columns:
            return candidate

    first_column = str(dataframe.columns[0])
    values = dataframe.iloc[:, 0].dropna().astype(str)
    peak_fraction = (
        values.map(lambda value: bool(PEAK_ID_PATTERN.match(value))).mean()
        if not values.empty
        else 0.0
    )

    if first_column.startswith("Unnamed:") or peak_fraction >= 0.90:
        return first_column

    raise ValueError(
        f"Could not identify the peak-ID column in {path}. The first "
        f"column is {first_column!r}."
    )


def is_excluded_sample(sample_name: str) -> bool:
    """Return whether a sample is excluded from the analysis."""
    lowered = sample_name.lower()
    return any(
        token.lower() in lowered
        for token in C.DESEQ2_EXCLUDED_SAMPLE_TOKENS
    )


def infer_condition(sample_name: str) -> str:
    """Infer Healthy or MASLD from a pseudobulk sample column name."""
    lowered = sample_name.lower()

    if any(keyword.lower() in lowered for keyword in C.HEALTHY_KEYWORDS):
        return HEALTHY_CONDITION

    if any(keyword.lower() in lowered for keyword in C.MASLD_KEYWORDS):
        return MASLD_CONDITION

    raise ValueError(
        f"Could not infer a condition from sample column {sample_name!r}."
    )


def replicate_number(sample_name: str) -> int:
    """Extract the replicate number for deterministic sample ordering."""
    match = re.search(r"rep[_-]?(\d+)", sample_name, flags=re.IGNORECASE)
    return int(match.group(1)) if match else 10_000


def short_sample_label(sample_name: str, condition: str) -> str:
    """Create a compact label such as H1 or M3."""
    prefix = "H" if condition == HEALTHY_CONDITION else "M"
    replicate = replicate_number(sample_name)
    return f"{prefix}{replicate}" if replicate < 10_000 else sample_name


def load_pseudobulk_counts(
    path: Path,
) -> tuple[pd.DataFrame, pd.Series]:
    """Load one peak-by-sample pseudobulk matrix and remove exclusions."""
    if not path.exists():
        raise FileNotFoundError(f"Missing pseudobulk file: {path}")

    raw = pd.read_csv(path)
    peak_column = detect_peak_column(raw, path)

    counts = raw.set_index(peak_column)
    counts.index = counts.index.astype(str)
    counts.index.name = "peak_id"

    if counts.index.duplicated().any():
        duplicates = counts.index[counts.index.duplicated()].unique()[:5]
        raise ValueError(
            f"{path} contains duplicate peak IDs, for example: "
            f"{', '.join(duplicates)}"
        )

    kept_columns = [
        str(column)
        for column in counts.columns
        if not is_excluded_sample(str(column))
    ]
    excluded_columns = [
        str(column)
        for column in counts.columns
        if is_excluded_sample(str(column))
    ]

    if not excluded_columns:
        raise ValueError(
            f"No excluded samples were found in {path}. Expected columns "
            f"matching {C.DESEQ2_EXCLUDED_SAMPLE_TOKENS}."
        )

    counts = counts.loc[:, kept_columns].copy()

    for column in counts.columns:
        counts[column] = pd.to_numeric(counts[column], errors="coerce")

    if counts.isna().any().any():
        bad_columns = counts.columns[counts.isna().any()].tolist()
        raise ValueError(
            f"{path} contains missing or nonnumeric counts in columns: "
            f"{', '.join(bad_columns)}"
        )

    if (counts.to_numpy(dtype=float) < 0).any():
        raise ValueError(f"{path} contains negative pseudobulk counts.")

    conditions = pd.Series(
        {
            sample: infer_condition(sample)
            for sample in counts.columns
        },
        name="condition",
    )

    healthy_samples = sorted(
        conditions.index[conditions == HEALTHY_CONDITION],
        key=replicate_number,
    )
    masld_samples = sorted(
        conditions.index[conditions == MASLD_CONDITION],
        key=replicate_number,
    )

    if len(healthy_samples) != C.DESEQ2_EXPECTED_HEALTHY_SAMPLE_COUNT:
        raise ValueError(
            f"Expected {C.DESEQ2_EXPECTED_HEALTHY_SAMPLE_COUNT} included "
            f"healthy samples in {path}, found {len(healthy_samples)}: "
            f"{healthy_samples}"
        )

    if len(masld_samples) != C.DESEQ2_EXPECTED_MASLD_SAMPLE_COUNT:
        raise ValueError(
            f"Expected {C.DESEQ2_EXPECTED_MASLD_SAMPLE_COUNT} included "
            f"MASLD samples in {path}, found {len(masld_samples)}: "
            f"{masld_samples}"
        )

    ordered_samples = healthy_samples + masld_samples
    counts = counts.loc[:, ordered_samples]
    conditions = conditions.loc[ordered_samples]

    return counts, conditions


def positive_count_median_ratio_normalization(
    counts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.Series]:
    """Normalize counts using a positive-count median-ratio strategy.

    This is used only for sample-level visualization and descriptive support
    metrics. Statistical significance continues to come from the saved DESeq2
    result table.
    """
    sample_by_peak = counts.T.astype(float)
    array = sample_by_peak.to_numpy(dtype=float)

    positive = array > 0
    positive_counts_per_peak = positive.sum(axis=0)

    log_counts = np.zeros_like(array, dtype=float)
    np.log(array, out=log_counts, where=positive)

    log_geometric_means = np.full(array.shape[1], np.nan, dtype=float)
    valid_peaks = positive_counts_per_peak > 0
    log_geometric_means[valid_peaks] = (
        log_counts[:, valid_peaks].sum(axis=0)
        / positive_counts_per_peak[valid_peaks]
    )
    geometric_means = np.exp(log_geometric_means)

    size_factors = np.full(array.shape[0], np.nan, dtype=float)

    for sample_index in range(array.shape[0]):
        valid_ratios = (
            positive[sample_index]
            & np.isfinite(geometric_means)
            & (geometric_means > 0)
        )
        ratios = (
            array[sample_index, valid_ratios]
            / geometric_means[valid_ratios]
        )
        ratios = ratios[np.isfinite(ratios) & (ratios > 0)]

        if ratios.size:
            size_factors[sample_index] = float(np.median(ratios))

    invalid = ~np.isfinite(size_factors) | (size_factors <= 0)
    if invalid.any():
        library_sizes = array.sum(axis=1)
        if (library_sizes <= 0).any():
            bad_samples = sample_by_peak.index[library_sizes <= 0].tolist()
            raise ValueError(
                "Could not calculate normalization factors because these "
                f"included samples have zero total counts: {bad_samples}"
            )
        geometric_library_size = float(
            np.exp(np.mean(np.log(library_sizes)))
        )
        fallback = library_sizes / geometric_library_size
        size_factors[invalid] = fallback[invalid]

    # Center size factors at geometric mean one for easier interpretation.
    size_factors /= np.exp(np.mean(np.log(size_factors)))

    normalized_array = array / size_factors[:, None]

    normalized = pd.DataFrame(
        normalized_array,
        index=sample_by_peak.index,
        columns=sample_by_peak.columns,
    ).T

    size_factor_series = pd.Series(
        size_factors,
        index=sample_by_peak.index,
        name="size_factor",
    )

    return normalized, size_factor_series


def select_candidate_peaks(results: pd.DataFrame) -> pd.DataFrame:
    """Select peaks using the configured FDR and effect-size thresholds."""
    mask = (
        results["padj"].notna()
        & (results["padj"] < C.DESEQ2_DEFAULT_PADJ_THRESHOLD)
        & results["log2FoldChange"].notna()
        & (
            results["log2FoldChange"].abs()
            >= C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
        )
    )

    return (
        results.loc[mask]
        .sort_values(
            ["padj", "pvalue", "log2FoldChange"],
            ascending=[True, True, False],
            na_position="last",
        )
        .reset_index(drop=True)
    )


def safe_standard_deviation(values: np.ndarray) -> float:
    """Return sample standard deviation, or NaN with fewer than two values."""
    return float(np.std(values, ddof=1)) if values.size >= 2 else float("nan")


def safe_coefficient_of_variation(values: np.ndarray) -> float:
    """Return sample SD divided by the mean when defined."""
    mean = float(np.mean(values))
    if mean <= 0:
        return float("nan")
    return safe_standard_deviation(values) / mean


def median_absolute_deviation(values: np.ndarray) -> float:
    """Return the unscaled median absolute deviation."""
    median = float(np.median(values))
    return float(np.median(np.abs(values - median)))


def maximum_sample_share(values: np.ndarray) -> float:
    """Return the largest single-sample share of a group's total signal."""
    total = float(np.sum(values))
    if total <= 0:
        return float("nan")
    return float(np.max(values) / total)


def pairwise_masld_greater_probability(
    healthy_values: np.ndarray,
    masld_values: np.ndarray,
) -> float:
    """Estimate P(random MASLD value > random Healthy value), ties as 0.5."""
    comparisons = masld_values[:, None] - healthy_values[None, :]
    wins = np.sum(comparisons > 0)
    ties = np.sum(comparisons == 0)
    return float((wins + 0.5 * ties) / comparisons.size)


def hedges_g(
    healthy_values: np.ndarray,
    masld_values: np.ndarray,
) -> float:
    """Calculate Hedges' g for MASLD minus Healthy."""
    n_healthy = healthy_values.size
    n_masld = masld_values.size

    if n_healthy < 2 or n_masld < 2:
        return float("nan")

    healthy_variance = float(np.var(healthy_values, ddof=1))
    masld_variance = float(np.var(masld_values, ddof=1))

    degrees_of_freedom = n_healthy + n_masld - 2
    pooled_variance = (
        (n_healthy - 1) * healthy_variance
        + (n_masld - 1) * masld_variance
    ) / degrees_of_freedom

    if pooled_variance <= 0:
        return float("nan")

    pooled_sd = math.sqrt(pooled_variance)
    cohen_d = (
        float(np.mean(masld_values))
        - float(np.mean(healthy_values))
    ) / pooled_sd

    correction = 1.0 - 3.0 / (4.0 * degrees_of_freedom - 1.0)
    return float(correction * cohen_d)


def descriptive_log2_ratio(
    healthy_values: np.ndarray,
    masld_values: np.ndarray,
) -> float:
    """Return the log2 ratio of group means using a display pseudocount."""
    pseudocount = C.DESEQ2_SAMPLE_SUPPORT_PSEUDOCOUNT
    return float(
        np.log2(
            (float(np.mean(masld_values)) + pseudocount)
            / (float(np.mean(healthy_values)) + pseudocount)
        )
    )


def leave_one_out_effects(
    healthy_values: np.ndarray,
    masld_values: np.ndarray,
) -> np.ndarray:
    """Calculate descriptive group-mean log2 ratios after omitting each sample."""
    effects: list[float] = []

    for index in range(healthy_values.size):
        effects.append(
            descriptive_log2_ratio(
                np.delete(healthy_values, index),
                masld_values,
            )
        )

    for index in range(masld_values.size):
        effects.append(
            descriptive_log2_ratio(
                healthy_values,
                np.delete(masld_values, index),
            )
        )

    return np.asarray(effects, dtype=float)


def classify_sample_support(
    direction_support: float,
    loo_direction_stability: float,
    high_single_sample_concentration: bool,
) -> str:
    """Assign a descriptive consistency class; this is not a new p-value."""
    if (
        direction_support >= C.DESEQ2_SAMPLE_SUPPORT_STRONG_THRESHOLD
        and loo_direction_stability
        >= C.DESEQ2_SAMPLE_SUPPORT_LOO_STRONG_THRESHOLD
        and not high_single_sample_concentration
    ):
        return STRONG_SUPPORT_LABEL

    if (
        direction_support >= C.DESEQ2_SAMPLE_SUPPORT_MODERATE_THRESHOLD
        and loo_direction_stability
        >= C.DESEQ2_SAMPLE_SUPPORT_LOO_MODERATE_THRESHOLD
    ):
        return MIXED_SUPPORT_LABEL

    return WEAK_SUPPORT_LABEL


def calculate_candidate_metrics(
    cell_type: str,
    candidates: pd.DataFrame,
    raw_counts: pd.DataFrame,
    normalized_counts: pd.DataFrame,
    conditions: pd.Series,
) -> pd.DataFrame:
    """Calculate sample-level descriptive metrics for candidate peaks."""
    healthy_samples = conditions.index[
        conditions == HEALTHY_CONDITION
    ].tolist()
    masld_samples = conditions.index[
        conditions == MASLD_CONDITION
    ].tolist()

    records: list[dict[str, object]] = []

    for _, result in candidates.iterrows():
        peak_id = str(result["peak_id"])

        if peak_id not in raw_counts.index:
            continue

        healthy = normalized_counts.loc[peak_id, healthy_samples].to_numpy(
            dtype=float
        )
        masld = normalized_counts.loc[peak_id, masld_samples].to_numpy(
            dtype=float
        )
        healthy_raw = raw_counts.loc[peak_id, healthy_samples].to_numpy(
            dtype=float
        )
        masld_raw = raw_counts.loc[peak_id, masld_samples].to_numpy(
            dtype=float
        )

        probability_masld_greater = pairwise_masld_greater_probability(
            healthy,
            masld,
        )

        deseq_positive = float(result["log2FoldChange"]) > 0
        direction_support = (
            probability_masld_greater
            if deseq_positive
            else 1.0 - probability_masld_greater
        )

        complete_separation = (
            float(np.min(masld)) > float(np.max(healthy))
            if deseq_positive
            else float(np.min(healthy)) > float(np.max(masld))
        )

        full_descriptive_effect = descriptive_log2_ratio(healthy, masld)
        loo_effects = leave_one_out_effects(healthy, masld)

        expected_sign = 1.0 if deseq_positive else -1.0
        loo_direction_stability = float(
            np.mean(np.sign(loo_effects) == expected_sign)
        )

        healthy_share = maximum_sample_share(healthy)
        masld_share = maximum_sample_share(masld)
        high_single_sample_concentration = bool(
            (
                np.isfinite(healthy_share)
                and healthy_share
                > C.DESEQ2_SAMPLE_SUPPORT_MAX_SAMPLE_SHARE_WARNING
            )
            or (
                np.isfinite(masld_share)
                and masld_share
                > C.DESEQ2_SAMPLE_SUPPORT_MAX_SAMPLE_SHARE_WARNING
            )
        )

        support_class = classify_sample_support(
            direction_support,
            loo_direction_stability,
            high_single_sample_concentration,
        )

        record: dict[str, object] = {
            "cell_type": cell_type,
            "peak_id": peak_id,
            "baseMean": result["baseMean"],
            "log2FoldChange": result["log2FoldChange"],
            "lfcSE": result["lfcSE"],
            "stat": result["stat"],
            "pvalue": result["pvalue"],
            "padj": result["padj"],
            "healthy_n": healthy.size,
            "masld_n": masld.size,
            "healthy_mean_normalized": float(np.mean(healthy)),
            "masld_mean_normalized": float(np.mean(masld)),
            "healthy_median_normalized": float(np.median(healthy)),
            "masld_median_normalized": float(np.median(masld)),
            "healthy_sd_normalized": safe_standard_deviation(healthy),
            "masld_sd_normalized": safe_standard_deviation(masld),
            "healthy_cv_normalized": safe_coefficient_of_variation(healthy),
            "masld_cv_normalized": safe_coefficient_of_variation(masld),
            "healthy_mad_normalized": median_absolute_deviation(healthy),
            "masld_mad_normalized": median_absolute_deviation(masld),
            "healthy_min_normalized": float(np.min(healthy)),
            "healthy_max_normalized": float(np.max(healthy)),
            "masld_min_normalized": float(np.min(masld)),
            "masld_max_normalized": float(np.max(masld)),
            "healthy_nonzero_samples": int(np.sum(healthy_raw > 0)),
            "masld_nonzero_samples": int(np.sum(masld_raw > 0)),
            "healthy_nonzero_fraction": float(np.mean(healthy_raw > 0)),
            "masld_nonzero_fraction": float(np.mean(masld_raw > 0)),
            "descriptive_log2_mean_ratio": full_descriptive_effect,
            "hedges_g_masld_minus_healthy": hedges_g(healthy, masld),
            "pairwise_probability_masld_greater": probability_masld_greater,
            "direction_support_fraction": direction_support,
            "complete_group_separation": complete_separation,
            "loo_direction_stability_fraction": loo_direction_stability,
            "loo_log2_ratio_min": float(np.min(loo_effects)),
            "loo_log2_ratio_max": float(np.max(loo_effects)),
            "loo_max_abs_change_from_full": float(
                np.max(np.abs(loo_effects - full_descriptive_effect))
            ),
            "healthy_max_single_sample_share": healthy_share,
            "masld_max_single_sample_share": masld_share,
            "high_single_sample_concentration": high_single_sample_concentration,
            "sample_support_class": support_class,
        }

        for sample in conditions.index:
            record[f"raw__{sample}"] = float(raw_counts.loc[peak_id, sample])
            record[f"normalized__{sample}"] = float(
                normalized_counts.loc[peak_id, sample]
            )

        records.append(record)

    return pd.DataFrame.from_records(records)


def sample_support_color(support_class: str) -> str:
    """Map a descriptive support class to a plotting color."""
    return {
        STRONG_SUPPORT_LABEL: STRONG_COLOR,
        MIXED_SUPPORT_LABEL: MIXED_COLOR,
        WEAK_SUPPORT_LABEL: WEAK_COLOR,
    }[support_class]


def plot_group_points(
    axis: Axes,
    values: np.ndarray,
    samples: list[str],
    condition: str,
    x_position: float,
) -> None:
    """Plot all samples in one condition with mean and median summaries."""
    color = HEALTHY_COLOR if condition == HEALTHY_CONDITION else MASLD_COLOR
    transformed = np.log2(values + 1.0)
    offsets = np.linspace(-0.12, 0.12, len(values))

    axis.scatter(
        np.full(len(values), x_position) + offsets,
        transformed,
        s=34,
        color=color,
        edgecolor="black",
        linewidth=0.35,
        zorder=3,
    )

    for offset, value, sample in zip(offsets, transformed, samples):
        axis.annotate(
            short_sample_label(sample, condition),
            (x_position + offset, value),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            fontsize=6,
        )

    median = float(np.median(transformed))
    mean = float(np.mean(transformed))

    axis.plot(
        [x_position - 0.20, x_position + 0.20],
        [median, median],
        color="black",
        linewidth=2.0,
        zorder=4,
    )
    axis.scatter(
        [x_position],
        [mean],
        marker="D",
        s=28,
        color="white",
        edgecolor="black",
        linewidth=0.8,
        zorder=5,
    )


def plot_top_peak_sample_distributions(
    cell_type: str,
    metrics: pd.DataFrame,
    normalized_counts: pd.DataFrame,
    conditions: pd.Series,
    output_dir: Path,
) -> None:
    """Create a grid showing every included sample for the top candidates."""
    if metrics.empty:
        return

    top = metrics.sort_values(
        ["padj", "pvalue"],
        ascending=[True, True],
    ).head(C.DESEQ2_SAMPLE_SUPPORT_TOP_PEAKS)

    healthy_samples = conditions.index[
        conditions == HEALTHY_CONDITION
    ].tolist()
    masld_samples = conditions.index[
        conditions == MASLD_CONDITION
    ].tolist()

    n_columns = 3
    n_rows = math.ceil(len(top) / n_columns)

    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(5.2 * n_columns, 4.2 * n_rows),
        squeeze=False,
        constrained_layout=True,
    )

    for axis, (_, row) in zip(axes.flat, top.iterrows()):
        peak_id = str(row["peak_id"])
        healthy = normalized_counts.loc[peak_id, healthy_samples].to_numpy(
            dtype=float
        )
        masld = normalized_counts.loc[peak_id, masld_samples].to_numpy(
            dtype=float
        )

        plot_group_points(
            axis,
            healthy,
            healthy_samples,
            HEALTHY_CONDITION,
            0.0,
        )
        plot_group_points(
            axis,
            masld,
            masld_samples,
            MASLD_CONDITION,
            1.0,
        )

        axis.set_xticks([0.0, 1.0])
        axis.set_xticklabels([HEALTHY_CONDITION, MASLD_CONDITION])
        axis.set_xlim(-0.35, 1.35)
        axis.set_ylabel("log2(normalized count + 1)")
        axis.grid(axis="y", alpha=0.18)

        axis.set_title(
            f"{peak_id}\n"
            f"padj={row['padj']:.2g}, "
            f"DESeq2 LFC={row['log2FoldChange']:.2f}\n"
            f"pairwise support={100 * row['direction_support_fraction']:.0f}%, "
            f"LOO={100 * row['loo_direction_stability_fraction']:.0f}%",
            fontsize=8,
            color=sample_support_color(str(row["sample_support_class"])),
        )

    for axis in axes.flat[len(top) :]:
        axis.set_visible(False)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=HEALTHY_COLOR,
            markeredgecolor="black",
            label="Healthy sample",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=MASLD_COLOR,
            markeredgecolor="black",
            label="MASLD sample",
        ),
        Line2D([0], [0], color="black", linewidth=2, label="Group median"),
        Line2D(
            [0],
            [0],
            marker="D",
            color="none",
            markerfacecolor="white",
            markeredgecolor="black",
            label="Group mean",
        ),
    ]

    figure.legend(
        handles=legend_handles,
        loc="outside lower center",
        ncols=4,
    )
    figure.suptitle(
        f"{cell_type}: sample-level accessibility of top candidate peaks\n"
        "All included biological samples are shown; values are for visualization",
        fontsize=14,
        fontweight="bold",
    )

    save_figure(
        figure,
        output_dir
        / f"{slugify(cell_type)}_top_peak_sample_distributions."
        f"{C.DESEQ2_PLOT_FORMAT}",
    )


def plot_top_peak_heatmap(
    cell_type: str,
    metrics: pd.DataFrame,
    normalized_counts: pd.DataFrame,
    conditions: pd.Series,
    output_dir: Path,
) -> None:
    """Create a row-standardized heatmap of top candidate peaks."""
    if metrics.empty:
        return

    top = metrics.sort_values(
        ["padj", "pvalue"],
        ascending=[True, True],
    ).head(C.DESEQ2_SAMPLE_SUPPORT_HEATMAP_PEAKS)

    peak_ids = top["peak_id"].astype(str).tolist()
    samples = conditions.index.tolist()

    log_counts = np.log2(
        normalized_counts.loc[peak_ids, samples].to_numpy(dtype=float) + 1.0
    )
    row_means = np.mean(log_counts, axis=1, keepdims=True)
    row_sds = np.std(log_counts, axis=1, ddof=1, keepdims=True)
    row_sds[~np.isfinite(row_sds) | (row_sds == 0)] = 1.0
    row_z_scores = (log_counts - row_means) / row_sds
    row_z_scores = np.clip(row_z_scores, -2.5, 2.5)

    figure, axis = plt.subplots(
        figsize=(10, max(5.0, 0.27 * len(peak_ids) + 2.3)),
        constrained_layout=True,
    )

    image = axis.imshow(
        row_z_scores,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        vmin=-2.5,
        vmax=2.5,
    )

    sample_labels = [
        short_sample_label(sample, str(conditions.loc[sample]))
        for sample in samples
    ]

    axis.set_xticks(np.arange(len(samples)))
    axis.set_xticklabels(sample_labels)
    axis.set_yticks(np.arange(len(peak_ids)))
    axis.set_yticklabels(peak_ids, fontsize=7)
    axis.set_xlabel("Biological sample")
    axis.set_ylabel("Candidate peak")

    n_healthy = int(np.sum(conditions == HEALTHY_CONDITION))
    axis.axvline(n_healthy - 0.5, color="black", linewidth=2)

    for tick, sample in zip(axis.get_xticklabels(), samples):
        tick.set_color(
            HEALTHY_COLOR
            if conditions.loc[sample] == HEALTHY_CONDITION
            else MASLD_COLOR
        )

    axis.set_title(
        f"{cell_type}: top candidate peaks across individual samples\n"
        "Row z-scores of log2(normalized count + 1)",
    )

    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Within-peak standardized accessibility")

    axis.legend(
        handles=[
            Patch(color=HEALTHY_COLOR, label="Healthy"),
            Patch(color=MASLD_COLOR, label="MASLD"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
    )

    save_figure(
        figure,
        output_dir
        / f"{slugify(cell_type)}_top_peak_sample_heatmap."
        f"{C.DESEQ2_PLOT_FORMAT}",
    )


def plot_effect_vs_sample_support(
    cell_type: str,
    metrics: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Compare DESeq2 effect size with descriptive sample-level support."""
    if metrics.empty:
        return

    figure, axis = plt.subplots(figsize=(9.5, 6), constrained_layout=True)

    for support_class in (
        STRONG_SUPPORT_LABEL,
        MIXED_SUPPORT_LABEL,
        WEAK_SUPPORT_LABEL,
    ):
        subset = metrics[metrics["sample_support_class"] == support_class]
        if subset.empty:
            continue

        axis.scatter(
            subset["log2FoldChange"],
            subset["direction_support_fraction"],
            s=48,
            color=sample_support_color(support_class),
            alpha=0.82,
            edgecolor="black",
            linewidth=0.35,
            label=f"{support_class} ({len(subset)})",
        )

    axis.axvline(0.0, color="black", linewidth=1)
    axis.axhline(
        C.DESEQ2_SAMPLE_SUPPORT_MODERATE_THRESHOLD,
        color="black",
        linestyle=":",
        alpha=0.6,
    )
    axis.axhline(
        C.DESEQ2_SAMPLE_SUPPORT_STRONG_THRESHOLD,
        color="black",
        linestyle="--",
        alpha=0.6,
    )

    top_to_label = metrics.sort_values(
        ["padj", "pvalue"],
        ascending=[True, True],
    ).head(C.DESEQ2_SAMPLE_SUPPORT_MAX_PEAK_LABELS)

    for _, row in top_to_label.iterrows():
        axis.annotate(
            str(row["peak_id"]),
            (
                float(row["log2FoldChange"]),
                float(row["direction_support_fraction"]),
            ),
            xytext=(3, 4),
            textcoords="offset points",
            fontsize=6,
        )

    axis.set_ylim(0.45, 1.03)
    axis.set_xlabel("DESeq2 log2 fold change (MASLD / Healthy)")
    axis.set_ylabel(
        "Pairwise sample support for the DESeq2 direction"
    )
    axis.set_title(
        f"{cell_type}: statistical effect versus sample-level consistency"
    )
    axis.grid(alpha=0.18)
    axis.legend(loc="best")

    save_figure(
        figure,
        output_dir
        / f"{slugify(cell_type)}_effect_vs_sample_support."
        f"{C.DESEQ2_PLOT_FORMAT}",
    )


def plot_all_cell_type_support_summary(
    all_metrics: pd.DataFrame,
    all_cell_types: Iterable[str],
    output_dir: Path,
) -> None:
    """Summarize descriptive support classes across all cell types."""
    cell_types = list(all_cell_types)
    support_order = [
        STRONG_SUPPORT_LABEL,
        MIXED_SUPPORT_LABEL,
        WEAK_SUPPORT_LABEL,
    ]

    counts = (
        all_metrics.groupby(["cell_type", "sample_support_class"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=cell_types, fill_value=0)
        .reindex(columns=support_order, fill_value=0)
    )

    totals = counts.sum(axis=1)
    ordered_cell_types = totals.sort_values(ascending=False).index.tolist()
    counts = counts.loc[ordered_cell_types]

    figure, axis = plt.subplots(
        figsize=(max(10.0, 1.25 * len(counts)), 5.8),
        constrained_layout=True,
    )

    x_positions = np.arange(len(counts))
    bottom = np.zeros(len(counts), dtype=float)

    for support_class in support_order:
        values = counts[support_class].to_numpy(dtype=float)
        axis.bar(
            x_positions,
            values,
            bottom=bottom,
            color=sample_support_color(support_class),
            label=support_class,
        )
        bottom += values

    for index, total in enumerate(bottom):
        axis.text(
            index,
            total + 0.8,
            str(int(total)),
            ha="center",
            va="bottom",
            fontsize=8,
        )

    axis.set_xticks(x_positions)
    axis.set_xticklabels(counts.index, rotation=30, ha="right")
    axis.set_ylabel("Number of DESeq2-selected peaks")
    axis.set_title(
        "Sample-level consistency of DESeq2-selected peaks\n"
        "Descriptive support classes do not replace FDR significance"
    )
    axis.legend()
    axis.grid(axis="y", alpha=0.18)

    save_figure(
        figure,
        output_dir
        / f"all_cell_types_candidate_sample_support."
        f"{C.DESEQ2_PLOT_FORMAT}",
    )


def process_cell_type(
    result_path: Path,
) -> tuple[pd.DataFrame, pd.Series]:
    """Calculate metrics and plots for one cell type."""
    cell_type = cell_type_from_result_path(result_path)
    pseudobulk_path = (
        C.PSEUDOBULK_OUTPUT_DIR / f"{cell_type}_pseudobulk.csv"
    )

    print(f"Processing {cell_type}...")

    results = load_deseq2_results(result_path)
    candidates = select_candidate_peaks(results)

    raw_counts, conditions = load_pseudobulk_counts(pseudobulk_path)
    normalized_counts, size_factors = (
        positive_count_median_ratio_normalization(raw_counts)
    )

    missing_peaks = candidates.loc[
        ~candidates["peak_id"].isin(raw_counts.index), "peak_id"
    ].tolist()
    if missing_peaks:
        print(
            f"  Warning: {len(missing_peaks)} candidate peak(s) were not "
            "found in the pseudobulk matrix and will be skipped."
        )

    candidates = candidates[
        candidates["peak_id"].isin(raw_counts.index)
    ].reset_index(drop=True)

    metrics = calculate_candidate_metrics(
        cell_type,
        candidates,
        raw_counts,
        normalized_counts,
        conditions,
    )

    cell_output_dir = (
        C.DESEQ2_SAMPLE_SUPPORT_PER_CELL_PLOTS_DIR / slugify(cell_type)
    )

    if metrics.empty:
        print("  No peaks pass the configured candidate criterion.")
    else:
        plot_top_peak_sample_distributions(
            cell_type,
            metrics,
            normalized_counts,
            conditions,
            cell_output_dir,
        )
        plot_top_peak_heatmap(
            cell_type,
            metrics,
            normalized_counts,
            conditions,
            cell_output_dir,
        )
        plot_effect_vs_sample_support(
            cell_type,
            metrics,
            cell_output_dir,
        )

    size_factors = size_factors.rename_axis("sample").reset_index()
    size_factors.insert(0, "cell_type", cell_type)
    size_factors["condition"] = size_factors["sample"].map(conditions)

    print(f"  Candidates analyzed: {len(metrics):,}")
    return metrics, size_factors


def main() -> None:
    """Create sample-level candidate summaries and plots."""
    configure_plotting()

    C.DESEQ2_SAMPLE_SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
    C.DESEQ2_SAMPLE_SUPPORT_PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    C.DESEQ2_SAMPLE_SUPPORT_PER_CELL_PLOTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    result_files = discover_result_files()

    metrics_frames: list[pd.DataFrame] = []
    size_factor_frames: list[pd.DataFrame] = []
    cell_types: list[str] = []

    for result_path in result_files:
        cell_type = cell_type_from_result_path(result_path)
        cell_types.append(cell_type)

        metrics, size_factors = process_cell_type(result_path)
        metrics_frames.append(metrics)
        size_factor_frames.append(size_factors)

    nonempty_metrics = [frame for frame in metrics_frames if not frame.empty]
    all_metrics = (
        pd.concat(nonempty_metrics, ignore_index=True)
        if nonempty_metrics
        else pd.DataFrame()
    )

    all_size_factors = pd.concat(size_factor_frames, ignore_index=True)
    all_size_factors.to_csv(
        C.DESEQ2_SAMPLE_SUPPORT_SIZE_FACTORS_CSV,
        index=False,
        float_format="%.10g",
    )

    if all_metrics.empty:
        print("No candidate peaks were found in any cell type.")
        return

    all_metrics = all_metrics.sort_values(
        ["cell_type", "padj", "pvalue"],
        ascending=[True, True, True],
    ).reset_index(drop=True)

    all_metrics.to_csv(
        C.DESEQ2_SAMPLE_SUPPORT_SUMMARY_CSV,
        index=False,
        float_format="%.10g",
    )

    plot_all_cell_type_support_summary(
        all_metrics,
        cell_types,
        C.DESEQ2_SAMPLE_SUPPORT_PLOTS_DIR,
    )

    print("Finished sample-support analysis.")
    print(
        f"Candidate metrics:\n  {C.DESEQ2_SAMPLE_SUPPORT_SUMMARY_CSV}"
    )
    print(
        f"Normalization factors:\n  "
        f"{C.DESEQ2_SAMPLE_SUPPORT_SIZE_FACTORS_CSV}"
    )
    print(f"Plots:\n  {C.DESEQ2_SAMPLE_SUPPORT_PLOTS_DIR}")


if __name__ == "__main__":
    main()
