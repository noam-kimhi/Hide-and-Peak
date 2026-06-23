#!/usr/bin/env python3
"""Create high-confidence differential-accessibility BED files per cell type.

For each configured cell type, this script combines:

1. the saved DESeq2 differential-accessibility results; and
2. the per-sample pseudobulk ATAC-seq count matrix.

A peak is written to ``<cell_type>_significant_peaks.bed.gz`` only when it
passes the DESeq2 adjusted-p-value and effect-size thresholds *and* the
individual Healthy and MASLD pseudobulk samples strongly support the DESeq2
estimated direction.

Healthy samples 4 and 6 are excluded from all sample-level calculations.
The output BED files contain exactly three columns with no header:
``chrom``, ``start``, and ``end``.

This script intentionally produces all configured cell-type files, including
valid empty ``.bed.gz`` files for cell types with no qualifying peaks.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

import numpy as np
import pandas as pd

import constants as C


HEALTHY_CONDITION = "Healthy"
MASLD_CONDITION = "MASLD"

REQUIRED_DESEQ2_COLUMNS = {
    "peak_id",
    "log2FoldChange",
    "pvalue",
    "padj",
}

NUMERIC_DESEQ2_COLUMNS = (
    "log2FoldChange",
    "pvalue",
    "padj",
)

PEAK_ID_PATTERN = re.compile(
    r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$"
)


@dataclass(frozen=True)
class PeakSupportMetrics:
    """Sample-level evidence used by the high-confidence peak filter."""

    median_direction_agrees: bool
    direction_support_fraction: float
    loo_direction_stability_fraction: float
    healthy_max_single_sample_share: float
    masld_max_single_sample_share: float
    high_single_sample_concentration: bool


@dataclass(frozen=True)
class CellTypeSelectionResult:
    """Selection result and audit counts for one cell type."""

    cell_type: str
    candidate_count: int
    median_direction_count: int
    pairwise_support_count: int
    loo_stability_count: int
    concentration_count: int
    final_peak_ids: tuple[str, ...]


@dataclass(frozen=True)
class BedInterval:
    """A genomic interval parsed from a peak identifier."""

    chrom: str
    start: int
    end: int


def validate_required_constants() -> None:
    """Fail early with a clear message when a required constant is missing."""
    required_names = (
        "UNFILTERED_SIG_PEAKS_DIR",
        "SIGNIFICANT_PEAK_CELL_TYPES",
        "DESEQ2_EXCLUDED_SAMPLE_TOKENS",
        "DESEQ2_EXPECTED_HEALTHY_SAMPLE_COUNT",
        "DESEQ2_EXPECTED_MASLD_SAMPLE_COUNT",
        "DESEQ2_SAMPLE_SUPPORT_PSEUDOCOUNT",
        "DESEQ2_SAMPLE_SUPPORT_STRONG_THRESHOLD",
        "DESEQ2_SAMPLE_SUPPORT_LOO_STRONG_THRESHOLD",
        "DESEQ2_SAMPLE_SUPPORT_MAX_SAMPLE_SHARE_WARNING",
    )

    missing = [name for name in required_names if not hasattr(C, name)]
    if missing:
        raise AttributeError(
            "constants.py is missing the following required constants: "
            + ", ".join(missing)
        )


def validate_cell_types(cell_types: Iterable[str]) -> tuple[str, ...]:
    """Return exactly seven unique, non-empty cell-type names."""
    normalized = tuple(str(cell_type).strip() for cell_type in cell_types)

    if len(normalized) != 7:
        raise ValueError(
            "SIGNIFICANT_PEAK_CELL_TYPES must contain exactly seven cell "
            f"types; found {len(normalized)}: {normalized}"
        )

    if any(not cell_type for cell_type in normalized):
        raise ValueError(
            "SIGNIFICANT_PEAK_CELL_TYPES contains an empty cell-type name."
        )

    if len(set(normalized)) != len(normalized):
        raise ValueError(
            "SIGNIFICANT_PEAK_CELL_TYPES contains duplicate cell types: "
            f"{normalized}"
        )

    return normalized


def detect_peak_column(dataframe: pd.DataFrame, path: Path) -> str:
    """Identify the peak-ID column in a pseudobulk CSV."""
    for candidate in ("peak_id", "peak", "region"):
        if candidate in dataframe.columns:
            return candidate

    if dataframe.empty or len(dataframe.columns) == 0:
        raise ValueError(f"{path} is empty or has no columns.")

    first_column = str(dataframe.columns[0])
    values = dataframe.iloc[:, 0].dropna().astype(str)
    peak_fraction = (
        values.map(lambda value: bool(PEAK_ID_PATTERN.fullmatch(value))).mean()
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
    """Return whether a sample matches one of the configured exclusions."""
    lowered = sample_name.lower()
    return any(
        str(token).lower() in lowered
        for token in C.DESEQ2_EXCLUDED_SAMPLE_TOKENS
    )


def infer_condition(sample_name: str) -> str:
    """Infer Healthy or MASLD from a pseudobulk sample column name."""
    lowered = sample_name.lower()

    healthy_match = any(
        str(keyword).lower() in lowered for keyword in C.HEALTHY_KEYWORDS
    )
    masld_match = any(
        str(keyword).lower() in lowered for keyword in C.MASLD_KEYWORDS
    )

    if healthy_match == masld_match:
        raise ValueError(
            "Could not uniquely infer a condition from sample column "
            f"{sample_name!r}."
        )

    return HEALTHY_CONDITION if healthy_match else MASLD_CONDITION


def replicate_number(sample_name: str) -> int:
    """Extract a replicate number for deterministic sample ordering."""
    match = re.search(r"rep[_-]?(\d+)", sample_name, flags=re.IGNORECASE)
    return int(match.group(1)) if match else 10_000


def load_deseq2_results(path: Path) -> pd.DataFrame:
    """Load and validate one cell type's DESeq2 result table."""
    if not path.exists():
        raise FileNotFoundError(f"Missing DESeq2 result file: {path}")

    dataframe = pd.read_csv(path)
    missing = REQUIRED_DESEQ2_COLUMNS.difference(dataframe.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: "
            f"{', '.join(sorted(missing))}"
        )

    dataframe = dataframe.copy()
    dataframe["peak_id"] = dataframe["peak_id"].astype(str)

    if dataframe["peak_id"].duplicated().any():
        examples = dataframe.loc[
            dataframe["peak_id"].duplicated(), "peak_id"
        ].head(5)
        raise ValueError(
            f"{path} contains duplicate peak IDs, for example: "
            f"{', '.join(examples)}"
        )

    for column in NUMERIC_DESEQ2_COLUMNS:
        dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")

    for column in ("pvalue", "padj"):
        valid = (
            np.isfinite(dataframe[column])
            & dataframe[column].between(0.0, 1.0)
        )
        dataframe.loc[~valid, column] = np.nan

    return dataframe


def load_pseudobulk_counts(
    path: Path,
) -> tuple[pd.DataFrame, pd.Series]:
    """Load a peak-by-sample pseudobulk matrix and remove excluded samples."""
    if not path.exists():
        raise FileNotFoundError(f"Missing pseudobulk file: {path}")

    raw = pd.read_csv(path)
    peak_column = detect_peak_column(raw, path)

    counts = raw.set_index(peak_column)
    counts.index = counts.index.astype(str)
    counts.index.name = "peak_id"

    if counts.index.duplicated().any():
        examples = counts.index[counts.index.duplicated()].unique()[:5]
        raise ValueError(
            f"{path} contains duplicate peak IDs, for example: "
            f"{', '.join(examples)}"
        )

    included_columns = [
        str(column)
        for column in counts.columns
        if not is_excluded_sample(str(column))
    ]
    counts = counts.loc[:, included_columns].copy()

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
        {sample: infer_condition(sample) for sample in counts.columns},
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

    expected_healthy = C.DESEQ2_EXPECTED_HEALTHY_SAMPLE_COUNT
    expected_masld = C.DESEQ2_EXPECTED_MASLD_SAMPLE_COUNT

    if len(healthy_samples) != expected_healthy:
        raise ValueError(
            f"Expected {expected_healthy} included Healthy samples in "
            f"{path}, found {len(healthy_samples)}: {healthy_samples}"
        )

    if len(masld_samples) != expected_masld:
        raise ValueError(
            f"Expected {expected_masld} included MASLD samples in {path}, "
            f"found {len(masld_samples)}: {masld_samples}"
        )

    ordered_samples = healthy_samples + masld_samples
    return counts.loc[:, ordered_samples], conditions.loc[ordered_samples]


def positive_count_median_ratio_normalization(
    counts: pd.DataFrame,
) -> pd.DataFrame:
    """Normalize counts using the same positive-count strategy as before.

    These normalized values are used only for descriptive sample-support
    filtering. Statistical significance comes exclusively from the saved
    DESeq2 result table.
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
        ratios = array[sample_index, valid_ratios] / geometric_means[
            valid_ratios
        ]
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

    size_factors /= np.exp(np.mean(np.log(size_factors)))
    normalized_array = array / size_factors[:, None]

    return pd.DataFrame(
        normalized_array,
        index=sample_by_peak.index,
        columns=sample_by_peak.columns,
    ).T


def select_deseq2_candidates(results: pd.DataFrame) -> pd.DataFrame:
    """Select peaks passing the configured FDR and effect-size thresholds."""
    mask = (
        results["padj"].notna()
        & (results["padj"] < C.DESEQ2_DEFAULT_PADJ_THRESHOLD)
        & results["log2FoldChange"].notna()
        & (
            results["log2FoldChange"].abs()
            >= C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
        )
    )

    return results.loc[mask].copy()


def pairwise_probability_masld_greater(
    healthy_values: np.ndarray,
    masld_values: np.ndarray,
) -> float:
    """Return P(MASLD > Healthy) over all sample pairs, with ties as 0.5."""
    comparisons = masld_values[:, None] - healthy_values[None, :]
    wins = int(np.sum(comparisons > 0))
    ties = int(np.sum(comparisons == 0))
    return float((wins + 0.5 * ties) / comparisons.size)


def descriptive_log2_mean_ratio(
    healthy_values: np.ndarray,
    masld_values: np.ndarray,
) -> float:
    """Return log2(mean MASLD / mean Healthy) using a small pseudocount."""
    pseudocount = float(C.DESEQ2_SAMPLE_SUPPORT_PSEUDOCOUNT)
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
    """Recalculate the descriptive effect after omitting each sample once."""
    effects: list[float] = []

    for index in range(healthy_values.size):
        effects.append(
            descriptive_log2_mean_ratio(
                np.delete(healthy_values, index),
                masld_values,
            )
        )

    for index in range(masld_values.size):
        effects.append(
            descriptive_log2_mean_ratio(
                healthy_values,
                np.delete(masld_values, index),
            )
        )

    return np.asarray(effects, dtype=float)


def maximum_sample_share(values: np.ndarray) -> float:
    """Return the largest single-sample share of a condition's total signal."""
    total = float(np.sum(values))
    if total <= 0:
        return float("nan")
    return float(np.max(values) / total)


def peak_support_metrics(
    log2_fold_change: float,
    healthy_values: np.ndarray,
    masld_values: np.ndarray,
) -> PeakSupportMetrics:
    """Calculate all sample-level criteria for one DESeq2 candidate peak."""
    positive_direction = log2_fold_change > 0

    healthy_median = float(np.median(healthy_values))
    masld_median = float(np.median(masld_values))
    median_direction_agrees = (
        masld_median > healthy_median
        if positive_direction
        else healthy_median > masld_median
    )

    probability_masld_greater = pairwise_probability_masld_greater(
        healthy_values,
        masld_values,
    )
    direction_support_fraction = (
        probability_masld_greater
        if positive_direction
        else 1.0 - probability_masld_greater
    )

    loo_effect_values = leave_one_out_effects(
        healthy_values,
        masld_values,
    )
    expected_sign = 1.0 if positive_direction else -1.0
    loo_direction_stability_fraction = float(
        np.mean(np.sign(loo_effect_values) == expected_sign)
    )

    healthy_share = maximum_sample_share(healthy_values)
    masld_share = maximum_sample_share(masld_values)
    concentration_threshold = float(
        C.DESEQ2_SAMPLE_SUPPORT_MAX_SAMPLE_SHARE_WARNING
    )
    high_single_sample_concentration = bool(
        (
            np.isfinite(healthy_share)
            and healthy_share > concentration_threshold
        )
        or (
            np.isfinite(masld_share)
            and masld_share > concentration_threshold
        )
    )

    return PeakSupportMetrics(
        median_direction_agrees=median_direction_agrees,
        direction_support_fraction=direction_support_fraction,
        loo_direction_stability_fraction=loo_direction_stability_fraction,
        healthy_max_single_sample_share=healthy_share,
        masld_max_single_sample_share=masld_share,
        high_single_sample_concentration=high_single_sample_concentration,
    )


def qualifies_as_high_confidence(metrics: PeakSupportMetrics) -> bool:
    """Return whether one candidate passes every sample-support condition."""
    return bool(
        metrics.median_direction_agrees
        and metrics.direction_support_fraction
        >= C.DESEQ2_SAMPLE_SUPPORT_STRONG_THRESHOLD
        and metrics.loo_direction_stability_fraction
        >= C.DESEQ2_SAMPLE_SUPPORT_LOO_STRONG_THRESHOLD
        and not metrics.high_single_sample_concentration
    )


def process_cell_type(cell_type: str) -> CellTypeSelectionResult:
    """Select high-confidence differential peaks for one cell type."""
    result_path = (
        C.DESEQ2_RESULTS_DIR
        / f"{cell_type}{C.DESEQ2_RESULTS_SUFFIX}"
    )
    pseudobulk_path = (
        C.PSEUDOBULK_OUTPUT_DIR / f"{cell_type}_pseudobulk.csv"
    )

    results = load_deseq2_results(result_path)
    counts, conditions = load_pseudobulk_counts(pseudobulk_path)
    normalized_counts = positive_count_median_ratio_normalization(counts)
    candidates = select_deseq2_candidates(results)

    missing_peaks = sorted(set(candidates["peak_id"]) - set(counts.index))
    if missing_peaks:
        examples = ", ".join(missing_peaks[:5])
        raise ValueError(
            f"{len(missing_peaks)} DESeq2 candidate peaks for {cell_type} "
            "are absent from the pseudobulk matrix, for example: "
            f"{examples}"
        )

    healthy_samples = conditions.index[
        conditions == HEALTHY_CONDITION
    ].tolist()
    masld_samples = conditions.index[
        conditions == MASLD_CONDITION
    ].tolist()

    median_direction_count = 0
    pairwise_support_count = 0
    loo_stability_count = 0
    concentration_count = 0
    final_peak_ids: list[str] = []

    for row in candidates.itertuples(index=False):
        peak_id = str(row.peak_id)
        log2_fold_change = float(row.log2FoldChange)

        healthy_values = normalized_counts.loc[
            peak_id, healthy_samples
        ].to_numpy(dtype=float)
        masld_values = normalized_counts.loc[
            peak_id, masld_samples
        ].to_numpy(dtype=float)

        metrics = peak_support_metrics(
            log2_fold_change,
            healthy_values,
            masld_values,
        )

        median_direction_count += int(metrics.median_direction_agrees)
        pairwise_support_count += int(
            metrics.direction_support_fraction
            >= C.DESEQ2_SAMPLE_SUPPORT_STRONG_THRESHOLD
        )
        loo_stability_count += int(
            metrics.loo_direction_stability_fraction
            >= C.DESEQ2_SAMPLE_SUPPORT_LOO_STRONG_THRESHOLD
        )
        concentration_count += int(
            not metrics.high_single_sample_concentration
        )

        if qualifies_as_high_confidence(metrics):
            final_peak_ids.append(peak_id)

    return CellTypeSelectionResult(
        cell_type=cell_type,
        candidate_count=len(candidates),
        median_direction_count=median_direction_count,
        pairwise_support_count=pairwise_support_count,
        loo_stability_count=loo_stability_count,
        concentration_count=concentration_count,
        final_peak_ids=tuple(final_peak_ids),
    )


def parse_peak_id(peak_id: str) -> BedInterval:
    """Parse ``chr:start-end`` without changing its BED coordinates."""
    match = PEAK_ID_PATTERN.fullmatch(peak_id)
    if match is None:
        raise ValueError(
            f"Peak ID {peak_id!r} does not match 'chrom:start-end'."
        )

    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))

    if start < 0:
        raise ValueError(f"Peak {peak_id!r} has a negative start coordinate.")
    if end <= start:
        raise ValueError(
            f"Peak {peak_id!r} has end <= start and is not a valid BED "
            "interval."
        )

    return BedInterval(chrom=chrom, start=start, end=end)


def chromosome_sort_key(chrom: str) -> tuple[int, int | str]:
    """Return a natural GRCh38-like chromosome ordering key."""
    lowered = chrom.lower()
    core = lowered[3:] if lowered.startswith("chr") else lowered

    if core.isdigit():
        return 0, int(core)
    if core == "x":
        return 1, 23
    if core == "y":
        return 1, 24
    if core in {"m", "mt"}:
        return 1, 25
    return 2, core


def sorted_intervals(peak_ids: Iterable[str]) -> list[BedInterval]:
    """Parse, deduplicate, and sort selected peak identifiers."""
    intervals = {parse_peak_id(peak_id) for peak_id in peak_ids}
    return sorted(
        intervals,
        key=lambda interval: (
            chromosome_sort_key(interval.chrom),
            interval.start,
            interval.end,
        ),
    )


def write_bed_gz(path: Path, intervals: Iterable[BedInterval]) -> None:
    """Write a three-column, headerless gzip-compressed BED file."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        bed_handle: TextIO = handle
        for interval in intervals:
            bed_handle.write(
                f"{interval.chrom}\t{interval.start}\t{interval.end}\n"
            )


def print_summary(results: Iterable[CellTypeSelectionResult]) -> None:
    """Print an audit-friendly selection summary."""
    rows = []
    for result in results:
        rows.append(
            {
                "cell_type": result.cell_type,
                "DESeq2_candidates": result.candidate_count,
                "median_direction": result.median_direction_count,
                "pairwise_support": result.pairwise_support_count,
                "LOO_stability": result.loo_stability_count,
                "no_single_sample_dominance": result.concentration_count,
                "final_significant": len(result.final_peak_ids),
            }
        )

    summary = pd.DataFrame(rows)
    print("\nSignificant-peak selection summary:")
    print(summary.to_string(index=False))
    print(
        "\nNote: intermediate columns count candidates passing each "
        "criterion individually; final_significant requires all criteria."
    )


def main() -> None:
    """Create all seven high-confidence significant-peak BED files."""
    validate_required_constants()
    cell_types = validate_cell_types(C.SIGNIFICANT_PEAK_CELL_TYPES)

    # Compute every result before writing anything, so an input error cannot
    # leave a misleading partially updated output directory.
    selections = [process_cell_type(cell_type) for cell_type in cell_types]

    output_dir = Path(C.UNFILTERED_SIG_PEAKS_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    for selection in selections:
        output_path = (
            output_dir
            / f"{selection.cell_type}_significant_peaks.bed.gz"
        )
        intervals = sorted_intervals(selection.final_peak_ids)
        write_bed_gz(output_path, intervals)
        print(
            f"Wrote {len(intervals):,} peaks to {output_path}"
        )

    print_summary(selections)
    print(f"\nFinished. Output directory: {output_dir}")


if __name__ == "__main__":
    main()
