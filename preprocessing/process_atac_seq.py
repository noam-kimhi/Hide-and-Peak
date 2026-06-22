from __future__ import annotations

import gzip
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import sparse
from scipy.io import mmread

from constants import (
    ATAC_SEQ_DIRS,
    BARCODES_SUFFIX,
    CONDITION_QC_OUTPUT_PATH,
    EXPECTED_SAMPLE_COUNT,
    EXPECTED_SAMPLES_PER_CONDITION,
    FRAGMENTS_SUFFIX,
    HEALTHY_KEYWORDS,
    MASLD_KEYWORDS,
    MATRIX_SUFFIX,
    PAIRWISE_PEAK_OVERLAP_OUTPUT_PATH,
    PEAKS_SUFFIX,
    SAMPLE_PEAK_SHARING_OUTPUT_PATH,
    SAMPLE_QC_OUTPUT_PATH,
)


@dataclass(frozen=True)
class SampleFiles:
    """Paths to the core snATAC-seq files belonging to one sample."""

    sample_dir: Path
    sample_id: str
    barcodes_path: Path
    peaks_path: Path
    matrix_path: Path
    fragments_path: Path
    fragments_index_path: Path


def find_single_file(sample_dir: Path, suffix: str) -> Path:
    """
    Find exactly one file whose name ends with the requested suffix.

    Args:
        sample_dir: Directory containing one snATAC-seq sample.
        suffix: Required filename suffix.

    Returns:
        Path to the matching file.

    Raises:
        FileNotFoundError: If no matching file is found.
        ValueError: If multiple matching files are found.
    """
    matches = sorted(sample_dir.glob(f"*{suffix}"))

    if not matches:
        raise FileNotFoundError(
            f"No file ending with {suffix!r} found in {sample_dir}"
        )

    if len(matches) > 1:
        raise ValueError(
            f"Expected one file ending with {suffix!r} in {sample_dir}, "
            f"but found {len(matches)}: {matches}"
        )

    return matches[0]


def collect_sample_files(sample_dir: Path) -> SampleFiles:
    """
    Collect all required files for one snATAC-seq sample.

    Args:
        sample_dir: Directory containing one sample.

    Returns:
        Paths and identifiers for the sample.

    Raises:
        FileNotFoundError: If a required file or fragment index is missing.
    """
    sample_dir = Path(sample_dir)
    sample_id = sample_dir.name

    fragments_path = find_single_file(sample_dir, FRAGMENTS_SUFFIX)
    fragments_index_path = Path(f"{fragments_path}.tbi")

    if not fragments_index_path.exists():
        raise FileNotFoundError(
            f"Missing Tabix index for {sample_id}: "
            f"{fragments_index_path}"
        )

    return SampleFiles(
        sample_dir=sample_dir,
        sample_id=sample_id,
        barcodes_path=find_single_file(sample_dir, BARCODES_SUFFIX),
        peaks_path=find_single_file(sample_dir, PEAKS_SUFFIX),
        matrix_path=find_single_file(sample_dir, MATRIX_SUFFIX),
        fragments_path=fragments_path,
        fragments_index_path=fragments_index_path,
    )


def collect_all_sample_files(
    sample_dirs: Iterable[Path],
) -> list[SampleFiles]:
    """
    Collect and validate file paths for all samples.

    Args:
        sample_dirs: Iterable of sample directories.

    Returns:
        List containing one SampleFiles object per sample.

    Raises:
        ValueError: If duplicate sample identifiers are found.
    """
    sample_files = [
        collect_sample_files(Path(sample_dir))
        for sample_dir in sample_dirs
    ]

    sample_ids = [sample.sample_id for sample in sample_files]

    if len(sample_ids) != len(set(sample_ids)):
        duplicated_ids = sorted(
            sample_id
            for sample_id in set(sample_ids)
            if sample_ids.count(sample_id) > 1
        )
        raise ValueError(
            f"Duplicate sample identifiers found: {duplicated_ids}"
        )

    return sample_files


def infer_condition(sample_id: str) -> str:
    """
    Infer the biological condition from a sample identifier.

    Args:
        sample_id: Sample directory name or sample identifier.

    Returns:
        One of 'healthy', 'MASLD', or 'unknown'.
    """
    sample_id_lower = sample_id.lower()

    matches_healthy = any(
        keyword.lower() in sample_id_lower
        for keyword in HEALTHY_KEYWORDS
    )
    matches_masld = any(
        keyword.lower() in sample_id_lower
        for keyword in MASLD_KEYWORDS
    )

    if matches_healthy and matches_masld:
        raise ValueError(
            f"Sample {sample_id!r} matches both healthy and MASLD keywords."
        )

    if matches_healthy:
        return "healthy"

    if matches_masld:
        return "MASLD"

    return "unknown"


def validate_cohort(sample_files: Sequence[SampleFiles]) -> None:
    """
    Validate the number of samples and their inferred conditions.

    Unexpected cohort composition produces warnings rather than stopping the
    analysis, because the data may intentionally differ from the paper.

    Args:
        sample_files: Files for all cohort samples.
    """
    if len(sample_files) != EXPECTED_SAMPLE_COUNT:
        warnings.warn(
            f"Expected {EXPECTED_SAMPLE_COUNT} samples, "
            f"but found {len(sample_files)}.",
            stacklevel=2,
        )

    conditions = pd.Series(
        [infer_condition(sample.sample_id) for sample in sample_files],
        dtype="string",
    )
    condition_counts = conditions.value_counts().to_dict()

    for condition, expected_count in EXPECTED_SAMPLES_PER_CONDITION.items():
        observed_count = int(condition_counts.get(condition, 0))

        if observed_count != expected_count:
            warnings.warn(
                f"Expected {expected_count} {condition} samples, "
                f"but found {observed_count}.",
                stacklevel=2,
            )

    unknown_samples = [
        sample.sample_id
        for sample in sample_files
        if infer_condition(sample.sample_id) == "unknown"
    ]

    if unknown_samples:
        warnings.warn(
            f"Could not infer condition for samples: {unknown_samples}",
            stacklevel=2,
        )


def read_barcodes(path: Path) -> pd.Series:
    """
    Read cell barcodes from a gzipped TSV file.

    Args:
        path: Path to a barcodes.tsv.gz file.

    Returns:
        Series containing one barcode per row.

    Raises:
        ValueError: If the file is empty or contains empty barcodes.
    """
    with gzip.open(path, "rt") as file:
        barcodes = pd.Series(
            [line.rstrip("\n\r") for line in file],
            dtype="string",
            name="barcode",
        )

    if barcodes.empty:
        raise ValueError(f"Barcode file is empty: {path}")

    empty_mask = barcodes.isna() | barcodes.str.strip().eq("")

    if empty_mask.any():
        raise ValueError(
            f"Barcode file {path} contains "
            f"{int(empty_mask.sum())} empty barcode rows."
        )

    return barcodes


def read_bed3(path: Path) -> pd.DataFrame:
    """
    Read and validate the first three columns of a BED file.

    BED intervals use zero-based, half-open coordinates:

        [start, end)

    Args:
        path: Path to a BED or BED.gz file.

    Returns:
        DataFrame with columns chrom, start, and end.

    Raises:
        ValueError: If coordinates are missing or invalid.
    """
    peaks = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={
            "chrom": "string",
            "start": "int64",
            "end": "int64",
        },
        compression="infer",
    )

    if peaks.empty:
        raise ValueError(f"Peak BED file is empty: {path}")

    if peaks.isna().any(axis=None):
        raise ValueError(f"Peak BED file contains missing values: {path}")

    negative_start_mask = peaks["start"] < 0
    invalid_length_mask = peaks["end"] <= peaks["start"]

    if negative_start_mask.any():
        raise ValueError(
            f"Peak BED file {path} contains "
            f"{int(negative_start_mask.sum())} negative start coordinates."
        )

    if invalid_length_mask.any():
        raise ValueError(
            f"Peak BED file {path} contains "
            f"{int(invalid_length_mask.sum())} intervals with end <= start."
        )

    return peaks


def load_sparse_matrix(path: Path) -> sparse.csc_matrix:
    """
    Load a Matrix Market matrix and convert it to CSC format.

    Explicitly stored zero values are removed before statistics are computed.

    Args:
        path: Path to a matrix.mtx or matrix.mtx.gz file.

    Returns:
        Sparse peak-by-cell matrix in CSC format.

    Raises:
        ValueError: If matrix entries are negative or non-finite.
    """
    matrix = mmread(path)

    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)

    matrix = matrix.tocsc()
    matrix.eliminate_zeros()

    if matrix.data.size > 0:
        if not np.isfinite(matrix.data).all():
            raise ValueError(
                f"Matrix contains non-finite values: {path}"
            )

        if np.any(matrix.data < 0):
            raise ValueError(
                f"Matrix contains negative counts: {path}"
            )

    return matrix


def are_matrix_values_integer_like(
    matrix: sparse.spmatrix,
    tolerance: float = 1e-8,
) -> bool:
    """
    Test whether all stored matrix values are approximately integers.

    Args:
        matrix: Sparse count matrix.
        tolerance: Maximum allowed distance from the nearest integer.

    Returns:
        True if all stored values are integer-like.
    """
    if matrix.data.size == 0:
        return True

    distance_from_integer = np.abs(
        matrix.data - np.rint(matrix.data)
    )
    return bool(np.all(distance_from_integer <= tolerance))


def chromosome_blocks_are_contiguous(peaks: pd.DataFrame) -> bool:
    """
    Check whether each chromosome appears in one contiguous BED block.

    This avoids requiring lexicographic chromosome ordering, which can differ
    from natural chromosome ordering.

    Args:
        peaks: BED3 DataFrame.

    Returns:
        True if no chromosome appears, disappears, and then appears again.
    """
    chromosomes = peaks["chrom"].to_numpy()

    if chromosomes.size <= 1:
        return True

    block_starts = np.concatenate(
        (
            np.array([True]),
            chromosomes[1:] != chromosomes[:-1],
        )
    )
    block_chromosomes = chromosomes[block_starts]

    return len(block_chromosomes) == len(set(block_chromosomes))


def starts_are_sorted_within_chromosomes(
    peaks: pd.DataFrame,
) -> bool:
    """
    Check whether starts are nondecreasing within every chromosome.

    Args:
        peaks: BED3 DataFrame.

    Returns:
        True if peak starts are sorted within every chromosome.
    """
    return all(
        chromosome_peaks["start"].is_monotonic_increasing
        for _, chromosome_peaks in peaks.groupby(
            "chrom",
            sort=False,
            observed=True,
        )
    )


def summarize_numeric_values(
    values: NDArray[np.number],
    prefix: str,
) -> dict[str, float]:
    """
    Compute a standard collection of descriptive statistics.

    Args:
        values: One-dimensional numeric array.
        prefix: Prefix used for output-column names.

    Returns:
        Dictionary containing min, mean, quantiles, and max.
    """
    values = np.asarray(values).reshape(-1)

    quantiles = np.quantile(
        values,
        [0.05, 0.25, 0.50, 0.75, 0.95],
    )

    return {
        f"{prefix}_min": float(np.min(values)),
        f"{prefix}_mean": float(np.mean(values)),
        f"{prefix}_p05": float(quantiles[0]),
        f"{prefix}_p25": float(quantiles[1]),
        f"{prefix}_median": float(quantiles[2]),
        f"{prefix}_p75": float(quantiles[3]),
        f"{prefix}_p95": float(quantiles[4]),
        f"{prefix}_max": float(np.max(values)),
    }


def summarize_sample(sample_files: SampleFiles) -> dict[str, object]:
    """
    Compute structural and count-based QC statistics for one sample.

    The matrix is expected to have peaks as rows and cells as columns.

    Args:
        sample_files: Paths belonging to one sample.

    Returns:
        Dictionary containing one row of sample-level QC information.

    Raises:
        ValueError: If files disagree about matrix dimensions.
    """
    barcodes = read_barcodes(sample_files.barcodes_path)
    peaks = read_bed3(sample_files.peaks_path)
    matrix = load_sparse_matrix(sample_files.matrix_path)

    n_peaks, n_cells = matrix.shape
    n_barcodes = len(barcodes)

    if n_cells != n_barcodes:
        raise ValueError(
            f"Cell mismatch in {sample_files.sample_id}: "
            f"matrix has {n_cells} columns, "
            f"but barcode file has {n_barcodes} rows."
        )

    if n_peaks != len(peaks):
        raise ValueError(
            f"Peak mismatch in {sample_files.sample_id}: "
            f"matrix has {n_peaks} rows, "
            f"but peak file has {len(peaks)} rows."
        )

    duplicate_barcode_count = int(barcodes.duplicated().sum())
    duplicate_peak_count = int(
        peaks.duplicated(
            subset=["chrom", "start", "end"]
        ).sum()
    )

    peak_lengths = (
        peaks["end"].to_numpy()
        - peaks["start"].to_numpy()
    )

    peak_counts_per_cell = np.asarray(
        matrix.sum(axis=0)
    ).reshape(-1)

    accessible_peaks_per_cell = np.asarray(
        matrix.getnnz(axis=0)
    ).reshape(-1)

    cells_per_peak = np.asarray(
        matrix.getnnz(axis=1)
    ).reshape(-1)

    total_entries = n_peaks * n_cells
    nonzero_entries = int(matrix.nnz)
    density = (
        nonzero_entries / total_entries
        if total_entries > 0
        else 0.0
    )
    sparsity_fraction = 1.0 - density

    result: dict[str, object] = {
        "sample_id": sample_files.sample_id,
        "condition": infer_condition(sample_files.sample_id),
        "n_cells": n_cells,
        "n_unique_barcodes": int(barcodes.nunique()),
        "duplicate_barcode_count": duplicate_barcode_count,
        "n_peaks": n_peaks,
        "n_unique_peak_intervals": int(
            peaks.drop_duplicates(
                subset=["chrom", "start", "end"]
            ).shape[0]
        ),
        "duplicate_peak_count": duplicate_peak_count,
        "n_chromosomes": int(peaks["chrom"].nunique()),
        "chromosome_blocks_contiguous": (
            chromosome_blocks_are_contiguous(peaks)
        ),
        "starts_sorted_within_chromosomes": (
            starts_are_sorted_within_chromosomes(peaks)
        ),
        "total_peak_bp": int(peak_lengths.sum()),
        "nonzero_entries": nonzero_entries,
        "matrix_density": float(density),
        "sparsity": float(sparsity_fraction),
        "matrix_values_integer_like": (
            are_matrix_values_integer_like(matrix)
        ),
        "matrix_min_stored_value": (
            float(matrix.data.min())
            if matrix.data.size > 0
            else 0.0
        ),
        "matrix_max_stored_value": (
            float(matrix.data.max())
            if matrix.data.size > 0
            else 0.0
        ),
        "barcodes_file_size_mib": (
            sample_files.barcodes_path.stat().st_size
            / (1024**2)
        ),
        "peaks_file_size_mib": (
            sample_files.peaks_path.stat().st_size
            / (1024**2)
        ),
        "matrix_file_size_mib": (
            sample_files.matrix_path.stat().st_size
            / (1024**2)
        ),
        "fragments_file_size_mib": (
            sample_files.fragments_path.stat().st_size
            / (1024**2)
        ),
        "fragments_index_exists": (
            sample_files.fragments_index_path.exists()
        ),
        "barcodes_path": sample_files.barcodes_path.name,
        "peaks_path": sample_files.peaks_path.name,
        "matrix_path": sample_files.matrix_path.name,
        "fragments_path": sample_files.fragments_path.name,
        "fragments_index_path": (
            sample_files.fragments_index_path.name
        ),
    }

    result.update(
        summarize_numeric_values(
            peak_lengths,
            prefix="peak_length",
        )
    )
    result.update(
        summarize_numeric_values(
            peak_counts_per_cell,
            prefix="peak_counts_per_cell",
        )
    )
    result.update(
        summarize_numeric_values(
            accessible_peaks_per_cell,
            prefix="accessible_peaks_per_cell",
        )
    )
    result.update(
        summarize_numeric_values(
            cells_per_peak,
            prefix="cells_per_peak",
        )
    )

    return result


def summarize_all_samples(
    sample_files: Sequence[SampleFiles],
) -> pd.DataFrame:
    """
    Compute QC statistics for all samples.

    Args:
        sample_files: Validated files for every sample.

    Returns:
        DataFrame with one row per sample.
    """
    summaries: list[dict[str, object]] = []

    for index, sample in enumerate(sample_files, start=1):
        print(
            f"[{index}/{len(sample_files)}] "
            f"Processing {sample.sample_id}..."
        )
        summaries.append(summarize_sample(sample))

    return (
        pd.DataFrame(summaries)
        .sort_values(
            by=["condition", "sample_id"],
            kind="stable",
        )
        .reset_index(drop=True)
    )


def create_condition_qc_summary(
    sample_qc: pd.DataFrame,
) -> pd.DataFrame:
    """
    Summarize major QC properties separately for each condition.

    Args:
        sample_qc: Sample-level QC table.

    Returns:
        DataFrame with one row per condition.
    """
    return (
        sample_qc.groupby(
            "condition",
            dropna=False,
            observed=True,
        )
        .agg(
            n_samples=("sample_id", "size"),
            total_cells=("n_cells", "sum"),
            median_cells_per_sample=("n_cells", "median"),
            min_cells_per_sample=("n_cells", "min"),
            max_cells_per_sample=("n_cells", "max"),
            median_peaks_per_sample=("n_peaks", "median"),
            min_peaks_per_sample=("n_peaks", "min"),
            max_peaks_per_sample=("n_peaks", "max"),
            median_peak_length=(
                "peak_length_median",
                "median",
            ),
            median_peak_counts_per_cell=(
                "peak_counts_per_cell_median",
                "median",
            ),
            median_accessible_peaks_per_cell=(
                "accessible_peaks_per_cell_median",
                "median",
            ),
            median_sparsity=("sparsity", "median"),
        )
        .reset_index()
    )


def prepare_peak_set_for_overlap(
    peaks: pd.DataFrame,
) -> pd.DataFrame:
    """
    Prepare a BED3 table for interval-overlap operations.

    Exact duplicate intervals are removed. The original peak BED files are not
    modified.

    Args:
        peaks: BED3 DataFrame.

    Returns:
        Sorted, duplicate-free BED3 DataFrame.
    """
    return (
        peaks.drop_duplicates(
            subset=["chrom", "start", "end"]
        )
        .sort_values(
            by=["chrom", "start", "end"],
            kind="stable",
        )
        .reset_index(drop=True)
    )


def peaks_overlapping_any_target(
    query_peaks: pd.DataFrame,
    target_peaks: pd.DataFrame,
) -> NDArray[np.bool_]:
    """
    Identify query peaks overlapping at least one target peak.

    An overlap requires at least one shared genomic base:

        query_start < target_end
        target_start < query_end

    The implementation is vectorized by chromosome. For each target
    chromosome, it stores the prefix maximum of target ends. This allows each
    query interval to be checked using binary search.

    Args:
        query_peaks: Intervals for which overlap status is requested.
        target_peaks: Intervals against which queries are tested.

    Returns:
        Boolean array aligned with query_peaks.
    """
    overlap_mask = np.zeros(
        len(query_peaks),
        dtype=bool,
    )

    target_by_chromosome: dict[
        str,
        tuple[
            NDArray[np.int64],
            NDArray[np.int64],
        ],
    ] = {}

    for chromosome, chromosome_peaks in target_peaks.groupby(
        "chrom",
        sort=False,
        observed=True,
    ):
        target_starts = chromosome_peaks[
            "start"
        ].to_numpy(dtype=np.int64)
        target_ends = chromosome_peaks[
            "end"
        ].to_numpy(dtype=np.int64)

        order = np.argsort(
            target_starts,
            kind="stable",
        )
        target_starts = target_starts[order]
        target_ends = target_ends[order]

        prefix_max_ends = np.maximum.accumulate(
            target_ends
        )

        target_by_chromosome[str(chromosome)] = (
            target_starts,
            prefix_max_ends,
        )

    query_chromosomes = query_peaks[
        "chrom"
    ].astype(str).to_numpy()

    for chromosome in np.unique(query_chromosomes):
        target_data = target_by_chromosome.get(chromosome)

        if target_data is None:
            continue

        target_starts, target_prefix_max_ends = target_data

        query_positions = np.flatnonzero(
            query_chromosomes == chromosome
        )
        chromosome_queries = query_peaks.iloc[
            query_positions
        ]

        query_starts = chromosome_queries[
            "start"
        ].to_numpy(dtype=np.int64)
        query_ends = chromosome_queries[
            "end"
        ].to_numpy(dtype=np.int64)

        last_candidate_indices = (
            np.searchsorted(
                target_starts,
                query_ends,
                side="left",
            )
            - 1
        )

        has_candidate = last_candidate_indices >= 0

        chromosome_overlap = np.zeros(
            len(query_positions),
            dtype=bool,
        )

        valid_candidate_indices = last_candidate_indices[
            has_candidate
        ]

        chromosome_overlap[has_candidate] = (
            target_prefix_max_ends[
                valid_candidate_indices
            ]
            > query_starts[has_candidate]
        )

        overlap_mask[query_positions] = chromosome_overlap

    return overlap_mask


def load_peak_sets(
    sample_files: Sequence[SampleFiles],
) -> dict[str, pd.DataFrame]:
    """
    Load duplicate-free peak sets for all samples.

    Args:
        sample_files: Files belonging to every sample.

    Returns:
        Mapping from sample ID to prepared BED3 DataFrame.
    """
    peak_sets: dict[str, pd.DataFrame] = {}

    for index, sample in enumerate(sample_files, start=1):
        print(
            f"[{index}/{len(sample_files)}] "
            f"Loading peaks for {sample.sample_id}..."
        )
        peaks = read_bed3(sample.peaks_path)
        peak_sets[sample.sample_id] = (
            prepare_peak_set_for_overlap(peaks)
        )

    return peak_sets


def calculate_peak_sharing(
    peak_sets: Mapping[str, pd.DataFrame],
    sample_conditions: Mapping[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate pairwise and cohort-wide peak sharing.

    Pairwise overlap is directional because sample A and sample B can contain
    different numbers of peaks. For example, the fraction of A overlapping B
    need not equal the fraction of B overlapping A.

    A sample-specific peak is defined as a peak that overlaps no peak from any
    other sample by even one base pair.

    Args:
        peak_sets: Mapping from sample ID to its prepared peaks.
        sample_conditions: Mapping from sample ID to condition.

    Returns:
        Tuple containing:
            1. Pairwise directional-overlap table.
            2. Per-sample shared-versus-unique peak summary.
    """
    sample_ids = sorted(peak_sets)

    shared_with_any_masks = {
        sample_id: np.zeros(
            len(peak_sets[sample_id]),
            dtype=bool,
        )
        for sample_id in sample_ids
    }

    pairwise_rows: list[dict[str, object]] = []

    for first_index, first_sample_id in enumerate(sample_ids):
        first_peaks = peak_sets[first_sample_id]

        for second_sample_id in sample_ids[first_index + 1:]:
            second_peaks = peak_sets[second_sample_id]

            first_overlaps_second = peaks_overlapping_any_target(
                query_peaks=first_peaks,
                target_peaks=second_peaks,
            )
            second_overlaps_first = peaks_overlapping_any_target(
                query_peaks=second_peaks,
                target_peaks=first_peaks,
            )

            shared_with_any_masks[
                first_sample_id
            ] |= first_overlaps_second
            shared_with_any_masks[
                second_sample_id
            ] |= second_overlaps_first

            n_first_overlapping_second = int(
                first_overlaps_second.sum()
            )
            n_second_overlapping_first = int(
                second_overlaps_first.sum()
            )

            pairwise_rows.append(
                {
                    "sample_a": first_sample_id,
                    "condition_a": sample_conditions[
                        first_sample_id
                    ],
                    "n_peaks_a": len(first_peaks),
                    "sample_b": second_sample_id,
                    "condition_b": sample_conditions[
                        second_sample_id
                    ],
                    "n_peaks_b": len(second_peaks),
                    "n_peaks_a_overlapping_b": (
                        n_first_overlapping_second
                    ),
                    "fraction_peaks_a_overlapping_b": (
                        n_first_overlapping_second
                        / len(first_peaks)
                    ),
                    "n_peaks_b_overlapping_a": (
                        n_second_overlapping_first
                    ),
                    "fraction_peaks_b_overlapping_a": (
                        n_second_overlapping_first
                        / len(second_peaks)
                    ),
                    "mean_directional_overlap_fraction": (
                        0.5
                        * (
                            n_first_overlapping_second
                            / len(first_peaks)
                            + n_second_overlapping_first
                            / len(second_peaks)
                        )
                    ),
                    "same_condition": (
                        sample_conditions[first_sample_id]
                        == sample_conditions[second_sample_id]
                    ),
                }
            )

    sharing_rows: list[dict[str, object]] = []

    for sample_id in sample_ids:
        n_peaks = len(peak_sets[sample_id])
        shared_mask = shared_with_any_masks[sample_id]

        n_shared = int(shared_mask.sum())
        n_unique = n_peaks - n_shared

        sharing_rows.append(
            {
                "sample_id": sample_id,
                "condition": sample_conditions[sample_id],
                "n_unique_peak_intervals": n_peaks,
                "n_peaks_overlapping_any_other_sample": n_shared,
                "fraction_peaks_overlapping_any_other_sample": (
                    n_shared / n_peaks
                ),
                "n_peaks_unique_to_sample": n_unique,
                "fraction_peaks_unique_to_sample": (
                    n_unique / n_peaks
                ),
            }
        )

    pairwise_overlap = pd.DataFrame(pairwise_rows)
    sharing_summary = pd.DataFrame(sharing_rows)

    return pairwise_overlap, sharing_summary


def save_dataframe(
    dataframe: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Save a DataFrame as a CSV file.

    Args:
        dataframe: Table to save.
        output_path: Output CSV path.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    dataframe.to_csv(output_path, index=False)


def print_sample_qc_summary(
    sample_qc: pd.DataFrame,
) -> None:
    """
    Print the most useful sample-level QC columns.

    Args:
        sample_qc: Sample-level QC table.
    """
    display_columns = [
        "sample_id",
        "condition",
        "n_cells",
        "n_peaks",
        "peak_counts_per_cell_median",
        "accessible_peaks_per_cell_median",
        "duplicate_barcode_count",
        "duplicate_peak_count",
        "sparsity",
    ]

    print("\nSample-level QC:")
    print(
        sample_qc[display_columns].to_string(
            index=False
        )
    )


def print_peak_sharing_summary(
    peak_sharing: pd.DataFrame,
) -> None:
    """
    Print sample-level peak-sharing results.

    Args:
        peak_sharing: Per-sample peak-sharing table.
    """
    display_columns = [
        "sample_id",
        "condition",
        "n_unique_peak_intervals",
        "n_peaks_overlapping_any_other_sample",
        "fraction_peaks_overlapping_any_other_sample",
        "n_peaks_unique_to_sample",
    ]

    print("\nPeak sharing across samples:")
    print(
        peak_sharing[display_columns].to_string(
            index=False
        )
    )


def main() -> None:
    """Run structural, sample-level, and cross-sample snATAC-seq QC."""
    sample_files = collect_all_sample_files(
        ATAC_SEQ_DIRS
    )
    validate_cohort(sample_files)

    sample_qc = summarize_all_samples(sample_files)
    condition_qc = create_condition_qc_summary(
        sample_qc
    )

    sample_conditions = {
        sample.sample_id: infer_condition(
            sample.sample_id
        )
        for sample in sample_files
    }

    peak_sets = load_peak_sets(sample_files)
    pairwise_peak_overlap, peak_sharing = (
        calculate_peak_sharing(
            peak_sets=peak_sets,
            sample_conditions=sample_conditions,
        )
    )

    save_dataframe(
        sample_qc,
        SAMPLE_QC_OUTPUT_PATH,
    )
    save_dataframe(
        condition_qc,
        CONDITION_QC_OUTPUT_PATH,
    )
    save_dataframe(
        pairwise_peak_overlap,
        PAIRWISE_PEAK_OVERLAP_OUTPUT_PATH,
    )
    save_dataframe(
        peak_sharing,
        SAMPLE_PEAK_SHARING_OUTPUT_PATH,
    )

    print_sample_qc_summary(sample_qc)
    print_peak_sharing_summary(peak_sharing)

    print("\nSaved outputs:")
    print(f"  Sample QC:       {SAMPLE_QC_OUTPUT_PATH}")
    print(f"  Condition QC:    {CONDITION_QC_OUTPUT_PATH}")
    print(
        "  Pairwise overlap: "
        f"{PAIRWISE_PEAK_OVERLAP_OUTPUT_PATH}"
    )
    print(
        "  Peak sharing:     "
        f"{SAMPLE_PEAK_SHARING_OUTPUT_PATH}"
    )


if __name__ == "__main__":
    main()