from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import bioframe as bf
from numpy.typing import NDArray

from constants import (
    ATAC_SEQ_DIRS,
    BLACKLIST_BED_PATH,
    CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
    CONSENSUS_PEAKS_OUTPUT_PATH,
    CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH,
    MERGE_BOOKENDED_PEAKS,
    MIN_REPRODUCIBLE_SAMPLE_SUPPORT,
    PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH,
    PEAKS_SUFFIX,
    REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
    REPRODUCIBLE_CONSENSUS_PEAKS_OUTPUT_PATH,
)


@dataclass(frozen=True)
class SamplePeakFile:
    """
    Store metadata for one sample-specific peak BED file.
    """

    sample_id: str
    peaks_path: Path
    sample_index: int


def find_single_file(sample_dir: Path, suffix: str) -> Path:
    """
    Find exactly one file in a directory whose name ends with a given suffix.

    :param sample_dir: Directory containing one snATAC-seq sample.
    :param suffix: Required filename suffix.
    :return: Path to the single matching file.
    """
    matches = sorted(sample_dir.glob(f"*{suffix}"))

    if not matches:
        raise FileNotFoundError(
            f"No file ending with {suffix!r} was found in {sample_dir}."
        )

    if len(matches) > 1:
        raise ValueError(
            f"Expected one file ending with {suffix!r} in {sample_dir}, "
            f"but found {len(matches)}: {matches}"
        )

    return matches[0]


def collect_sample_peak_files(
    sample_dirs: Iterable[Path],
) -> list[SamplePeakFile]:
    """
    Collect the peak BED file and sample identifier for every sample directory.

    :param sample_dirs: Directories containing the snATAC-seq samples.
    :return: One SamplePeakFile object per sample.
    """
    sample_peak_files: list[SamplePeakFile] = []

    for sample_index, sample_dir_value in enumerate(sample_dirs):
        sample_dir = Path(sample_dir_value)

        sample_peak_files.append(
            SamplePeakFile(
                sample_id=sample_dir.name,
                peaks_path=find_single_file(
                    sample_dir=sample_dir,
                    suffix=PEAKS_SUFFIX,
                ),
                sample_index=sample_index,
            )
        )

    if not sample_peak_files:
        raise ValueError("No sample directories were provided.")

    sample_ids = [
        sample_peak_file.sample_id
        for sample_peak_file in sample_peak_files
    ]

    if len(sample_ids) != len(set(sample_ids)):
        duplicated_sample_ids = sorted(
            sample_id
            for sample_id in set(sample_ids)
            if sample_ids.count(sample_id) > 1
        )

        raise ValueError(
            "Duplicate sample identifiers were found: "
            f"{duplicated_sample_ids}"
        )

    return sample_peak_files


def read_bed3(path: Path) -> pd.DataFrame:
    """
    Read and validate the first three columns of a BED file.

    BED intervals are interpreted as zero-based, half-open intervals:
    [start, end).

    :param path: Path to a BED or BED.gz file.
    :return: Validated table with chrom, start, and end columns.
    """
    intervals = pd.read_csv(
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

    if intervals.empty:
        raise ValueError(f"BED file is empty: {path}")

    if intervals.isna().any(axis=None):
        raise ValueError(
            f"BED file contains missing values: {path}"
        )

    negative_start_mask = intervals["start"] < 0
    invalid_interval_mask = (
        intervals["end"] <= intervals["start"]
    )

    if negative_start_mask.any():
        raise ValueError(
            f"{path} contains "
            f"{int(negative_start_mask.sum())} "
            "negative start coordinates."
        )

    if invalid_interval_mask.any():
        raise ValueError(
            f"{path} contains "
            f"{int(invalid_interval_mask.sum())} "
            "intervals with end <= start."
        )

    return intervals


def read_peak_bed(
    sample_peak_file: SamplePeakFile,
) -> pd.DataFrame:
    """
    Read one sample's peak BED file and attach its sample index.

    :param sample_peak_file: Metadata and BED path for one sample.
    :return: Peak table with chrom, start, end, and sample_index columns.
    """
    peaks = read_bed3(sample_peak_file.peaks_path)
    peaks["sample_index"] = np.int16(
        sample_peak_file.sample_index
    )

    return peaks


def read_blacklist(path: Path) -> pd.DataFrame:
    """
    Read, deduplicate, and sort the ENCODE blacklist intervals.

    :param path: Path to the ENCODE blacklist BED file.
    :return: Sorted, duplicate-free blacklist intervals.
    """
    blacklist = read_bed3(path)

    return (
        blacklist.drop_duplicates(
            subset=["chrom", "start", "end"]
        )
        .sort_values(
            by=["chrom", "start", "end"],
            kind="stable",
        )
        .reset_index(drop=True)
    )


def chromosome_sort_key(
    chromosome: str,
) -> tuple[int, int, str]:
    """
    Create a natural genomic ordering key for a chromosome or contig name.

    Canonical chromosomes are ordered as chr1-chr22, chrX, chrY, and chrM.
    Non-canonical contigs are placed afterward in lexicographic order.

    :param chromosome: Chromosome or contig name.
    :return: Tuple suitable for deterministic chromosome sorting.
    """
    normalized = chromosome.removeprefix("chr")

    if normalized.isdigit():
        chromosome_number = int(normalized)

        if 1 <= chromosome_number <= 22:
            return 0, chromosome_number, ""

    special_chromosome_order = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }

    normalized_upper = normalized.upper()

    if normalized_upper in special_chromosome_order:
        return (
            0,
            special_chromosome_order[normalized_upper],
            "",
        )

    return 1, 0, chromosome


def sort_genomic_intervals(
    intervals: pd.DataFrame,
) -> pd.DataFrame:
    """
    Sort genomic intervals using natural chromosome and coordinate order.

    All additional columns in the input table are preserved.

    :param intervals: Table containing chrom, start, and end columns.
    :return: Genomically sorted interval table.
    """
    sorted_intervals = intervals.copy()

    ordered_chromosomes = sorted(
        sorted_intervals["chrom"]
        .drop_duplicates()
        .astype(str)
        .tolist(),
        key=chromosome_sort_key,
    )

    sorted_intervals["chrom"] = pd.Categorical(
        sorted_intervals["chrom"],
        categories=ordered_chromosomes,
        ordered=True,
    )

    sorted_intervals.sort_values(
        by=["chrom", "start", "end"],
        kind="stable",
        inplace=True,
        ignore_index=True,
    )

    return sorted_intervals


def load_and_sort_all_peaks(
    sample_peak_files: Sequence[SamplePeakFile],
) -> pd.DataFrame:
    """
    Load all sample-specific peaks and sort them in genomic coordinate order.

    Exact duplicate intervals are retained because they represent support from
    potentially different samples.

    :param sample_peak_files: Peak BED metadata for all samples.
    :return: Combined peak table sorted by chromosome, start, and end.
    """
    peak_tables: list[pd.DataFrame] = []

    for position, sample_peak_file in enumerate(
        sample_peak_files,
        start=1,
    ):
        print(
            f"[{position}/{len(sample_peak_files)}] "
            f"Loading peaks from {sample_peak_file.sample_id}..."
        )

        peak_tables.append(
            read_peak_bed(sample_peak_file)
        )

    all_peaks = pd.concat(
        peak_tables,
        axis=0,
        ignore_index=True,
        copy=False,
    )

    return sort_genomic_intervals(all_peaks)


def intervals_overlapping_any_target(
    query_intervals: pd.DataFrame,
    target_intervals: pd.DataFrame,
) -> NDArray[np.bool_]:
    """
    Identify query intervals that overlap at least one target interval.

    Two BED intervals overlap when they share at least one genomic base:
    query_start < target_end and target_start < query_end.

    :param query_intervals: Intervals whose overlap status should be calculated.
    :param target_intervals: Intervals against which queries should be tested.
    :return: Boolean array aligned with the query interval rows.
    """
    overlap_mask = np.zeros(
        len(query_intervals),
        dtype=bool,
    )

    target_by_chromosome: dict[
        str,
        tuple[
            NDArray[np.int64],
            NDArray[np.int64],
        ],
    ] = {}

    for chromosome, chromosome_targets in target_intervals.groupby(
        "chrom",
        sort=False,
        observed=True,
    ):
        target_starts = chromosome_targets[
            "start"
        ].to_numpy(dtype=np.int64)

        target_ends = chromosome_targets[
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

    query_chromosomes = (
        query_intervals["chrom"]
        .astype(str)
        .to_numpy()
    )

    for chromosome in np.unique(query_chromosomes):
        target_data = target_by_chromosome.get(chromosome)

        if target_data is None:
            continue

        target_starts, target_prefix_max_ends = target_data

        query_positions = np.flatnonzero(
            query_chromosomes == chromosome
        )

        chromosome_queries = query_intervals.iloc[
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

        chromosome_overlap_mask = np.zeros(
            len(query_positions),
            dtype=bool,
        )

        valid_candidate_indices = (
            last_candidate_indices[has_candidate]
        )

        chromosome_overlap_mask[has_candidate] = (
            target_prefix_max_ends[
                valid_candidate_indices
            ]
            > query_starts[has_candidate]
        )

        overlap_mask[
            query_positions
        ] = chromosome_overlap_mask

    return overlap_mask


def remove_blacklisted_intervals(
    intervals: pd.DataFrame,
    blacklist: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Remove intervals that overlap any ENCODE blacklist region.

    The removed intervals are returned separately for auditing and summary
    calculations.

    :param intervals: Genomic intervals to filter.
    :param blacklist: ENCODE blacklist intervals.
    :return: Tuple containing retained intervals and removed intervals.
    """
    shared_chromosomes = (
        set(intervals["chrom"].astype(str))
        & set(blacklist["chrom"].astype(str))
    )

    if not shared_chromosomes:
        raise ValueError(
            "The peak and blacklist files have no chromosome names in common. "
            "Check whether one file uses names such as 'chr1' and the other "
            "uses names such as '1'."
        )

    blacklist_overlap_mask = intervals_overlapping_any_target(
        query_intervals=intervals,
        target_intervals=blacklist,
    )

    retained_intervals = (
        intervals.loc[~blacklist_overlap_mask]
        .copy()
        .reset_index(drop=True)
    )

    removed_intervals = (
        intervals.loc[blacklist_overlap_mask]
        .copy()
        .reset_index(drop=True)
    )

    return retained_intervals, removed_intervals


def merge_peaks(
    sorted_peaks: pd.DataFrame,
    merge_bookended_peaks: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merge sample-specific peaks into consensus intervals using Bioframe.

    Bioframe assigns every original peak to a transitive interval cluster.
    Peaks in the same cluster are represented by the common cluster_start and
    cluster_end coordinates. The cluster assignments are then summarized to
    calculate the number of original peaks and distinct contributing samples.

    When merge_bookended_peaks is True, directly adjacent BED intervals are
    merged using min_dist=0. Otherwise, min_dist=None requires intervals to
    overlap by at least one genomic base.

    :param sorted_peaks: Input peaks sorted by chromosome, start, and end.
    :param merge_bookended_peaks: Whether directly adjacent intervals are merged.
    :return: Tuple of (consensus peaks with interval length and sample-support
        statistics, clustered peaks DataFrame retaining the original per-peak
        rows with their assigned cluster coordinates).
    """
    if sorted_peaks.empty:
        raise ValueError("Cannot merge an empty peak table.")

    required_columns = {
        "chrom",
        "start",
        "end",
        "sample_index",
    }
    missing_columns = required_columns - set(sorted_peaks.columns)

    if missing_columns:
        raise ValueError(
            "The peak table is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    bioframe_input = sorted_peaks.copy()

    # Bioframe works most reliably with ordinary string chromosome names,
    # rather than pandas categorical values.
    bioframe_input["chrom"] = (
        bioframe_input["chrom"]
        .astype("string")
    )

    minimum_distance: int | None = (
        0 if merge_bookended_peaks else None
    )

    clustered_peaks = bf.cluster(
        bioframe_input,
        min_dist=minimum_distance,
        cols=("chrom", "start", "end"),
        return_input=True,
        return_cluster_ids=True,
        return_cluster_intervals=True,
    )

    expected_bioframe_columns = {
        "cluster",
        "cluster_start",
        "cluster_end",
    }
    missing_bioframe_columns = (
        expected_bioframe_columns
        - set(clustered_peaks.columns)
    )

    if missing_bioframe_columns:
        raise RuntimeError(
            "Bioframe did not return the expected clustering columns: "
            f"{sorted(missing_bioframe_columns)}"
        )

    consensus_peaks = (
        clustered_peaks.groupby(
            ["chrom", "cluster"],
            sort=False,
            observed=True,
        )
        .agg(
            start=("cluster_start", "first"),
            end=("cluster_end", "first"),
            n_source_peaks=("cluster", "size"),
            n_contributing_samples=(
                "sample_index",
                "nunique",
            ),
        )
        .reset_index()
    )

    consensus_peaks["start"] = (
        consensus_peaks["start"].astype("int64")
    )
    consensus_peaks["end"] = (
        consensus_peaks["end"].astype("int64")
    )

    consensus_peaks["length_bp"] = (
        consensus_peaks["end"]
        - consensus_peaks["start"]
    )

    consensus_peaks = consensus_peaks[
        [
            "chrom",
            "start",
            "end",
            "length_bp",
            "n_source_peaks",
            "n_contributing_samples",
        ]
    ]

    return sort_genomic_intervals(consensus_peaks), clustered_peaks


def create_reproducible_consensus_peaks(
    consensus_peaks: pd.DataFrame,
    minimum_sample_support: int,
) -> pd.DataFrame:
    """
    Retain consensus peaks supported by a minimum number of distinct samples.

    Support is based on n_contributing_samples, not on the number of original
    peak rows. Therefore, several overlapping peaks from one sample count as
    support from only one donor.

    :param consensus_peaks: Blacklist-filtered consensus peak table.
    :param minimum_sample_support: Minimum number of contributing samples.
    :return: Reproducible subset of the consensus peaks.
    """
    if minimum_sample_support < 1:
        raise ValueError(
            "minimum_sample_support must be at least 1."
        )

    reproducible_peaks = consensus_peaks.loc[
        consensus_peaks["n_contributing_samples"]
        >= minimum_sample_support
    ].copy()

    reproducible_peaks.reset_index(
        drop=True,
        inplace=True,
    )

    if reproducible_peaks.empty:
        raise ValueError(
            "No reproducible consensus peaks remained after applying "
            f"a minimum support of {minimum_sample_support} samples."
        )

    return reproducible_peaks


def create_peak_to_consensus_map(
    clustered_peaks: pd.DataFrame,
    sample_peak_files: Sequence[SamplePeakFile],
    reproducible_consensus_peaks: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build a mapping from each original sample-specific peak to the reproducible
    consensus peak it was merged into.

    Only input peaks whose cluster survived the reproducibility filter are
    included. Peaks whose merged cluster was removed by the blacklist safety
    check or did not reach the minimum sample support are excluded.

    :param clustered_peaks: Pre-aggregation output of bf.cluster, containing
        the original chrom/start/end/sample_index columns plus
        cluster/cluster_start/cluster_end columns.
    :param sample_peak_files: Ordered list of sample metadata used to recover
        sample identifiers from integer indices.
    :param reproducible_consensus_peaks: Filtered consensus peak table used to
        select which clusters appear in the final reproducible set.
    :return: DataFrame with columns sample_id, original_chrom, original_start,
        original_end, consensus_chrom, consensus_start, consensus_end.
    """
    sample_id_by_index: dict[int, str] = {
        spf.sample_index: spf.sample_id
        for spf in sample_peak_files
    }

    # Build a lookup keyed by the cluster's merged coordinates so we can
    # filter clustered_peaks to only reproducible clusters.
    reproducible_keys = set(
        zip(
            reproducible_consensus_peaks["chrom"].astype(str),
            reproducible_consensus_peaks["start"],
            reproducible_consensus_peaks["end"],
        )
    )

    cluster_key = list(
        zip(
            clustered_peaks["chrom"].astype(str),
            clustered_peaks["cluster_start"],
            clustered_peaks["cluster_end"],
        )
    )

    reproducible_mask = np.array(
        [k in reproducible_keys for k in cluster_key],
        dtype=bool,
    )

    filtered = clustered_peaks.loc[reproducible_mask].copy()

    filtered["sample_id"] = (
        filtered["sample_index"]
        .map(sample_id_by_index)
        .astype("string")
    )

    filtered = filtered.rename(
        columns={
            "chrom": "original_chrom",
            "start": "original_start",
            "end": "original_end",
            "cluster_start": "consensus_start",
            "cluster_end": "consensus_end",
        }
    )

    filtered["consensus_chrom"] = filtered["original_chrom"]

    result = filtered[
        [
            "sample_id",
            "original_chrom",
            "original_start",
            "original_end",
            "consensus_chrom",
            "consensus_start",
            "consensus_end",
        ]
    ].reset_index(drop=True)

    return result


def create_peak_set_statistics(
    peaks: pd.DataFrame,
    prefix: str,
) -> dict[str, int | float]:
    """
    Calculate descriptive statistics for one consensus peak set.

    :param peaks: Consensus peak table containing a length_bp column.
    :param prefix: Prefix added to every returned statistic name.
    :return: Dictionary of peak count and interval-length statistics.
    """
    if peaks.empty:
        return {
            f"{prefix}_n_peaks": 0,
            f"{prefix}_total_bp": 0,
            f"{prefix}_peak_length_min": np.nan,
            f"{prefix}_peak_length_mean": np.nan,
            f"{prefix}_peak_length_median": np.nan,
            f"{prefix}_peak_length_p95": np.nan,
            f"{prefix}_peak_length_max": np.nan,
        }

    return {
        f"{prefix}_n_peaks": int(len(peaks)),
        f"{prefix}_total_bp": int(
            peaks["length_bp"].sum()
        ),
        f"{prefix}_peak_length_min": int(
            peaks["length_bp"].min()
        ),
        f"{prefix}_peak_length_mean": float(
            peaks["length_bp"].mean()
        ),
        f"{prefix}_peak_length_median": float(
            peaks["length_bp"].median()
        ),
        f"{prefix}_peak_length_p95": float(
            peaks["length_bp"].quantile(0.95)
        ),
        f"{prefix}_peak_length_max": int(
            peaks["length_bp"].max()
        ),
    }


def create_consensus_summary(
    all_input_peaks: pd.DataFrame,
    blacklisted_input_peaks: pd.DataFrame,
    merged_consensus_peaks: pd.DataFrame,
    blacklisted_consensus_peaks: pd.DataFrame,
    consensus_peaks: pd.DataFrame,
    reproducible_consensus_peaks: pd.DataFrame,
    blacklist: pd.DataFrame,
    n_samples: int,
    merge_bookended_peaks: bool,
    minimum_sample_support: int,
) -> pd.DataFrame:
    """
    Create a one-row summary of blacklist filtering and consensus construction.

    :param all_input_peaks: All original sample-specific peaks.
    :param blacklisted_input_peaks: Input peaks removed by blacklist overlap.
    :param merged_consensus_peaks: Consensus peaks before the final safety filter.
    :param blacklisted_consensus_peaks: Merged peaks removed by the safety filter.
    :param consensus_peaks: Final blacklist-filtered full consensus set.
    :param reproducible_consensus_peaks: Consensus peaks passing support filtering.
    :param blacklist: ENCODE blacklist intervals.
    :param n_samples: Number of samples contributing peaks.
    :param merge_bookended_peaks: Whether directly adjacent peaks were merged.
    :param minimum_sample_support: Support threshold for reproducible peaks.
    :return: One-row DataFrame containing construction statistics.
    """
    n_input_peaks = len(all_input_peaks)
    n_blacklisted_input_peaks = len(
        blacklisted_input_peaks
    )
    n_consensus_peaks = len(consensus_peaks)
    n_reproducible_peaks = len(
        reproducible_consensus_peaks
    )

    summary_row: dict[str, object] = {
        "n_samples": n_samples,
        "n_blacklist_intervals": len(blacklist),
        "n_input_peaks_before_blacklist": n_input_peaks,
        "n_exact_unique_input_intervals": int(
            all_input_peaks[
                ["chrom", "start", "end"]
            ]
            .drop_duplicates()
            .shape[0]
        ),
        "n_input_peaks_removed_by_blacklist": (
            n_blacklisted_input_peaks
        ),
        "fraction_input_peaks_removed_by_blacklist": (
            n_blacklisted_input_peaks / n_input_peaks
        ),
        "n_input_peaks_after_blacklist": (
            n_input_peaks - n_blacklisted_input_peaks
        ),
        "n_merged_consensus_peaks_before_safety_filter": (
            len(merged_consensus_peaks)
        ),
        "n_consensus_peaks_removed_by_blacklist_safety_filter": (
            len(blacklisted_consensus_peaks)
        ),
        "merge_bookended_peaks": merge_bookended_peaks,
        "minimum_reproducible_sample_support": (
            minimum_sample_support
        ),
        "n_reproducible_consensus_peaks": (
            n_reproducible_peaks
        ),
        "fraction_consensus_peaks_reproducible": (
            n_reproducible_peaks / n_consensus_peaks
        ),
        "median_source_peaks_per_consensus_peak": float(
            consensus_peaks[
                "n_source_peaks"
            ].median()
        ),
        "median_contributing_samples_per_consensus_peak": float(
            consensus_peaks[
                "n_contributing_samples"
            ].median()
        ),
        "fraction_supported_by_all_samples": float(
            (
                consensus_peaks[
                    "n_contributing_samples"
                ]
                == n_samples
            ).mean()
        ),
    }

    summary_row.update(
        create_peak_set_statistics(
            peaks=consensus_peaks,
            prefix="full_consensus",
        )
    )

    summary_row.update(
        create_peak_set_statistics(
            peaks=reproducible_consensus_peaks,
            prefix="reproducible_consensus",
        )
    )

    return pd.DataFrame([summary_row])


def save_bed3(
    peaks: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Save the chrom, start, and end columns as a headerless BED3 file.

    :param peaks: Peak table to save.
    :param output_path: Destination BED file.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    peaks[["chrom", "start", "end"]].to_csv(
        output_path,
        sep="\t",
        header=False,
        index=False,
    )


def save_consensus_outputs(
    consensus_peaks: pd.DataFrame,
    reproducible_consensus_peaks: pd.DataFrame,
    summary: pd.DataFrame,
    peak_to_consensus_map: pd.DataFrame,
) -> None:
    """
    Save the full consensus, reproducible consensus, annotations, summary,
    and the per-input-peak-to-reproducible-consensus mapping.

    :param consensus_peaks: Final blacklist-filtered full consensus peak table.
    :param reproducible_consensus_peaks: Reproducible consensus peak subset.
    :param summary: One-row construction summary.
    :param peak_to_consensus_map: Mapping from original sample peaks to their
        reproducible consensus peak.
    """
    output_paths = (
        CONSENSUS_PEAKS_OUTPUT_PATH,
        CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
        CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH,
        REPRODUCIBLE_CONSENSUS_PEAKS_OUTPUT_PATH,
        REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
        PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH,
    )

    for output_path in output_paths:
        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    save_bed3(
        peaks=consensus_peaks,
        output_path=CONSENSUS_PEAKS_OUTPUT_PATH,
    )

    consensus_peaks.to_csv(
        CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
        index=False,
    )

    save_bed3(
        peaks=reproducible_consensus_peaks,
        output_path=(
            REPRODUCIBLE_CONSENSUS_PEAKS_OUTPUT_PATH
        ),
    )

    reproducible_consensus_peaks.to_csv(
        REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
        index=False,
    )

    summary.to_csv(
        CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH,
        index=False,
    )

    peak_to_consensus_map.to_csv(
        PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH,
        index=False,
    )


def print_consensus_summary(
    summary: pd.DataFrame,
) -> None:
    """
    Print the most important blacklist and consensus-construction statistics.

    :param summary: One-row consensus construction summary.
    """
    row = summary.iloc[0]

    print("\nConsensus peak construction completed.")

    print(
        "Input peaks: "
        f"{int(row['n_input_peaks_before_blacklist']):,}"
    )

    print(
        "Input peaks removed by blacklist: "
        f"{int(row['n_input_peaks_removed_by_blacklist']):,} "
        f"({100.0 * float(row['fraction_input_peaks_removed_by_blacklist']):.2f}%)"
    )

    print(
        "Full blacklist-filtered consensus peaks: "
        f"{int(row['full_consensus_n_peaks']):,}"
    )

    print(
        "Reproducible consensus peaks: "
        f"{int(row['n_reproducible_consensus_peaks']):,} "
        f"({100.0 * float(row['fraction_consensus_peaks_reproducible']):.2f}%)"
    )

    print(
        "Minimum reproducible support: "
        f"{int(row['minimum_reproducible_sample_support'])} samples"
    )

    print(
        "Median full-consensus length: "
        f"{float(row['full_consensus_peak_length_median']):,.1f} bp"
    )

    print(
        "Median reproducible-consensus length: "
        f"{float(row['reproducible_consensus_peak_length_median']):,.1f} bp"
    )


def main() -> None:
    """
    Create blacklist-filtered full and reproducible consensus peak sets.
    """
    sample_peak_files = collect_sample_peak_files(
        ATAC_SEQ_DIRS
    )

    all_input_peaks = load_and_sort_all_peaks(
        sample_peak_files
    )

    print(
        f"\nReading ENCODE blacklist from "
        f"{BLACKLIST_BED_PATH}..."
    )

    blacklist = read_blacklist(
        BLACKLIST_BED_PATH
    )

    print(
        f"Removing input peaks overlapping "
        f"{len(blacklist):,} blacklist intervals..."
    )

    (
        non_blacklisted_input_peaks,
        blacklisted_input_peaks,
    ) = remove_blacklisted_intervals(
        intervals=all_input_peaks,
        blacklist=blacklist,
    )

    print(
        f"Merging "
        f"{len(non_blacklisted_input_peaks):,} "
        "non-blacklisted input peaks..."
    )

    merged_consensus_peaks, clustered_peaks = merge_peaks(
        sorted_peaks=non_blacklisted_input_peaks,
        merge_bookended_peaks=MERGE_BOOKENDED_PEAKS,
    )

    # This second filtering step is a safety check. It should usually remove
    # zero intervals because blacklisted input peaks were already excluded.
    (
        consensus_peaks,
        blacklisted_consensus_peaks,
    ) = remove_blacklisted_intervals(
        intervals=merged_consensus_peaks,
        blacklist=blacklist,
    )

    reproducible_consensus_peaks = (
        create_reproducible_consensus_peaks(
            consensus_peaks=consensus_peaks,
            minimum_sample_support=(
                MIN_REPRODUCIBLE_SAMPLE_SUPPORT
            ),
        )
    )

    peak_to_consensus_map = create_peak_to_consensus_map(
        clustered_peaks=clustered_peaks,
        sample_peak_files=sample_peak_files,
        reproducible_consensus_peaks=reproducible_consensus_peaks,
    )

    summary = create_consensus_summary(
        all_input_peaks=all_input_peaks,
        blacklisted_input_peaks=blacklisted_input_peaks,
        merged_consensus_peaks=merged_consensus_peaks,
        blacklisted_consensus_peaks=(
            blacklisted_consensus_peaks
        ),
        consensus_peaks=consensus_peaks,
        reproducible_consensus_peaks=(
            reproducible_consensus_peaks
        ),
        blacklist=blacklist,
        n_samples=len(sample_peak_files),
        merge_bookended_peaks=(
            MERGE_BOOKENDED_PEAKS
        ),
        minimum_sample_support=(
            MIN_REPRODUCIBLE_SAMPLE_SUPPORT
        ),
    )

    save_consensus_outputs(
        consensus_peaks=consensus_peaks,
        reproducible_consensus_peaks=(
            reproducible_consensus_peaks
        ),
        summary=summary,
        peak_to_consensus_map=peak_to_consensus_map,
    )

    print_consensus_summary(summary)

    print("\nSaved outputs:")
    print(
        f"  Full consensus BED: "
        f"{CONSENSUS_PEAKS_OUTPUT_PATH}"
    )
    print(
        f"  Full annotations: "
        f"{CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH}"
    )
    print(
        f"  Reproducible consensus BED: "
        f"{REPRODUCIBLE_CONSENSUS_PEAKS_OUTPUT_PATH}"
    )
    print(
        f"  Reproducible annotations: "
        f"{REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH}"
    )
    print(
        f"  Construction summary: "
        f"{CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH}"
    )
    print(
        f"  Peak-to-consensus map: "
        f"{PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH}"
    )


if __name__ == "__main__":
    main()