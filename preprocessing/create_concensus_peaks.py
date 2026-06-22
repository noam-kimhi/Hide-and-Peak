from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from constants import (
    ATAC_SEQ_DIRS,
    CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
    CONSENSUS_PEAKS_OUTPUT_PATH,
    CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH,
    MERGE_BOOKENDED_PEAKS,
    PEAKS_SUFFIX,
)


@dataclass(frozen=True)
class SamplePeakFile:
    """Metadata for one sample-specific peak BED file."""

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
                peaks_path=find_single_file(sample_dir, PEAKS_SUFFIX),
                sample_index=sample_index,
            )
        )

    sample_ids = [sample.sample_id for sample in sample_peak_files]

    if len(sample_ids) != len(set(sample_ids)):
        duplicated_sample_ids = sorted(
            sample_id
            for sample_id in set(sample_ids)
            if sample_ids.count(sample_id) > 1
        )
        raise ValueError(
            f"Duplicate sample identifiers were found: {duplicated_sample_ids}"
        )

    if not sample_peak_files:
        raise ValueError("No sample directories were provided.")

    return sample_peak_files


def read_peak_bed(sample_peak_file: SamplePeakFile) -> pd.DataFrame:
    """
    Read and validate the BED3 peak intervals belonging to one sample.

    BED coordinates are interpreted as zero-based, half-open intervals:
    [start, end).

    :param sample_peak_file: Metadata and BED path for one sample.
    :return: Validated peak table with chrom, start, end, and sample_index columns.
    """
    peaks = pd.read_csv(
        sample_peak_file.peaks_path,
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
        raise ValueError(
            f"Peak BED file is empty: {sample_peak_file.peaks_path}"
        )

    if peaks.isna().any(axis=None):
        raise ValueError(
            f"Peak BED file contains missing values: "
            f"{sample_peak_file.peaks_path}"
        )

    negative_start_mask = peaks["start"] < 0
    invalid_interval_mask = peaks["end"] <= peaks["start"]

    if negative_start_mask.any():
        raise ValueError(
            f"{sample_peak_file.peaks_path} contains "
            f"{int(negative_start_mask.sum())} negative start coordinates."
        )

    if invalid_interval_mask.any():
        raise ValueError(
            f"{sample_peak_file.peaks_path} contains "
            f"{int(invalid_interval_mask.sum())} intervals with end <= start."
        )

    peaks["sample_index"] = np.int16(sample_peak_file.sample_index)

    return peaks


def chromosome_sort_key(chromosome: str) -> tuple[int, int, str]:
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
        return 0, special_chromosome_order[normalized_upper], ""

    return 1, 0, chromosome


def load_and_sort_all_peaks(
    sample_peak_files: Sequence[SamplePeakFile],
) -> pd.DataFrame:
    """
    Load all sample-specific peaks and sort them in genomic coordinate order.

    Exact duplicate intervals are retained because the number of source peaks
    and contributing samples is recorded during merging.

    :param sample_peak_files: Peak BED metadata for all samples.
    :return: Combined peak table sorted by chromosome, start, and end.
    """
    peak_tables: list[pd.DataFrame] = []

    for position, sample_peak_file in enumerate(sample_peak_files, start=1):
        print(
            f"[{position}/{len(sample_peak_files)}] "
            f"Loading peaks from {sample_peak_file.sample_id}..."
        )
        peak_tables.append(read_peak_bed(sample_peak_file))

    all_peaks = pd.concat(
        peak_tables,
        axis=0,
        ignore_index=True,
        copy=False,
    )

    ordered_chromosomes = sorted(
        all_peaks["chrom"].drop_duplicates().astype(str).tolist(),
        key=chromosome_sort_key,
    )

    all_peaks["chrom"] = pd.Categorical(
        all_peaks["chrom"],
        categories=ordered_chromosomes,
        ordered=True,
    )

    all_peaks.sort_values(
        by=["chrom", "start", "end"],
        kind="stable",
        inplace=True,
        ignore_index=True,
    )

    return all_peaks


def intervals_should_merge(
    current_chrom: str,
    current_end: int,
    next_chrom: str,
    next_start: int,
    merge_bookended_peaks: bool,
) -> bool:
    """
    Determine whether a sorted interval should be merged into the current one.

    When book-ended merging is enabled, intervals such as [100, 200) and
    [200, 300) are merged even though they do not share a genomic base.

    :param current_chrom: Chromosome of the currently accumulated interval.
    :param current_end: End coordinate of the currently accumulated interval.
    :param next_chrom: Chromosome of the next interval.
    :param next_start: Start coordinate of the next interval.
    :param merge_bookended_peaks: Whether directly adjacent intervals are merged.
    :return: True when the intervals should belong to the same consensus peak.
    """
    if current_chrom != next_chrom:
        return False

    if merge_bookended_peaks:
        return next_start <= current_end

    return next_start < current_end


def append_consensus_interval(
    output_rows: list[dict[str, object]],
    chrom: str,
    start: int,
    end: int,
    source_peak_count: int,
    contributing_sample_indices: set[int],
) -> None:
    """
    Append one completed consensus interval to the output collection.

    :param output_rows: Mutable collection receiving consensus interval records.
    :param chrom: Chromosome of the consensus interval.
    :param start: Start coordinate of the consensus interval.
    :param end: End coordinate of the consensus interval.
    :param source_peak_count: Number of original peaks merged into the interval.
    :param contributing_sample_indices: Samples contributing at least one peak.
    """
    output_rows.append(
        {
            "chrom": chrom,
            "start": start,
            "end": end,
            "length_bp": end - start,
            "n_source_peaks": source_peak_count,
            "n_contributing_samples": len(contributing_sample_indices),
        }
    )


def merge_peaks(
    sorted_peaks: pd.DataFrame,
    merge_bookended_peaks: bool,
) -> pd.DataFrame:
    """
    Merge overlapping sample-specific peaks into shared consensus intervals.

    Merging is transitive. For example, if peak A overlaps peak B and peak B
    overlaps peak C, all three peaks become one consensus interval even when
    peak A and peak C do not directly overlap.

    :param sorted_peaks: All input peaks sorted by chromosome, start, and end.
    :param merge_bookended_peaks: Whether directly adjacent intervals are merged.
    :return: Consensus peaks with interval length and source-support statistics.
    """
    if sorted_peaks.empty:
        raise ValueError("Cannot merge an empty peak table.")

    rows = sorted_peaks.itertuples(index=False)

    first_peak = next(rows)
    current_chrom = str(first_peak.chrom)
    current_start = int(first_peak.start)
    current_end = int(first_peak.end)
    current_source_peak_count = 1
    current_sample_indices = {int(first_peak.sample_index)}

    consensus_rows: list[dict[str, object]] = []

    for peak in rows:
        next_chrom = str(peak.chrom)
        next_start = int(peak.start)
        next_end = int(peak.end)
        next_sample_index = int(peak.sample_index)

        if intervals_should_merge(
            current_chrom=current_chrom,
            current_end=current_end,
            next_chrom=next_chrom,
            next_start=next_start,
            merge_bookended_peaks=merge_bookended_peaks,
        ):
            current_end = max(current_end, next_end)
            current_source_peak_count += 1
            current_sample_indices.add(next_sample_index)
            continue

        append_consensus_interval(
            output_rows=consensus_rows,
            chrom=current_chrom,
            start=current_start,
            end=current_end,
            source_peak_count=current_source_peak_count,
            contributing_sample_indices=current_sample_indices,
        )

        current_chrom = next_chrom
        current_start = next_start
        current_end = next_end
        current_source_peak_count = 1
        current_sample_indices = {next_sample_index}

    append_consensus_interval(
        output_rows=consensus_rows,
        chrom=current_chrom,
        start=current_start,
        end=current_end,
        source_peak_count=current_source_peak_count,
        contributing_sample_indices=current_sample_indices,
    )

    return pd.DataFrame(consensus_rows)


def create_consensus_summary(
    all_peaks: pd.DataFrame,
    consensus_peaks: pd.DataFrame,
    n_samples: int,
    merge_bookended_peaks: bool,
) -> pd.DataFrame:
    """
    Create a one-row summary of the consensus-peak construction.

    :param all_peaks: Combined table of all original sample-specific peaks.
    :param consensus_peaks: Merged consensus peak table.
    :param n_samples: Number of samples contributing peaks.
    :param merge_bookended_peaks: Whether directly adjacent peaks were merged.
    :return: One-row DataFrame containing major construction statistics.
    """
    n_input_peaks = len(all_peaks)
    n_consensus_peaks = len(consensus_peaks)

    return pd.DataFrame(
        [
            {
                "n_samples": n_samples,
                "n_input_peaks": n_input_peaks,
                "n_exact_unique_input_intervals": int(
                    all_peaks[["chrom", "start", "end"]]
                    .drop_duplicates()
                    .shape[0]
                ),
                "n_consensus_peaks": n_consensus_peaks,
                "reduction_fraction": (
                    1.0 - n_consensus_peaks / n_input_peaks
                ),
                "merge_bookended_peaks": merge_bookended_peaks,
                "consensus_total_bp": int(
                    consensus_peaks["length_bp"].sum()
                ),
                "consensus_peak_length_min": int(
                    consensus_peaks["length_bp"].min()
                ),
                "consensus_peak_length_mean": float(
                    consensus_peaks["length_bp"].mean()
                ),
                "consensus_peak_length_median": float(
                    consensus_peaks["length_bp"].median()
                ),
                "consensus_peak_length_p95": float(
                    consensus_peaks["length_bp"].quantile(0.95)
                ),
                "consensus_peak_length_max": int(
                    consensus_peaks["length_bp"].max()
                ),
                "median_source_peaks_per_consensus_peak": float(
                    consensus_peaks["n_source_peaks"].median()
                ),
                "median_contributing_samples_per_consensus_peak": float(
                    consensus_peaks["n_contributing_samples"].median()
                ),
                "fraction_supported_by_at_least_2_samples": float(
                    (
                        consensus_peaks["n_contributing_samples"] >= 2
                    ).mean()
                ),
                "fraction_supported_by_all_samples": float(
                    (
                        consensus_peaks["n_contributing_samples"] == n_samples
                    ).mean()
                ),
            }
        ]
    )


def save_consensus_outputs(
    consensus_peaks: pd.DataFrame,
    summary: pd.DataFrame,
    consensus_bed_path: Path,
    annotation_output_path: Path,
    summary_output_path: Path,
) -> None:
    """
    Save the BED3 consensus set, its annotations, and a construction summary.

    The BED output contains only chrom, start, and end so it can be used
    directly by downstream interval-counting tools.

    :param consensus_peaks: Merged consensus peak table.
    :param summary: One-row construction summary.
    :param consensus_bed_path: Destination path for the BED3 consensus set.
    :param annotation_output_path: Destination CSV for per-peak annotations.
    :param summary_output_path: Destination CSV for global summary statistics.
    """
    for output_path in (
        consensus_bed_path,
        annotation_output_path,
        summary_output_path,
    ):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    consensus_peaks[["chrom", "start", "end"]].to_csv(
        consensus_bed_path,
        sep="\t",
        header=False,
        index=False,
    )

    consensus_peaks.to_csv(
        annotation_output_path,
        index=False,
    )

    summary.to_csv(
        summary_output_path,
        index=False,
    )


def print_consensus_summary(summary: pd.DataFrame) -> None:
    """
    Print the most important consensus-peak construction statistics.

    :param summary: One-row consensus construction summary.
    """
    row = summary.iloc[0]

    print("\nConsensus peak construction completed.")
    print(f"Input samples: {int(row['n_samples']):,}")
    print(f"Input peaks: {int(row['n_input_peaks']):,}")
    print(f"Consensus peaks: {int(row['n_consensus_peaks']):,}")
    print(
        "Reduction: "
        f"{100.0 * float(row['reduction_fraction']):.2f}%"
    )
    print(
        "Median consensus length: "
        f"{float(row['consensus_peak_length_median']):,.1f} bp"
    )
    print(
        "95th percentile consensus length: "
        f"{float(row['consensus_peak_length_p95']):,.1f} bp"
    )
    print(
        "Fraction supported by at least two samples: "
        f"{100.0 * float(row['fraction_supported_by_at_least_2_samples']):.2f}%"
    )


def main() -> None:
    """Create and save a consensus peak set from all sample-specific BED files."""
    sample_peak_files = collect_sample_peak_files(ATAC_SEQ_DIRS)
    all_peaks = load_and_sort_all_peaks(sample_peak_files)

    print(f"\nMerging {len(all_peaks):,} input peaks...")

    consensus_peaks = merge_peaks(
        sorted_peaks=all_peaks,
        merge_bookended_peaks=MERGE_BOOKENDED_PEAKS,
    )

    summary = create_consensus_summary(
        all_peaks=all_peaks,
        consensus_peaks=consensus_peaks,
        n_samples=len(sample_peak_files),
        merge_bookended_peaks=MERGE_BOOKENDED_PEAKS,
    )

    save_consensus_outputs(
        consensus_peaks=consensus_peaks,
        summary=summary,
        consensus_bed_path=CONSENSUS_PEAKS_OUTPUT_PATH,
        annotation_output_path=CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
        summary_output_path=CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH,
    )

    print_consensus_summary(summary)

    print("\nSaved outputs:")
    print(f"  Consensus BED: {CONSENSUS_PEAKS_OUTPUT_PATH}")
    print(f"  Peak annotations: {CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH}")
    print(f"  Construction summary: {CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH}")


if __name__ == "__main__":
    main()
