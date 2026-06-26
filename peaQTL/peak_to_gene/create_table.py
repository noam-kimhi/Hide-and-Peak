"""Build pseudo-bulk count tables for each cell type.

For each cell type, this script:

1. Loads the reproducible consensus peak set (the shared row space).
2. Loads the peak-to-consensus mapping produced by create_consensus_peaks_bioframe.py.
3. For every sample and cell type, loads the binarized sub-matrix, sums
   accessibility across all cells (columns) to obtain a per-original-peak
   open-cell count, then re-indexes those counts into the consensus peak space.
4. Stacks the 12 per-sample vectors into a (n_consensus_peaks × 12) DataFrame
   and saves it as a CSV.

Samples that are entirely missing (e.g. GSM8619372_Normal_rep4, which had no
matched barcodes during splitting) produce a zero column.  Cell types that are
absent for a particular sample also produce a zero column.

Output per cell type:
    data/derived/pseudobulk/{cell_type}_pseudobulk.csv

Rows  : consensus peak IDs in 'chrom:start-end' format, same order as
        reproducible_consensus_peaks_annotations.csv.
Columns : Normal_rep1 … Normal_rep6, MASH_rep1 … MASH_rep6.
Values  : integer sum of binarized open-chromatin events from all valid cells
          belonging to that replicate.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.io import mmread

from constants import (
    BINARIZED_CELL_TYPE_MATRICES_DIR,
    CELL_TYPE_STANDARDIZATION,
    MATRIX_SUFFIX,
    PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH,
    PEAKS_SUFFIX,
    PSEUDOBULK_OUTPUT_DIR,
    PSEUDOBULK_SUMMARY_PATH,
    REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH,
)

# Ordered sample names for output columns.
# Normal replicates appear first, then MASH replicates (following pipeline.md).
ORDERED_SAMPLE_NAMES: list[str] = [
    "GSM8619369_Normal_rep1",
    "GSM8619370_Normal_rep2",
    "GSM8619371_Normal_rep3",
    "GSM8619372_Normal_rep4",   # missing — produces a zero column
    "GSM8619373_Normal_rep5",
    "GSM8619374_Normal_rep6",
    "GSM8619363_MASH_rep1",
    "GSM8619364_MASH_rep2",
    "GSM8619365_MASH_rep3",
    "GSM8619366_MASH_rep4",
    "GSM8619367_MASH_rep5",
    "GSM8619368_MASH_rep6",
]


def parse_column_name(sample_name: str) -> str:
    """
    Derive a readable column label from a sample directory name.

    Strips the leading GSM accession prefix and returns the remainder.

    :param sample_name: Sample directory name, e.g. 'GSM8619363_MASH_rep1'.
    :return: Column label, e.g. 'MASH_rep1'.
    """
    _, _, remainder = sample_name.partition("_")
    return remainder


def load_consensus_peaks(
    annotations_path: Path,
) -> tuple[list[str], dict[tuple[str, int, int], int]]:
    """
    Load the reproducible consensus peak set and build a coordinate lookup.

    :param annotations_path: Path to reproducible_consensus_peaks_annotations.csv.
    :return: Tuple of (ordered peak-id strings in 'chrom:start-end' format,
        dict mapping (chrom, start, end) to the peak's integer row index).
    """
    peaks = pd.read_csv(
        annotations_path,
        usecols=["chrom", "start", "end"],
        dtype={"chrom": str, "start": "int64", "end": "int64"},
    )

    peak_ids: list[str] = [
        f"{row.chrom}:{row.start}-{row.end}"
        for row in peaks.itertuples(index=False)
    ]

    lookup: dict[tuple[str, int, int], int] = {
        (str(row.chrom), int(row.start), int(row.end)): idx
        for idx, row in enumerate(peaks.itertuples(index=False))
    }

    return peak_ids, lookup


def build_sample_peak_lookup(
    sample_mapping: pd.DataFrame,
    consensus_lookup: dict[tuple[str, int, int], int],
) -> dict[tuple[str, int, int], int]:
    """
    Build a mapping from original peak coordinates to consensus row indices
    for one sample.

    Only original peaks whose merged cluster appears in the reproducible
    consensus set are included; blacklisted or non-reproducible peaks are
    silently omitted.

    :param sample_mapping: Rows from the full mapping table for one sample.
    :param consensus_lookup: Coordinate-keyed dict for the reproducible
        consensus peaks.
    :return: Dict mapping (original_chrom, original_start, original_end) to
        the consensus peak's integer row index.
    """
    result: dict[tuple[str, int, int], int] = {}

    for row in sample_mapping.itertuples(index=False):
        consensus_idx = consensus_lookup.get(
            (
                str(row.consensus_chrom),
                int(row.consensus_start),
                int(row.consensus_end),
            )
        )

        if consensus_idx is not None:
            result[
                (
                    str(row.original_chrom),
                    int(row.original_start),
                    int(row.original_end),
                )
            ] = consensus_idx

    return result


def find_single_file(directory: Path, suffix: str) -> Path:
    """Return the unique file in *directory* whose name ends with *suffix*."""
    matches = sorted(directory.glob(f"*{suffix}"))

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one file ending with {suffix!r} in "
            f"{directory}, but found {len(matches)}."
        )

    return matches[0]


def read_peaks_bed(path: Path) -> list[tuple[str, int, int]]:
    """
    Read a compressed BED3 file and return an ordered list of coordinates.

    The list is aligned with the row order of the accompanying matrix file.

    :param path: Path to a gzip-compressed BED file.
    :return: List of (chrom, start, end) tuples, one per matrix row.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": "int64", "end": "int64"},
        compression="gzip",
    )

    return list(
        zip(df["chrom"], df["start"].astype(int), df["end"].astype(int))
    )


def aggregate_sample_cell_type(
    sample_name: str,
    cell_type: str,
    sample_peak_lookup: dict[tuple[str, int, int], int],
    n_consensus_peaks: int,
) -> np.ndarray:
    """
    Aggregate one sample's binarized cell-type matrix into a pseudo-bulk vector.

    Sums the binarized accessibility values across all cells for each original
    peak, then re-indexes those sums into the reproducible consensus peak space.
    Original peaks not present in *sample_peak_lookup* are skipped.

    Returns a zero vector if the sub-matrix directory does not exist (e.g. the
    sample was skipped during matrix splitting or has no cells of this type).

    :param sample_name: Sample directory name.
    :param cell_type: Standardised cell-type label.
    :param sample_peak_lookup: Maps (chrom, start, end) → consensus row index.
    :param n_consensus_peaks: Length of the output vector.
    :return: Integer vector of length *n_consensus_peaks*.
    """
    cell_type_dir = (
        BINARIZED_CELL_TYPE_MATRICES_DIR / sample_name / cell_type
    )

    result = np.zeros(n_consensus_peaks, dtype=np.int64)

    if not cell_type_dir.is_dir():
        return result

    matrix_path = find_single_file(cell_type_dir, MATRIX_SUFFIX)
    peaks_path = find_single_file(cell_type_dir, PEAKS_SUFFIX)

    original_peaks = read_peaks_bed(peaks_path)

    matrix = mmread(matrix_path).tocsr()
    per_peak_sums = (
        np.asarray(matrix.sum(axis=1)).flatten().astype(np.int64)
    )

    for row_i, coords in enumerate(original_peaks):
        consensus_idx = sample_peak_lookup.get(coords)

        if consensus_idx is not None:
            result[consensus_idx] += per_peak_sums[row_i]

    return result


def create_pseudobulk_table(
    cell_type: str,
    mapping_by_sample: dict[str, pd.DataFrame],
    peak_ids: list[str],
    consensus_lookup: dict[tuple[str, int, int], int],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Build the pseudo-bulk count table for one cell type across all samples.

    :param cell_type: Standardised cell-type label.
    :param mapping_by_sample: Pre-grouped mapping table keyed by sample name.
    :param peak_ids: Ordered list of consensus peak IDs (row labels).
    :param consensus_lookup: Coordinate-keyed dict for the consensus peaks.
    :return: Tuple of (DataFrame with peak IDs as index and sample columns,
        summary dict with basic statistics).
    """
    n_consensus_peaks = len(peak_ids)
    columns: dict[str, np.ndarray] = {}
    samples_with_data = 0

    for sample_name in ORDERED_SAMPLE_NAMES:
        col_name = parse_column_name(sample_name)

        sample_mapping = mapping_by_sample.get(
            sample_name, pd.DataFrame()
        )

        sample_peak_lookup = build_sample_peak_lookup(
            sample_mapping=sample_mapping,
            consensus_lookup=consensus_lookup,
        )

        vec = aggregate_sample_cell_type(
            sample_name=sample_name,
            cell_type=cell_type,
            sample_peak_lookup=sample_peak_lookup,
            n_consensus_peaks=n_consensus_peaks,
        )

        if vec.any():
            samples_with_data += 1

        columns[col_name] = vec

        print(
            f"    {col_name}: "
            f"{int(vec.sum()):,} total open-chromatin events, "
            f"{int((vec > 0).sum()):,} covered consensus peaks"
        )

    table = pd.DataFrame(
        columns,
        index=pd.Index(peak_ids, name="peak_id"),
    )

    summary: dict[str, Any] = {
        "cell_type": cell_type,
        "n_consensus_peaks": n_consensus_peaks,
        "n_samples": len(ORDERED_SAMPLE_NAMES),
        "n_samples_with_data": samples_with_data,
        "total_nonzero_entries": int((table > 0).values.sum()),
        "total_open_chromatin_events": int(table.values.sum()),
    }

    return table, summary


def main() -> None:
    """
    Create one pseudo-bulk count table per cell type.
    """
    PSEUDOBULK_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(
        f"Loading consensus peaks from "
        f"{REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH.name}..."
    )

    peak_ids, consensus_lookup = load_consensus_peaks(
        REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH
    )

    print(f"  {len(peak_ids):,} reproducible consensus peaks.")

    print(
        f"Loading peak-to-consensus map from "
        f"{PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH.name}..."
    )

    mapping_df = pd.read_csv(
        PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH,
        dtype={
            "sample_id": str,
            "original_chrom": str,
            "original_start": "int64",
            "original_end": "int64",
            "consensus_chrom": str,
            "consensus_start": "int64",
            "consensus_end": "int64",
        },
    )

    print(f"  {len(mapping_df):,} mapping entries across all samples.")

    # Pre-group by sample to avoid repeatedly filtering the full table.
    mapping_by_sample: dict[str, pd.DataFrame] = {
        sample_id: df.reset_index(drop=True)
        for sample_id, df in mapping_df.groupby("sample_id", sort=False)
    }

    cell_types = sorted(set(CELL_TYPE_STANDARDIZATION.values()))
    summary_records: list[dict[str, Any]] = []

    for cell_type in cell_types:
        print(f"\nAggregating cell type: {cell_type}")

        table, summary = create_pseudobulk_table(
            cell_type=cell_type,
            mapping_by_sample=mapping_by_sample,
            peak_ids=peak_ids,
            consensus_lookup=consensus_lookup,
        )

        output_path = PSEUDOBULK_OUTPUT_DIR / f"{cell_type}_pseudobulk.csv"
        table.to_csv(output_path)

        summary["output_path"] = str(output_path.name)
        summary_records.append(summary)

        print(
            f"  Saved {table.shape[0]:,} peaks × {table.shape[1]} samples "
            f"→ {output_path.name}"
        )

    summary_df = pd.DataFrame.from_records(summary_records)
    summary_df.to_csv(PSEUDOBULK_SUMMARY_PATH, index=False)

    print("\nPseudo-bulk aggregation complete.")
    print(f"Output directory : {PSEUDOBULK_OUTPUT_DIR}")
    print(f"Summary          : {PSEUDOBULK_SUMMARY_PATH}")


if __name__ == "__main__":
    main()
