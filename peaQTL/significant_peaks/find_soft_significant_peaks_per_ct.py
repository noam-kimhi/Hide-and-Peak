#!/usr/bin/env python3
"""Create adjusted-p-value-filtered differential-accessibility BED files.

For each configured cell type, this script reads the saved DESeq2
differential-accessibility results and selects peaks using only:

    padj < DESEQ2_DEFAULT_PADJ_THRESHOLD

No log2-fold-change threshold and no sample-level support filters are
applied.

For every configured cell type, the script writes:

    <cell_type>_significant_peaks.bed.gz

to ``SOFT_UNFILTERED_SIG_PEAKS_DIR``.

Each output file is a gzip-compressed, headerless BED file containing
exactly three columns:

    chrom, start, end

Valid empty BED files are produced for cell types with no qualifying
peaks.
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


REQUIRED_DESEQ2_COLUMNS = {
    "peak_id",
    "padj",
}

PEAK_ID_PATTERN = re.compile(
    r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$"
)


@dataclass(frozen=True)
class BedInterval:
    """A genomic interval parsed from a DESeq2 peak identifier."""

    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class CellTypeSelectionResult:
    """Adjusted-p-value selection result for one cell type."""

    cell_type: str
    total_result_rows: int
    valid_padj_count: int
    missing_padj_count: int
    significant_peak_ids: tuple[str, ...]


def validate_required_constants() -> None:
    """Fail early when constants.py is missing a required value."""
    required_names = (
        "SOFT_UNFILTERED_SIG_PEAKS_DIR",
        "SIGNIFICANT_PEAK_CELL_TYPES",
        "SIGNIFICANT_PEAK_EXPECTED_CELL_TYPE_COUNT",
        "DESEQ2_RESULTS_DIR",
        "DESEQ2_RESULTS_SUFFIX",
        "DESEQ2_DEFAULT_PADJ_THRESHOLD",
    )

    missing = [name for name in required_names if not hasattr(C, name)]
    if missing:
        raise AttributeError(
            "constants.py is missing the following required constants: "
            + ", ".join(missing)
        )

    threshold = float(C.DESEQ2_DEFAULT_PADJ_THRESHOLD)
    if not np.isfinite(threshold) or not 0.0 < threshold <= 1.0:
        raise ValueError(
            "DESEQ2_DEFAULT_PADJ_THRESHOLD must be a finite value in "
            f"(0, 1]; received {threshold!r}."
        )


def validate_cell_types(cell_types: Iterable[str]) -> tuple[str, ...]:
    """Validate and return the configured cell-type names."""
    normalized = tuple(str(cell_type).strip() for cell_type in cell_types)
    expected_count = int(C.SIGNIFICANT_PEAK_EXPECTED_CELL_TYPE_COUNT)

    if len(normalized) != expected_count:
        raise ValueError(
            "SIGNIFICANT_PEAK_CELL_TYPES must contain exactly "
            f"{expected_count} cell types; found {len(normalized)}: "
            f"{normalized}"
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


def load_deseq2_results(path: Path) -> pd.DataFrame:
    """Load and validate one cell type's DESeq2 result table."""
    if not path.exists():
        raise FileNotFoundError(f"Missing DESeq2 result file: {path}")

    dataframe = pd.read_csv(path)

    missing_columns = REQUIRED_DESEQ2_COLUMNS.difference(dataframe.columns)
    if missing_columns:
        raise ValueError(
            f"{path} is missing required columns: "
            f"{', '.join(sorted(missing_columns))}"
        )

    dataframe = dataframe.copy()
    dataframe["peak_id"] = dataframe["peak_id"].astype(str)

    if dataframe["peak_id"].duplicated().any():
        examples = (
            dataframe.loc[
                dataframe["peak_id"].duplicated(keep=False),
                "peak_id",
            ]
            .drop_duplicates()
            .head(5)
            .tolist()
        )
        raise ValueError(
            f"{path} contains duplicate peak IDs, for example: "
            f"{', '.join(examples)}"
        )

    dataframe["padj"] = pd.to_numeric(
        dataframe["padj"],
        errors="coerce",
    )

    valid_padj = (
        np.isfinite(dataframe["padj"])
        & dataframe["padj"].between(0.0, 1.0)
    )
    dataframe.loc[~valid_padj, "padj"] = np.nan

    return dataframe


def select_significant_peaks(results: pd.DataFrame) -> pd.DataFrame:
    """Select peaks using only the configured adjusted-P-value cutoff."""
    threshold = float(C.DESEQ2_DEFAULT_PADJ_THRESHOLD)

    mask = (
        results["padj"].notna()
        & (results["padj"] < threshold)
    )

    return results.loc[mask, ["peak_id", "padj"]].copy()


def process_cell_type(cell_type: str) -> CellTypeSelectionResult:
    """Select adjusted-P-value-significant peaks for one cell type."""
    result_path = (
        Path(C.DESEQ2_RESULTS_DIR)
        / f"{cell_type}{C.DESEQ2_RESULTS_SUFFIX}"
    )

    results = load_deseq2_results(result_path)
    selected = select_significant_peaks(results)

    valid_padj_count = int(results["padj"].notna().sum())
    missing_padj_count = int(results["padj"].isna().sum())

    return CellTypeSelectionResult(
        cell_type=cell_type,
        total_result_rows=len(results),
        valid_padj_count=valid_padj_count,
        missing_padj_count=missing_padj_count,
        significant_peak_ids=tuple(selected["peak_id"].tolist()),
    )


def parse_peak_id(peak_id: str) -> BedInterval:
    """Parse a ``chrom:start-end`` peak identifier."""
    match = PEAK_ID_PATTERN.fullmatch(peak_id)
    if match is None:
        raise ValueError(
            f"Peak ID {peak_id!r} does not match 'chrom:start-end'."
        )

    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))

    if start < 0:
        raise ValueError(
            f"Peak {peak_id!r} has a negative start coordinate."
        )

    if end <= start:
        raise ValueError(
            f"Peak {peak_id!r} has end <= start and is not a valid "
            "BED interval."
        )

    return BedInterval(
        chrom=chrom,
        start=start,
        end=end,
    )


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
    intervals = {
        parse_peak_id(peak_id)
        for peak_id in peak_ids
    }

    return sorted(
        intervals,
        key=lambda interval: (
            chromosome_sort_key(interval.chrom),
            interval.start,
            interval.end,
        ),
    )


def write_bed_gz(
    path: Path,
    intervals: Iterable[BedInterval],
) -> None:
    """Write a three-column, headerless gzip-compressed BED file."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(
        path,
        mode="wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        bed_handle: TextIO = handle

        for interval in intervals:
            bed_handle.write(
                f"{interval.chrom}\t"
                f"{interval.start}\t"
                f"{interval.end}\n"
            )


def print_summary(
    results: Iterable[CellTypeSelectionResult],
) -> None:
    """Print an audit-friendly adjusted-P-value filtering summary."""
    threshold = float(C.DESEQ2_DEFAULT_PADJ_THRESHOLD)

    rows = []
    for result in results:
        significant_count = len(result.significant_peak_ids)

        rows.append(
            {
                "cell_type": result.cell_type,
                "total_DESeq2_rows": result.total_result_rows,
                "valid_padj": result.valid_padj_count,
                "missing_or_invalid_padj": result.missing_padj_count,
                f"padj_below_{threshold:g}": significant_count,
                "selected_percent": (
                    100.0
                    * significant_count
                    / result.valid_padj_count
                    if result.valid_padj_count
                    else 0.0
                ),
            }
        )

    summary = pd.DataFrame(rows)

    print("\nSoft significant-peak selection summary:")
    print(
        summary.to_string(
            index=False,
            formatters={
                "selected_percent": lambda value: f"{value:.2f}%"
            },
        )
    )

    print(
        "\nSelection criterion: "
        f"padj < {threshold:g}. "
        "No log2FoldChange or sample-support filtering was applied."
    )


def main() -> None:
    """Create adjusted-P-value-filtered BED files for all cell types."""
    validate_required_constants()

    cell_types = validate_cell_types(
        C.SIGNIFICANT_PEAK_CELL_TYPES
    )

    # Process every input before writing outputs. This prevents an error in
    # one cell type from leaving a partially updated output directory.
    selections = [
        process_cell_type(cell_type)
        for cell_type in cell_types
    ]

    output_dir = Path(C.SOFT_UNFILTERED_SIG_PEAKS_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    for selection in selections:
        output_path = (
            output_dir
            / f"{selection.cell_type}_significant_peaks.bed.gz"
        )

        intervals = sorted_intervals(
            selection.significant_peak_ids
        )
        write_bed_gz(output_path, intervals)

        print(
            f"Wrote {len(intervals):,} peaks to {output_path}"
        )

    print_summary(selections)
    print(f"\nFinished. Output directory: {output_dir}")


if __name__ == "__main__":
    main()