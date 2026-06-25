#!/usr/bin/env python3
"""Map filtered non-eQTL peaks to candidate target genes using ABC.

Default mode
------------
For every

    <cell_type>_significant_peaks_without_eqtl.bed.gz

inside ``FILTERED_SIG_PEAKS_DIR``, write

    <cell_type>_peak2gene.csv

to ``PEAK2GENE_OUTPUT_DIR``.

Soft mode
---------
When ``--soft`` is supplied, read peaks from
``SOFT_FILTERED_SIG_PEAKS_DIR`` and write outputs to
``SOFT_PEAK2GENE_OUTPUT_DIR``.

Output columns
--------------
    chrom, start, end, TargetGene, ABC.Score

The output coordinates are the input peak coordinates. Intervals use BED
semantics: zero-based, half-open [start, end). When multiple ABC records map
the same peak to the same gene, the largest ABC.Score is retained.

Examples
--------
Default mode:

    python -m peaQTL.peak_to_gene.peak2gene

Soft mode:

    python -m peaQTL.peak_to_gene.peak2gene --soft
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd

from constants import (
    ABC_DICT_PATH,
    FILTERED_SIG_PEAKS_DIR,
    FILTERED_SIGNIFICANT_PEAKS_SUFFIX,
    PEAK2GENE_OUTPUT_DIR,
    SIGNIFICANT_PEAK_EXPECTED_CELL_TYPE_COUNT,
    SOFT_FILTERED_SIG_PEAKS_DIR,
    SOFT_PEAK2GENE_OUTPUT_DIR,
)


LOGGER: Final[logging.Logger] = logging.getLogger(__name__)

BED_COLUMNS: Final[tuple[str, str, str]] = (
    "chrom",
    "start",
    "end",
)

ABC_COLUMNS: Final[tuple[str, str, str, str, str]] = (
    "chrom",
    "start",
    "end",
    "TargetGene",
    "ABC.Score",
)

OUTPUT_COLUMNS: Final[list[str]] = list(ABC_COLUMNS)

SUMMARY_FILENAME: Final[str] = "peak2gene_mapping_summary.csv"

SUMMARY_COLUMNS: Final[list[str]] = [
    "cell_type",
    "total_peaks",
    "mapped_peaks",
    "mapped_peak_percent",
    "peak_gene_pairs",
]


def configure_logging() -> None:
    """Configure command-line logging."""

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Map filtered non-eQTL significant peaks to candidate target "
            "genes using the ABC dictionary."
        )
    )

    parser.add_argument(
        "--soft",
        action="store_true",
        help=(
            "Read peaks from SOFT_FILTERED_SIG_PEAKS_DIR and write outputs "
            "to SOFT_PEAK2GENE_OUTPUT_DIR."
        ),
    )

    return parser.parse_args()


def resolve_mode_paths(
    soft: bool,
) -> tuple[Path, Path]:
    """Return the input and output directories for the selected mode."""

    if soft:
        return (
            SOFT_FILTERED_SIG_PEAKS_DIR,
            SOFT_PEAK2GENE_OUTPUT_DIR,
        )

    return (
        FILTERED_SIG_PEAKS_DIR,
        PEAK2GENE_OUTPUT_DIR,
    )


def chromosome_key(value: object) -> str:
    """Normalize chromosome labels only for matching."""

    chromosome = str(value).strip()

    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]

    chromosome = chromosome.upper()

    return "M" if chromosome == "MT" else chromosome


def validate_intervals(
    dataframe: pd.DataFrame,
    source: Path,
) -> pd.DataFrame:
    """Validate BED-like coordinates and add an internal chromosome key."""

    result = dataframe.copy()

    if result.empty:
        result["_chrom_key"] = pd.Series(dtype="object")
        return result

    if result[list(BED_COLUMNS)].isna().any().any():
        raise ValueError(
            f"{source} contains missing interval values."
        )

    result["chrom"] = (
        result["chrom"]
        .astype(str)
        .str.strip()
    )

    if result["chrom"].eq("").any():
        raise ValueError(
            f"{source} contains empty chromosome names."
        )

    for column in ("start", "end"):
        values = pd.to_numeric(
            result[column],
            errors="raise",
        )

        if not np.isfinite(
            values.to_numpy(dtype=float)
        ).all():
            raise ValueError(
                f"{source} contains non-finite {column} values."
            )

        if not np.equal(
            values,
            np.floor(values),
        ).all():
            raise ValueError(
                f"{source} contains non-integer {column} values."
            )

        result[column] = values.astype("int64")

    invalid = (
        (result["start"] < 0)
        | (result["end"] <= result["start"])
    )

    if invalid.any():
        example = (
            result.loc[
                invalid,
                list(BED_COLUMNS),
            ]
            .iloc[0]
            .to_dict()
        )

        raise ValueError(
            f"{source} contains an invalid interval: {example}"
        )

    result["_chrom_key"] = result["chrom"].map(
        chromosome_key
    )

    return result


def read_abc_dictionary(path: Path) -> pd.DataFrame:
    """Read and validate the ABC dictionary once for all cell types."""

    if not path.is_file():
        raise FileNotFoundError(
            f"ABC dictionary not found: {path}"
        )

    available_columns = set(
        pd.read_csv(path, nrows=0).columns
    )

    missing = sorted(
        set(ABC_COLUMNS) - available_columns
    )

    if missing:
        raise ValueError(
            f"{path} is missing required columns: {missing}"
        )

    abc = pd.read_csv(
        path,
        usecols=OUTPUT_COLUMNS,
    )

    abc = validate_intervals(
        abc,
        path,
    )

    if abc.empty:
        raise ValueError(
            f"ABC dictionary is empty: {path}"
        )

    if abc[
        ["TargetGene", "ABC.Score"]
    ].isna().any().any():
        raise ValueError(
            f"{path} contains missing ABC annotations."
        )

    abc["TargetGene"] = (
        abc["TargetGene"]
        .astype(str)
        .str.strip()
    )

    abc["ABC.Score"] = pd.to_numeric(
        abc["ABC.Score"],
        errors="raise",
    )

    if abc["TargetGene"].eq("").any():
        raise ValueError(
            f"{path} contains empty TargetGene values."
        )

    if not np.isfinite(
        abc["ABC.Score"].to_numpy(dtype=float)
    ).all():
        raise ValueError(
            f"{path} contains non-finite ABC.Score values."
        )

    abc = abc.reset_index(drop=True)
    abc["_abc_id"] = range(len(abc))

    return abc


def read_bed3(path: Path) -> pd.DataFrame:
    """Read the first three columns of a plain or gzipped BED file."""

    try:
        peaks = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=list(BED_COLUMNS),
            usecols=[0, 1, 2],
            comment="#",
            compression="infer",
        )

    except pd.errors.EmptyDataError:
        peaks = pd.DataFrame(
            columns=BED_COLUMNS
        )

    peaks = validate_intervals(
        peaks,
        path,
    )

    duplicate_count = int(
        peaks.duplicated(
            list(BED_COLUMNS)
        ).sum()
    )

    if duplicate_count:
        LOGGER.warning(
            "%s contains %d duplicate peaks; keeping one copy",
            path.name,
            duplicate_count,
        )

        peaks = peaks.drop_duplicates(
            list(BED_COLUMNS)
        )

    peaks = peaks.reset_index(drop=True)
    peaks["_peak_id"] = range(len(peaks))

    return peaks


def find_overlap_pairs(
    peaks: pd.DataFrame,
    abc: pd.DataFrame,
) -> list[tuple[int, int]]:
    """Return all overlapping ``(peak_id, abc_id)`` pairs.

    A chromosome-wise sweep line is used. End events are handled before
    start events at the same coordinate, so book-ended intervals do not
    overlap.
    """

    if peaks.empty:
        return []

    abc_by_chromosome = {
        chromosome: group
        for chromosome, group in abc.groupby(
            "_chrom_key",
            sort=False,
        )
    }

    pairs: list[tuple[int, int]] = []

    # Event order:
    # 0 = peak end
    # 1 = ABC end
    # 2 = peak start
    # 3 = ABC start
    for chromosome, chromosome_peaks in peaks.groupby(
        "_chrom_key",
        sort=False,
    ):
        chromosome_abc = abc_by_chromosome.get(
            chromosome
        )

        if chromosome_abc is None:
            continue

        events: list[tuple[int, int, int]] = []

        for peak_id, start, end in chromosome_peaks[
            ["_peak_id", "start", "end"]
        ].itertuples(index=False, name=None):
            events.extend(
                (
                    (
                        int(end),
                        0,
                        int(peak_id),
                    ),
                    (
                        int(start),
                        2,
                        int(peak_id),
                    ),
                )
            )

        for abc_id, start, end in chromosome_abc[
            ["_abc_id", "start", "end"]
        ].itertuples(index=False, name=None):
            events.extend(
                (
                    (
                        int(end),
                        1,
                        int(abc_id),
                    ),
                    (
                        int(start),
                        3,
                        int(abc_id),
                    ),
                )
            )

        events.sort()

        active_peaks: set[int] = set()
        active_abc: set[int] = set()

        for _position, event_type, interval_id in events:
            if event_type == 0:
                active_peaks.discard(interval_id)

            elif event_type == 1:
                active_abc.discard(interval_id)

            elif event_type == 2:
                pairs.extend(
                    (interval_id, abc_id)
                    for abc_id in active_abc
                )

                active_peaks.add(interval_id)

            else:
                pairs.extend(
                    (peak_id, interval_id)
                    for peak_id in active_peaks
                )

                active_abc.add(interval_id)

    return pairs


def build_output(
    peaks: pd.DataFrame,
    abc: pd.DataFrame,
    pairs: list[tuple[int, int]],
) -> pd.DataFrame:
    """Build one row per unique peak/target-gene mapping."""

    if not pairs:
        return pd.DataFrame(
            columns=OUTPUT_COLUMNS
        )

    peak_ids = [
        peak_id
        for peak_id, _ in pairs
    ]

    abc_ids = [
        abc_id
        for _, abc_id in pairs
    ]

    matched_peaks = (
        peaks
        .set_index("_peak_id")
        .loc[peak_ids]
    )

    matched_abc = (
        abc
        .set_index("_abc_id")
        .loc[abc_ids]
    )

    output = pd.DataFrame(
        {
            "_peak_id": peak_ids,
            "chrom": matched_peaks[
                "chrom"
            ].to_numpy(),
            "start": matched_peaks[
                "start"
            ].to_numpy(),
            "end": matched_peaks[
                "end"
            ].to_numpy(),
            "TargetGene": matched_abc[
                "TargetGene"
            ].to_numpy(),
            "ABC.Score": matched_abc[
                "ABC.Score"
            ].to_numpy(),
        }
    )

    # Multiple ABC records can map the same peak to the same gene.
    # Retain the strongest ABC prediction.
    output = output.sort_values(
        [
            "_peak_id",
            "TargetGene",
            "ABC.Score",
        ],
        ascending=[
            True,
            True,
            False,
        ],
        kind="mergesort",
    )

    output = output.drop_duplicates(
        [
            "_peak_id",
            "TargetGene",
        ],
        keep="first",
    )

    output = output.sort_values(
        [
            "_peak_id",
            "TargetGene",
        ],
        kind="mergesort",
    )

    return (
        output[OUTPUT_COLUMNS]
        .reset_index(drop=True)
    )


def cell_type_from_filename(path: Path) -> str:
    """Extract ``<cell_type>`` from the configured input filename."""

    if not path.name.endswith(
        FILTERED_SIGNIFICANT_PEAKS_SUFFIX
    ):
        raise ValueError(
            f"Unexpected peak filename: {path.name}"
        )

    cell_type = path.name[
        : -len(
            FILTERED_SIGNIFICANT_PEAKS_SUFFIX
        )
    ]

    if not cell_type:
        raise ValueError(
            f"Could not extract cell type from {path.name}"
        )

    return cell_type


def process_peak_file(
    path: Path,
    abc: pd.DataFrame,
    output_dir: Path,
) -> dict[str, object]:
    """Create one cell type's output and return mapping statistics."""

    cell_type = cell_type_from_filename(
        path
    )

    peaks = read_bed3(
        path
    )

    pairs = find_overlap_pairs(
        peaks,
        abc,
    )

    output = build_output(
        peaks,
        abc,
        pairs,
    )

    output_path = (
        output_dir
        / f"{cell_type}_peak2gene.csv"
    )

    output.to_csv(
        output_path,
        index=False,
    )

    total_peaks = len(peaks)

    mapped_peaks = (
        output[
            list(BED_COLUMNS)
        ]
        .drop_duplicates()
        .shape[0]
        if not output.empty
        else 0
    )

    mapped_peak_percent = (
        100.0 * mapped_peaks / total_peaks
        if total_peaks > 0
        else 0.0
    )

    LOGGER.info(
        (
            "%s: %d input peaks, %d mapped peaks (%.2f%%), "
            "%d peak-gene pairs -> %s"
        ),
        cell_type,
        total_peaks,
        mapped_peaks,
        mapped_peak_percent,
        len(output),
        output_path,
    )

    return {
        "cell_type": cell_type,
        "total_peaks": total_peaks,
        "mapped_peaks": mapped_peaks,
        "mapped_peak_percent": mapped_peak_percent,
        "peak_gene_pairs": len(output),
    }


def run(soft: bool = False) -> pd.DataFrame:
    """Run peak-to-gene mapping in default or soft mode."""

    input_dir, output_dir = resolve_mode_paths(
        soft
    )

    mode_name = (
        "soft"
        if soft
        else "default"
    )

    LOGGER.info(
        "Running peak2gene in %s mode",
        mode_name,
    )

    LOGGER.info(
        "Input directory: %s",
        input_dir,
    )

    LOGGER.info(
        "Output directory: %s",
        output_dir,
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    abc = read_abc_dictionary(
        ABC_DICT_PATH
    )

    peak_paths = sorted(
        path
        for path in input_dir.glob(
            f"*{FILTERED_SIGNIFICANT_PEAKS_SUFFIX}"
        )
        if path.is_file()
    )

    if not peak_paths:
        raise FileNotFoundError(
            f"No *{FILTERED_SIGNIFICANT_PEAKS_SUFFIX} files found in "
            f"{input_dir}"
        )

    if (
        len(peak_paths)
        != SIGNIFICANT_PEAK_EXPECTED_CELL_TYPE_COUNT
    ):
        LOGGER.warning(
            "Expected %d cell-type files but found %d",
            SIGNIFICANT_PEAK_EXPECTED_CELL_TYPE_COUNT,
            len(peak_paths),
        )

    LOGGER.info(
        "Loaded %d ABC predictions",
        len(abc),
    )

    summary_rows = [
        process_peak_file(
            path=peak_path,
            abc=abc,
            output_dir=output_dir,
        )
        for peak_path in peak_paths
    ]

    summary = pd.DataFrame(
        summary_rows,
        columns=SUMMARY_COLUMNS,
    )

    summary = summary.sort_values(
        "cell_type",
        kind="mergesort",
    )

    summary_path = (
        output_dir
        / SUMMARY_FILENAME
    )

    summary.to_csv(
        summary_path,
        index=False,
        float_format="%.4f",
    )

    LOGGER.info(
        "Wrote %s-mode mapping summary to %s",
        mode_name,
        summary_path,
    )

    return summary


def main() -> None:
    """Command-line entry point."""

    configure_logging()

    args = parse_arguments()

    try:
        run(
            soft=args.soft
        )

    except Exception:
        LOGGER.exception(
            "Peak-to-gene mapping failed"
        )
        raise


if __name__ == "__main__":
    main()