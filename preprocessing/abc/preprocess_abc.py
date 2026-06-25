#!/usr/bin/env python3
"""Create an hg38 liver ABC dictionary used by peak-to-gene analyses.

The source ABC prediction table is distributed in hg19 coordinates. This
script loads the required prediction columns, converts each unique enhancer
interval to hg38 with the UCSC ``liftOver`` executable, and saves the complete
lifted prediction table to ``ABC_38_PRED_FILE``. The original hg19 coordinates
are retained in explicit ``*_hg19`` columns.

Only after coordinate conversion does the script retain liver- or
hepatocyte-related predictions, apply the configured ABC-score threshold,
rename ``chr`` to ``chrom`` for bioframe compatibility, and write the reduced
liver dictionary to ``ABC_DICT_PATH``.
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Final

import pandas as pd

from constants import (
    ABC_38_PRED_FILE,
    ABC_DICT_PATH,
    ABC_HG19_TO_HG38_CHAIN_PATH,
    ABC_LIFTOVER_EXECUTABLE_PATH,
    ABC_PRED_FILE,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
LOGGER: Final[logging.Logger] = logging.getLogger(__name__)

REQUIRED_COLUMNS: Final[list[str]] = [
    "chr",
    "start",
    "end",
    "TargetGene",
    "CellType",
    "ABC.Score",
]
LIFTOVER_ID_COLUMN: Final[str] = "liftover_id"
ROW_ORDER_COLUMN: Final[str] = "abc_row_order"
ABC_DICTIONARY_COLUMNS: Final[list[str]] = [
    "chrom",
    "start",
    "end",
    "TargetGene",
    "CellType",
    "ABC.Score",
]
ABC_SCORE_THRESHOLD: Final[float] = 0.02
LIVER_KEYWORDS: Final[tuple[str, ...]] = ("liver", "hepatocyte", "HepG2")


def validate_input_file(path: Path, description: str) -> None:
    """Validate that a required input file exists.

    :param path: Path expected to reference an existing file.
    :param description: Human-readable description used in error messages.
    """
    if not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")


def validate_liftover_executable(path: Path) -> None:
    """Validate that the configured UCSC liftOver executable can be run.

    :param path: Path to the UCSC ``liftOver`` executable.
    """
    validate_input_file(path, "UCSC liftOver executable")
    if not os.access(path, os.X_OK):
        raise PermissionError(
            f"UCSC liftOver is not executable: {path}. "
            f"Run: chmod +x {path}"
        )


def load_abc_predictions(path: Path) -> pd.DataFrame:
    """Load the ABC columns required by the downstream pipeline.

    :param path: Path to the gzip-compressed, tab-separated ABC prediction file.
    :return: DataFrame containing the required ABC prediction columns.
    """
    validate_input_file(path, "ABC prediction file")
    LOGGER.info("Loading ABC predictions from: %s", path)
    predictions = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        usecols=REQUIRED_COLUMNS,
        low_memory=False,
    )
    LOGGER.info("Loaded %d ABC prediction rows.", len(predictions))
    return predictions


def validate_abc_coordinates(predictions: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize the hg19 enhancer coordinates.

    :param predictions: ABC prediction rows containing ``chr``, ``start``, and
        ``end`` columns.
    :return: Copy of the input with integer ``start`` and ``end`` coordinates.
    """
    validated = predictions.copy()

    if validated["chr"].isna().any():
        raise ValueError("The ABC prediction table contains missing chromosomes.")

    validated["chr"] = validated["chr"].astype(str)
    validated["start"] = pd.to_numeric(validated["start"], errors="raise")
    validated["end"] = pd.to_numeric(validated["end"], errors="raise")

    non_integer_start = validated["start"] % 1 != 0
    non_integer_end = validated["end"] % 1 != 0
    if non_integer_start.any() or non_integer_end.any():
        raise ValueError("ABC start and end coordinates must be integers.")

    validated["start"] = validated["start"].astype("int64")
    validated["end"] = validated["end"].astype("int64")

    invalid_coordinates = (validated["start"] < 0) | (
        validated["end"] <= validated["start"]
    )
    if invalid_coordinates.any():
        invalid_count = int(invalid_coordinates.sum())
        raise ValueError(
            f"The ABC prediction table contains {invalid_count} invalid intervals."
        )

    return validated


def create_unique_liftover_intervals(predictions: pd.DataFrame) -> pd.DataFrame:
    """Create one BED record for each unique hg19 enhancer interval.

    :param predictions: Validated ABC prediction rows in hg19 coordinates.
    :return: Unique hg19 intervals with a stable identifier for joining the
        lifted coordinates back to all prediction rows.
    """
    intervals = (
        predictions[["chr", "start", "end"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    intervals[LIFTOVER_ID_COLUMN] = [
        f"abc_interval_{index}" for index in range(len(intervals))
    ]
    LOGGER.info(
        "Prepared %d unique hg19 intervals from %d ABC rows.",
        len(intervals),
        len(predictions),
    )
    return intervals


def write_liftover_input(intervals: pd.DataFrame, path: Path) -> None:
    """Write unique hg19 intervals as a BED4 file for UCSC liftOver.

    :param intervals: Unique hg19 intervals with stable liftOver identifiers.
    :param path: Destination temporary BED4 path.
    """
    intervals[["chr", "start", "end", LIFTOVER_ID_COLUMN]].to_csv(
        path,
        sep="\t",
        header=False,
        index=False,
    )


def run_liftover(
    input_bed_path: Path,
    output_bed_path: Path,
    unmapped_bed_path: Path,
) -> None:
    """Run UCSC liftOver from hg19 to hg38.

    :param input_bed_path: Temporary BED4 file containing hg19 intervals.
    :param output_bed_path: Destination BED4 file for successfully lifted hg38
        intervals.
    :param unmapped_bed_path: Destination file for intervals that could not be
        lifted.
    """
    validate_liftover_executable(ABC_LIFTOVER_EXECUTABLE_PATH)
    validate_input_file(
        ABC_HG19_TO_HG38_CHAIN_PATH,
        "hg19-to-hg38 liftOver chain file",
    )

    command = [
        str(ABC_LIFTOVER_EXECUTABLE_PATH),
        str(input_bed_path),
        str(ABC_HG19_TO_HG38_CHAIN_PATH),
        str(output_bed_path),
        str(unmapped_bed_path),
    ]
    LOGGER.info("Lifting ABC enhancer intervals from hg19 to hg38...")
    subprocess.run(command, check=True)


def load_lifted_coordinates(path: Path) -> pd.DataFrame:
    """Load successfully converted hg38 BED4 intervals.

    :param path: BED4 file written by UCSC liftOver.
    :return: DataFrame containing hg38 coordinates and liftOver identifiers.
    """
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError("UCSC liftOver produced no mapped hg38 intervals.")

    lifted = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chr_hg38", "start_hg38", "end_hg38", LIFTOVER_ID_COLUMN],
        dtype={
            "chr_hg38": str,
            "start_hg38": "int64",
            "end_hg38": "int64",
            LIFTOVER_ID_COLUMN: str,
        },
    )

    duplicated_ids = lifted[LIFTOVER_ID_COLUMN].duplicated(keep=False)
    if duplicated_ids.any():
        duplicated_count = int(
            lifted.loc[duplicated_ids, LIFTOVER_ID_COLUMN].nunique()
        )
        raise RuntimeError(
            "UCSC liftOver returned multiple mappings for "
            f"{duplicated_count} input intervals; ambiguous mappings are not "
            "accepted by this preprocessing step."
        )

    return lifted


def apply_lifted_coordinates(
    predictions: pd.DataFrame,
    intervals: pd.DataFrame,
    lifted_coordinates: pd.DataFrame,
) -> pd.DataFrame:
    """Replace hg19 coordinates with their unambiguous hg38 mappings.

    Rows whose enhancer interval could not be lifted are excluded. Their
    original hg19 coordinates remain available in ``chr_hg19``, ``start_hg19``,
    and ``end_hg19``.

    :param predictions: Validated ABC prediction rows in hg19 coordinates.
    :param intervals: Unique hg19 intervals and their stable identifiers.
    :param lifted_coordinates: Successful hg38 mappings returned by liftOver.
    :return: ABC prediction rows using hg38 coordinates, with original hg19
        coordinates retained in separate columns.
    """
    coordinate_map = intervals.merge(
        lifted_coordinates,
        on=LIFTOVER_ID_COLUMN,
        how="inner",
        validate="one_to_one",
        sort=False,
    )

    predictions_with_order = predictions.copy()
    predictions_with_order[ROW_ORDER_COLUMN] = range(len(predictions_with_order))

    lifted_predictions = predictions_with_order.merge(
        coordinate_map[
            [
                "chr",
                "start",
                "end",
                "chr_hg38",
                "start_hg38",
                "end_hg38",
            ]
        ],
        on=["chr", "start", "end"],
        how="inner",
        validate="many_to_one",
        sort=False,
    )
    lifted_predictions = lifted_predictions.sort_values(
        ROW_ORDER_COLUMN,
        kind="stable",
    ).drop(columns=[ROW_ORDER_COLUMN])

    lifted_predictions = lifted_predictions.rename(
        columns={
            "chr": "chr_hg19",
            "start": "start_hg19",
            "end": "end_hg19",
            "chr_hg38": "chr",
            "start_hg38": "start",
            "end_hg38": "end",
        }
    )

    output_columns = [
        "chr",
        "start",
        "end",
        "chr_hg19",
        "start_hg19",
        "end_hg19",
        "TargetGene",
        "CellType",
        "ABC.Score",
    ]
    lifted_predictions = lifted_predictions[output_columns]

    mapped_unique_intervals = len(coordinate_map)
    total_unique_intervals = len(intervals)
    mapped_rows = len(lifted_predictions)
    total_rows = len(predictions)
    LOGGER.info(
        "liftOver mapped %d/%d unique intervals (%.2f%%).",
        mapped_unique_intervals,
        total_unique_intervals,
        100.0 * mapped_unique_intervals / total_unique_intervals,
    )
    LOGGER.info(
        "Retained %d/%d ABC rows after liftOver (%.2f%%).",
        mapped_rows,
        total_rows,
        100.0 * mapped_rows / total_rows,
    )
    return lifted_predictions


def lift_abc_predictions_to_hg38(predictions: pd.DataFrame) -> pd.DataFrame:
    """Lift all unique ABC enhancer intervals from hg19 to hg38.

    :param predictions: Validated ABC prediction rows in hg19 coordinates.
    :return: ABC prediction rows in hg38 coordinates.
    """
    intervals = create_unique_liftover_intervals(predictions)

    with TemporaryDirectory(prefix="abc_liftover_") as temporary_directory:
        temporary_path = Path(temporary_directory)
        input_bed_path = temporary_path / "abc_hg19.bed"
        output_bed_path = temporary_path / "abc_hg38.bed"
        unmapped_bed_path = temporary_path / "abc_hg19_unmapped.bed"

        write_liftover_input(intervals, input_bed_path)
        run_liftover(
            input_bed_path,
            output_bed_path,
            unmapped_bed_path,
        )
        lifted_coordinates = load_lifted_coordinates(output_bed_path)

    return apply_lifted_coordinates(
        predictions,
        intervals,
        lifted_coordinates,
    )


def save_lifted_predictions(predictions: pd.DataFrame, path: Path) -> None:
    """Save the complete lifted ABC prediction table as a compressed TSV.

    :param predictions: ABC prediction rows in hg38 coordinates.
    :param path: Destination ``.txt.gz`` path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Saving %d lifted ABC rows to: %s", len(predictions), path)
    predictions.to_csv(
        path,
        sep="\t",
        compression="gzip",
        index=False,
    )


def filter_liver_cell_types(predictions: pd.DataFrame) -> pd.DataFrame:
    """Retain predictions whose cell-type label is liver related.

    :param predictions: Lifted ABC prediction rows.
    :return: Rows whose ``CellType`` contains a configured liver keyword.
    """
    LOGGER.info("Filtering for liver/hepatocyte cell types...")
    pattern = "|".join(LIVER_KEYWORDS)
    mask = predictions["CellType"].str.contains(
        pattern,
        case=False,
        na=False,
        regex=True,
    )
    filtered = predictions.loc[mask].copy()
    LOGGER.info(
        "Rows after tissue filter: %d (removed %d).",
        len(filtered),
        len(predictions) - len(filtered),
    )
    return filtered


def filter_by_abc_score(predictions: pd.DataFrame) -> pd.DataFrame:
    """Retain predictions meeting the configured ABC-score threshold.

    :param predictions: Liver-related ABC prediction rows.
    :return: Rows whose ``ABC.Score`` is at least ``ABC_SCORE_THRESHOLD``.
    """
    LOGGER.info(
        "Filtering rows with ABC.Score >= %.4f...",
        ABC_SCORE_THRESHOLD,
    )
    filtered = predictions.loc[
        predictions["ABC.Score"] >= ABC_SCORE_THRESHOLD
    ].copy()
    LOGGER.info(
        "Rows after score filter: %d (removed %d).",
        len(filtered),
        len(predictions) - len(filtered),
    )
    return filtered


def rename_columns_for_bioframe(predictions: pd.DataFrame) -> pd.DataFrame:
    """Rename the hg38 chromosome column for bioframe compatibility.

    :param predictions: Filtered ABC prediction rows in hg38 coordinates.
    :return: Copy with ``chr`` renamed to ``chrom``.
    """
    LOGGER.info("Renaming 'chr' to 'chrom' for bioframe compatibility...")
    renamed = predictions.rename(columns={"chr": "chrom"})
    return renamed[ABC_DICTIONARY_COLUMNS].copy()


def save_abc_dictionary(predictions: pd.DataFrame, path: Path) -> None:
    """Write the filtered hg38 liver ABC dictionary as CSV.

    :param predictions: Filtered liver ABC prediction rows in hg38 coordinates.
    :param path: Destination CSV path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Saving %d liver ABC rows to: %s", len(predictions), path)
    predictions.to_csv(path, index=False)
    LOGGER.info("ABC preprocessing completed successfully.")


def main() -> None:
    """Lift ABC predictions to hg38 and create the liver dictionary."""
    predictions = load_abc_predictions(ABC_PRED_FILE)
    predictions = validate_abc_coordinates(predictions)
    lifted_predictions = lift_abc_predictions_to_hg38(predictions)
    save_lifted_predictions(lifted_predictions, ABC_38_PRED_FILE)

    liver_predictions = filter_liver_cell_types(lifted_predictions)
    liver_predictions = filter_by_abc_score(liver_predictions)
    liver_predictions = rename_columns_for_bioframe(liver_predictions)
    save_abc_dictionary(liver_predictions, ABC_DICT_PATH)


if __name__ == "__main__":
    main()
