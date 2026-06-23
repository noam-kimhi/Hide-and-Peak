"""Remove significant peaks that overlap exact GTEx liver eQTL positions.

For every file in ``UNFILTERED_SIG_PEAKS_DIR`` named

    <cell_type>_significant_peaks.bed.gz

the script writes

    <cell_type>_significant_peaks_without_eqtl.bed.gz

to ``FILTERED_SIG_PEAKS_DIR``.

A complete input peak is removed when it overlaps at least one significant
GTEx liver eQTL position. GTEx b38 variant positions are converted from
one-based coordinates to one-base, zero-based BED intervals [position-1,
position). Peaks are never shortened or split.

Run from the project root with:

    python -m preprocessing.drop_eqtl_from_significant_peaks

Use ``--overwrite`` to replace existing outputs.
"""

from __future__ import annotations

import argparse
import gzip
import logging
from pathlib import Path
from typing import Final, Sequence

import bioframe
import pandas as pd

from constants import (
    FILTERED_SIG_PEAKS_DIR,
    FILTERED_SIGNIFICANT_PEAKS_SUFFIX,
    LIVER_EQTL_SIGNIFICANT_PAIRS_PATH,
    SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH,
    UNFILTERED_SIG_PEAKS_DIR,
    UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX,
)


LOGGER: Final[logging.Logger] = logging.getLogger(__name__)

BED_MIN_COLUMNS: Final[int] = 3
CHROM_COLUMN_INDEX: Final[int] = 0
START_COLUMN_INDEX: Final[int] = 1
END_COLUMN_INDEX: Final[int] = 2

VARIANT_ID_COLUMN: Final[str] = "variant_id"
EXPECTED_GTEX_BUILD: Final[str] = "b38"

BIOFRAME_CHROM_COLUMN: Final[str] = "chrom"
BIOFRAME_START_COLUMN: Final[str] = "start"
BIOFRAME_END_COLUMN: Final[str] = "end"
BIOFRAME_ROW_ID_COLUMN: Final[str] = "_original_row_id"

SUMMARY_COLUMNS: Final[Sequence[str]] = (
    "cell_type",
    "input_file",
    "output_file",
    "n_input_significant_peaks",
    "n_removed_eqtl_overlapping_peaks",
    "n_significant_peaks_without_eqtl",
    "fraction_removed",
)


def configure_logging(log_level: str) -> None:
    """Configure command-line logging."""

    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def normalize_chromosome(chromosome: object) -> str:
    """Convert common chromosome spellings to a canonical ``chr*`` form.

    The normalization is used only for interval matching. Original BED
    chromosome values are preserved in the written output.
    """

    value = str(chromosome).strip()
    if not value:
        raise ValueError("Encountered an empty chromosome value.")

    suffix = value[3:] if value.lower().startswith("chr") else value
    if not suffix:
        raise ValueError(f"Invalid chromosome value: {value!r}")

    upper_suffix = suffix.upper()
    if upper_suffix in {"M", "MT"}:
        return "chrM"
    if upper_suffix in {"X", "Y"}:
        return f"chr{upper_suffix}"
    if suffix.isdigit():
        return f"chr{int(suffix)}"

    # Preserve non-standard contig text after ensuring one ``chr`` prefix.
    return f"chr{suffix}"


def parse_gtex_variant_ids(variant_ids: pd.Series) -> pd.DataFrame:
    """Parse GTEx variant IDs and return unique one-base GRCh38 intervals.

    Expected variant format:

        chr1_14677_G_A_b38

    The split is performed from the right so that chromosome/contig names
    containing underscores remain parseable.
    """

    ids = variant_ids.astype("string").str.strip()
    if ids.isna().any() or (ids == "").any():
        raise ValueError("The GTEx eQTL table contains empty variant IDs.")

    parsed = ids.str.rsplit("_", n=4, expand=True)
    if parsed.shape[1] != 5:
        examples = ids.head(5).tolist()
        raise ValueError(
            "Could not parse GTEx variant IDs into "
            "chromosome, position, reference, alternate, and build fields. "
            f"Examples: {examples}"
        )

    parsed.columns = ("chrom_raw", "position_raw", "ref", "alt", "build")
    positions = pd.to_numeric(parsed["position_raw"], errors="coerce")

    invalid_mask = (
        positions.isna()
        | (positions < 1)
        | parsed["chrom_raw"].isna()
        | (parsed["chrom_raw"].astype("string").str.strip() == "")
    )
    if invalid_mask.any():
        invalid_examples = ids.loc[invalid_mask].head(10).tolist()
        raise ValueError(
            "Found malformed GTEx variant IDs. "
            f"Examples: {invalid_examples}"
        )

    unexpected_builds = sorted(
        set(parsed.loc[parsed["build"] != EXPECTED_GTEX_BUILD, "build"].dropna())
    )
    if unexpected_builds:
        raise ValueError(
            "The GTEx eQTL file contains variants outside GRCh38. "
            f"Expected only {EXPECTED_GTEX_BUILD!r}; found {unexpected_builds}."
        )

    positions = positions.astype("int64")
    intervals = pd.DataFrame(
        {
            BIOFRAME_CHROM_COLUMN: parsed["chrom_raw"].map(
                normalize_chromosome
            ),
            BIOFRAME_START_COLUMN: positions - 1,
            BIOFRAME_END_COLUMN: positions,
        }
    )

    # Several significant variant-gene pairs can refer to the same variant,
    # and multiple alleles can share a genomic position. For peak removal,
    # one interval per genomic position is sufficient.
    intervals = intervals.drop_duplicates(
        subset=(
            BIOFRAME_CHROM_COLUMN,
            BIOFRAME_START_COLUMN,
            BIOFRAME_END_COLUMN,
        ),
        keep="first",
    ).reset_index(drop=True)

    return intervals


def load_eqtl_intervals(eqtl_path: Path) -> tuple[pd.DataFrame, int]:
    """Load significant liver eQTLs and return unique genomic positions."""

    if not eqtl_path.is_file():
        raise FileNotFoundError(f"GTEx liver eQTL file not found: {eqtl_path}")

    LOGGER.info("Reading GTEx liver eQTLs from %s", eqtl_path)
    eqtl_table = pd.read_csv(
        eqtl_path,
        sep="\t",
        compression="infer",
        usecols=[VARIANT_ID_COLUMN],
        dtype={VARIANT_ID_COLUMN: "string"},
    )

    n_unique_variant_ids = int(eqtl_table[VARIANT_ID_COLUMN].nunique())
    unique_variant_ids = eqtl_table[VARIANT_ID_COLUMN].drop_duplicates()
    eqtl_intervals = parse_gtex_variant_ids(unique_variant_ids)

    LOGGER.info(
        "Loaded %d unique variant IDs at %d unique genomic positions",
        n_unique_variant_ids,
        len(eqtl_intervals),
    )
    return eqtl_intervals, n_unique_variant_ids


def read_bed_preserving_values(bed_path: Path) -> pd.DataFrame:
    """Read a headerless BED-like file while preserving original values."""

    try:
        bed = pd.read_csv(
            bed_path,
            sep="\t",
            header=None,
            comment="#",
            dtype="string",
            keep_default_na=False,
            compression="infer",
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    if bed.shape[1] < BED_MIN_COLUMNS:
        raise ValueError(
            f"{bed_path} has {bed.shape[1]} columns; expected at least "
            f"{BED_MIN_COLUMNS} BED columns."
        )

    starts = pd.to_numeric(
        bed.iloc[:, START_COLUMN_INDEX],
        errors="coerce",
    )
    ends = pd.to_numeric(
        bed.iloc[:, END_COLUMN_INDEX],
        errors="coerce",
    )

    invalid_coordinate_mask = (
        starts.isna()
        | ends.isna()
        | (starts < 0)
        | (ends <= starts)
    )
    if invalid_coordinate_mask.any():
        invalid_rows = (
            invalid_coordinate_mask[invalid_coordinate_mask]
            .index[:10]
            .tolist()
        )
        raise ValueError(
            f"{bed_path} contains invalid BED coordinates at zero-based "
            f"row indices {invalid_rows}."
        )

    return bed


def build_bioframe_peak_table(peaks: pd.DataFrame) -> pd.DataFrame:
    """Create a normalized interval table linked to original BED rows."""

    return pd.DataFrame(
        {
            BIOFRAME_CHROM_COLUMN: peaks.iloc[:, CHROM_COLUMN_INDEX].map(
                normalize_chromosome
            ),
            BIOFRAME_START_COLUMN: pd.to_numeric(
                peaks.iloc[:, START_COLUMN_INDEX],
                errors="raise",
            ).astype("int64"),
            BIOFRAME_END_COLUMN: pd.to_numeric(
                peaks.iloc[:, END_COLUMN_INDEX],
                errors="raise",
            ).astype("int64"),
            BIOFRAME_ROW_ID_COLUMN: range(len(peaks)),
        }
    )


def write_bed_atomic(bed: pd.DataFrame, output_path: Path) -> None:
    """Write a gzipped, headerless BED file atomically."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f"{output_path.name}.tmp")

    try:
        with gzip.open(temporary_path, mode="wt", newline="") as handle:
            if not bed.empty:
                bed.to_csv(
                    handle,
                    sep="\t",
                    header=False,
                    index=False,
                    lineterminator="\n",
                )
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def write_summary_atomic(
    summary: pd.DataFrame,
    summary_path: Path,
) -> None:
    """Write the CSV summary atomically."""

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = summary_path.with_name(f"{summary_path.name}.tmp")

    try:
        summary.to_csv(temporary_path, index=False)
        temporary_path.replace(summary_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def discover_input_files(input_dir: Path) -> list[Path]:
    """Find significant-peak files and validate cell-type names."""

    if not input_dir.is_dir():
        raise FileNotFoundError(
            f"Unfiltered significant-peaks directory not found: {input_dir}"
        )

    input_paths = sorted(
        path
        for path in input_dir.glob(
            f"*{UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX}"
        )
        if path.is_file()
    )
    if not input_paths:
        raise FileNotFoundError(
            "No significant-peak files were found in "
            f"{input_dir}. Expected files named "
            f"<cell_type>{UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX}."
        )

    cell_types = [
        path.name[: -len(UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX)]
        for path in input_paths
    ]
    empty_cell_types = [
        path.name
        for path, cell_type in zip(input_paths, cell_types)
        if not cell_type
    ]
    if empty_cell_types:
        raise ValueError(
            f"Could not infer cell types from filenames: {empty_cell_types}"
        )
    if len(cell_types) != len(set(cell_types)):
        raise ValueError(
            "Multiple input files resolve to the same cell-type name."
        )

    return input_paths


def make_output_path(input_path: Path, output_dir: Path) -> tuple[str, Path]:
    """Infer the cell type and corresponding filtered output path."""

    cell_type = input_path.name[
        : -len(UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX)
    ]
    output_path = (
        output_dir
        / f"{cell_type}{FILTERED_SIGNIFICANT_PEAKS_SUFFIX}"
    )
    return cell_type, output_path


def validate_output_paths(
    input_paths: Sequence[Path],
    output_dir: Path,
    summary_path: Path,
    overwrite: bool,
) -> None:
    """Refuse accidental replacement unless ``--overwrite`` is supplied."""

    if overwrite:
        return

    existing_paths = [
        make_output_path(path, output_dir)[1]
        for path in input_paths
        if make_output_path(path, output_dir)[1].exists()
    ]
    if summary_path.exists():
        existing_paths.append(summary_path)

    if existing_paths:
        preview = "\n".join(f"  - {path}" for path in existing_paths[:10])
        suffix = (
            "\n  - ..."
            if len(existing_paths) > 10
            else ""
        )
        raise FileExistsError(
            "Output files already exist. Use --overwrite to replace them:\n"
            f"{preview}{suffix}"
        )


def process_peak_file(
    input_path: Path,
    output_dir: Path,
    eqtl_intervals: pd.DataFrame,
) -> dict[str, object]:
    """Remove complete peaks overlapping exact eQTL positions."""

    cell_type, output_path = make_output_path(input_path, output_dir)
    LOGGER.info("Processing %s", input_path)

    original_peaks = read_bed_preserving_values(input_path)
    n_input = len(original_peaks)

    if original_peaks.empty:
        filtered_original = original_peaks
    else:
        peak_intervals = build_bioframe_peak_table(original_peaks)

        # setdiff removes each whole interval from the first dataframe if it
        # overlaps any interval in the second dataframe. It does not split or
        # shorten peaks.
        filtered_intervals = bioframe.setdiff(
            peak_intervals,
            eqtl_intervals,
            cols1=(
                BIOFRAME_CHROM_COLUMN,
                BIOFRAME_START_COLUMN,
                BIOFRAME_END_COLUMN,
            ),
            cols2=(
                BIOFRAME_CHROM_COLUMN,
                BIOFRAME_START_COLUMN,
                BIOFRAME_END_COLUMN,
            ),
        )

        retained_row_ids = (
            filtered_intervals[BIOFRAME_ROW_ID_COLUMN]
            .astype("int64")
            .sort_values()
            .to_numpy()
        )
        filtered_original = original_peaks.iloc[retained_row_ids].copy()

    n_output = len(filtered_original)
    n_removed = n_input - n_output

    if n_input != n_removed + n_output:
        raise RuntimeError(
            f"Peak-count integrity check failed for {input_path}."
        )

    write_bed_atomic(filtered_original, output_path)

    fraction_removed = (
        float(n_removed / n_input)
        if n_input > 0
        else 0.0
    )
    LOGGER.info(
        "%s: retained %d/%d peaks; removed %d (%.2f%%)",
        cell_type,
        n_output,
        n_input,
        n_removed,
        100.0 * fraction_removed,
    )

    return {
        "cell_type": cell_type,
        "input_file": str(input_path),
        "output_file": str(output_path),
        "n_input_significant_peaks": n_input,
        "n_removed_eqtl_overlapping_peaks": n_removed,
        "n_significant_peaks_without_eqtl": n_output,
        "fraction_removed": fraction_removed,
    }


def run(overwrite: bool = False) -> pd.DataFrame:
    """Run exact-eQTL filtering for every significant-peak file."""

    input_paths = discover_input_files(UNFILTERED_SIG_PEAKS_DIR)
    FILTERED_SIG_PEAKS_DIR.mkdir(parents=True, exist_ok=True)

    validate_output_paths(
        input_paths=input_paths,
        output_dir=FILTERED_SIG_PEAKS_DIR,
        summary_path=SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH,
        overwrite=overwrite,
    )

    eqtl_intervals, _ = load_eqtl_intervals(
        LIVER_EQTL_SIGNIFICANT_PAIRS_PATH
    )

    records = [
        process_peak_file(
            input_path=input_path,
            output_dir=FILTERED_SIG_PEAKS_DIR,
            eqtl_intervals=eqtl_intervals,
        )
        for input_path in input_paths
    ]

    summary = pd.DataFrame.from_records(
        records,
        columns=SUMMARY_COLUMNS,
    ).sort_values("cell_type", kind="stable")

    write_summary_atomic(
        summary,
        SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH,
    )

    total_input = int(summary["n_input_significant_peaks"].sum())
    total_removed = int(
        summary["n_removed_eqtl_overlapping_peaks"].sum()
    )
    total_output = int(
        summary["n_significant_peaks_without_eqtl"].sum()
    )

    if total_input != total_removed + total_output:
        raise RuntimeError("Dataset-wide peak-count integrity check failed.")

    LOGGER.info(
        "Completed %d cell types: retained %d/%d peaks; removed %d "
        "(%.2f%%). Summary: %s",
        len(summary),
        total_output,
        total_input,
        total_removed,
        100.0 * total_removed / total_input if total_input else 0.0,
        SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH,
    )
    return summary


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Remove significant peaks that overlap exact GTEx liver "
            "eQTL positions."
        )
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing filtered BED files and summary.",
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        help="Logging verbosity. Default: INFO.",
    )
    return parser.parse_args()


def main() -> None:
    """Command-line entry point."""

    args = parse_arguments()
    configure_logging(args.log_level)

    try:
        run(overwrite=args.overwrite)
    except Exception:
        LOGGER.exception("Significant-peak eQTL filtering failed")
        raise


if __name__ == "__main__":
    main()
