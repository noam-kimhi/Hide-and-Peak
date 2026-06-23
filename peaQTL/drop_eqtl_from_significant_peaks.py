"""Remove significant peaks that overlap exact GTEx liver eQTL positions.

For every file in ``UNFILTERED_SIG_PEAKS_DIR`` named

    <cell_type>_significant_peaks.bed.gz

the script writes

    <cell_type>_significant_peaks_without_eqtl.bed.gz

to ``FILTERED_SIG_PEAKS_DIR`` and also creates:

1. One pie chart per cell type showing the fraction of significant peaks that
   overlap at least one exact liver eQTL position versus the fraction that do
   not.
2. A global bar plot with one bar per cell type plus one bar for the liver
   eQTL set size.
3. A second global bar plot with the same bars, plus a horizontal line inside
   each cell-type bar indicating the number of eQTL-overlapping significant
   peaks for that cell type.

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
import matplotlib.pyplot as plt
import pandas as pd

from constants import (
    FILTERED_SIG_PEAKS_DIR,
    FILTERED_SIGNIFICANT_PEAKS_SUFFIX,
    LIVER_EQTL_SIGNIFICANT_PAIRS_PATH,
    PLOT_DPI,
    SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH,
    SIGNIFICANT_PEAKS_EQTL_PLOTS_DIR,
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

PIE_CHART_SUFFIX: Final[str] = "_significant_peaks_eqtl_pie.png"
TOTAL_COUNTS_BAR_PLOT_FILENAME: Final[str] = (
    "significant_peaks_total_counts_and_eqtl_bar_plot.png"
)
TOTAL_COUNTS_WITH_OVERLAP_FILENAME: Final[str] = (
    "significant_peaks_total_counts_with_eqtl_overlap_lines.png"
)

SUMMARY_COLUMNS: Final[Sequence[str]] = (
    "cell_type",
    "input_file",
    "output_file",
    "pie_chart_file",
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
    """Convert common chromosome spellings to a canonical ``chr*`` form."""

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
    return f"chr{suffix}"


def parse_gtex_variant_ids(variant_ids: pd.Series) -> pd.DataFrame:
    """Parse GTEx variant IDs and return unique one-base GRCh38 intervals."""

    ids = variant_ids.astype("string").str.strip()
    if ids.isna().any() or (ids == "").any():
        raise ValueError("The GTEx eQTL table contains empty variant IDs.")

    parsed = ids.str.rsplit("_", n=4, expand=True)
    if parsed.shape[1] != 5:
        examples = ids.head(5).tolist()
        raise ValueError(
            "Could not parse GTEx variant IDs into chromosome, position, "
            "reference, alternate, and build fields. "
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
            BIOFRAME_CHROM_COLUMN: parsed["chrom_raw"].map(normalize_chromosome),
            BIOFRAME_START_COLUMN: positions - 1,
            BIOFRAME_END_COLUMN: positions,
        }
    )

    intervals = intervals.drop_duplicates(
        subset=(BIOFRAME_CHROM_COLUMN, BIOFRAME_START_COLUMN, BIOFRAME_END_COLUMN),
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
    return eqtl_intervals, len(eqtl_intervals)


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

    starts = pd.to_numeric(bed.iloc[:, START_COLUMN_INDEX], errors="coerce")
    ends = pd.to_numeric(bed.iloc[:, END_COLUMN_INDEX], errors="coerce")

    invalid_coordinate_mask = starts.isna() | ends.isna() | (starts < 0) | (ends <= starts)
    if invalid_coordinate_mask.any():
        invalid_rows = invalid_coordinate_mask[invalid_coordinate_mask].index[:10].tolist()
        raise ValueError(
            f"{bed_path} contains invalid BED coordinates at zero-based row "
            f"indices {invalid_rows}."
        )

    return bed


def build_bioframe_peak_table(peaks: pd.DataFrame) -> pd.DataFrame:
    """Create a normalized interval table linked to original BED rows."""

    return pd.DataFrame(
        {
            BIOFRAME_CHROM_COLUMN: peaks.iloc[:, CHROM_COLUMN_INDEX].map(normalize_chromosome),
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


def write_summary_atomic(summary: pd.DataFrame, summary_path: Path) -> None:
    """Write the CSV summary atomically."""

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = summary_path.with_name(f"{summary_path.name}.tmp")

    try:
        summary.to_csv(temporary_path, index=False)
        temporary_path.replace(summary_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def make_temporary_plot_path(output_path: Path) -> Path:
    """Create a temporary plot path that preserves the real image suffix."""

    return output_path.with_name(f"{output_path.stem}.tmp{output_path.suffix}")


def save_figure_atomic(figure: plt.Figure, output_path: Path) -> None:
    """Save a Matplotlib figure atomically as PNG."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = make_temporary_plot_path(output_path)

    try:
        figure.savefig(
            temporary_path,
            dpi=PLOT_DPI,
            bbox_inches="tight",
            format="png",
        )
        temporary_path.replace(output_path)
    finally:
        plt.close(figure)
        if temporary_path.exists():
            temporary_path.unlink()


def create_pie_chart(
    pie_chart_path: Path,
    *,
    cell_type: str,
    n_eqtl_overlap: int,
    n_without_eqtl: int,
) -> None:
    """Create a pie chart summarizing eQTL overlap for one cell type."""

    figure = plt.figure(figsize=(6, 6))
    if n_eqtl_overlap + n_without_eqtl == 0:
        plt.text(0.5, 0.5, "No significant peaks", ha="center", va="center", fontsize=14)
        plt.axis("off")
    else:
        plt.pie(
            [n_eqtl_overlap, n_without_eqtl],
            labels=["Overlaps eQTL", "No eQTL overlap"],
            autopct="%1.1f%%",
            startangle=90,
        )
        plt.axis("equal")

    plt.title(f"{cell_type}: significant peaks vs. liver eQTL overlap")
    plt.tight_layout()
    save_figure_atomic(figure, pie_chart_path)


def add_bar_value_labels(ax: plt.Axes, heights: Sequence[int]) -> None:
    """Add value labels above bars."""

    upper_limit = max(heights) if heights else 0
    offset = max(upper_limit * 0.01, 1.0)
    for x_position, height in enumerate(heights):
        ax.text(
            x_position,
            height + offset,
            f"{int(height)}",
            ha="center",
            va="bottom",
            fontsize=9,
            rotation=0,
        )


def create_total_counts_bar_plot(
    summary: pd.DataFrame,
    *,
    n_unique_eqtl_positions: int,
    output_path: Path,
) -> None:
    """Create the global total-count bar plot requested by the user."""

    ordered = summary.sort_values("cell_type", kind="stable").reset_index(drop=True)
    categories = ordered["cell_type"].tolist() + ["eQTL"]
    heights = ordered["n_input_significant_peaks"].astype(int).tolist() + [
        int(n_unique_eqtl_positions)
    ]

    figure, ax = plt.subplots(figsize=(max(8, len(categories) * 1.1), 6))
    ax.bar(range(len(categories)), heights)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=45, ha="right")
    ax.set_ylabel("Total interval count")
    ax.set_xlabel("Cell type / reference set")
    ax.set_title("Total significant-peak counts by cell type and liver eQTL set size")
    add_bar_value_labels(ax, heights)
    figure.tight_layout()
    save_figure_atomic(figure, output_path)


def create_total_counts_with_overlap_lines_plot(
    summary: pd.DataFrame,
    *,
    n_unique_eqtl_positions: int,
    output_path: Path,
) -> None:
    """Create the second global bar plot with internal eQTL-overlap lines."""

    ordered = summary.sort_values("cell_type", kind="stable").reset_index(drop=True)
    categories = ordered["cell_type"].tolist() + ["eQTL"]
    total_heights = ordered["n_input_significant_peaks"].astype(int).tolist() + [
        int(n_unique_eqtl_positions)
    ]
    overlap_counts = ordered["n_removed_eqtl_overlapping_peaks"].astype(int).tolist()

    figure, ax = plt.subplots(figsize=(max(8, len(categories) * 1.1), 6))
    bars = ax.bar(range(len(categories)), total_heights)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=45, ha="right")
    ax.set_ylabel("Total interval count")
    ax.set_xlabel("Cell type / reference set")
    ax.set_title(
        "Total significant-peak counts with eQTL-overlapping peak counts per cell type"
    )

    add_bar_value_labels(ax, total_heights)

    added_legend_label = False
    for index, overlap_count in enumerate(overlap_counts):
        bar = bars[index]
        left_x = bar.get_x()
        right_x = left_x + bar.get_width()
        label = "eQTL-overlapping peaks" if not added_legend_label else None
        ax.hlines(
            y=overlap_count,
            xmin=left_x,
            xmax=right_x,
            linewidth=2.0,
            label=label,
        )
        added_legend_label = True

        vertical_offset = max(max(total_heights) * 0.01, 1.0)
        text_y = min(overlap_count + vertical_offset, bar.get_height() - vertical_offset)
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            text_y,
            str(int(overlap_count)),
            ha="center",
            va="bottom" if text_y >= overlap_count else "top",
            fontsize=9,
        )

    if added_legend_label:
        ax.legend()

    figure.tight_layout()
    save_figure_atomic(figure, output_path)


def discover_input_files(input_dir: Path) -> list[Path]:
    """Find significant-peak files and validate cell-type names."""

    if not input_dir.is_dir():
        raise FileNotFoundError(
            f"Unfiltered significant-peaks directory not found: {input_dir}"
        )

    input_paths = sorted(
        path
        for path in input_dir.glob(f"*{UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX}")
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
        path.name for path, cell_type in zip(input_paths, cell_types) if not cell_type
    ]
    if empty_cell_types:
        raise ValueError(f"Could not infer cell types from filenames: {empty_cell_types}")
    if len(cell_types) != len(set(cell_types)):
        raise ValueError("Multiple input files resolve to the same cell-type name.")
    return input_paths


def make_output_path(input_path: Path, output_dir: Path) -> tuple[str, Path]:
    """Infer the cell type and corresponding filtered output path."""

    cell_type = input_path.name[: -len(UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX)]
    output_path = output_dir / f"{cell_type}{FILTERED_SIGNIFICANT_PEAKS_SUFFIX}"
    return cell_type, output_path


def make_pie_chart_path(cell_type: str, plots_dir: Path) -> Path:
    """Construct the pie-chart path for one cell type."""

    return plots_dir / f"{cell_type}{PIE_CHART_SUFFIX}"


def get_global_plot_paths(plots_dir: Path) -> tuple[Path, Path]:
    """Return the two requested global plot paths."""

    return (
        plots_dir / TOTAL_COUNTS_BAR_PLOT_FILENAME,
        plots_dir / TOTAL_COUNTS_WITH_OVERLAP_FILENAME,
    )


def validate_output_paths(
    input_paths: Sequence[Path],
    output_dir: Path,
    plots_dir: Path,
    summary_path: Path,
    overwrite: bool,
) -> None:
    """Refuse accidental replacement unless ``--overwrite`` is supplied."""

    if overwrite:
        return

    existing_paths: list[Path] = []
    for path in input_paths:
        cell_type, output_path = make_output_path(path, output_dir)
        pie_chart_path = make_pie_chart_path(cell_type, plots_dir)
        if output_path.exists():
            existing_paths.append(output_path)
        if pie_chart_path.exists():
            existing_paths.append(pie_chart_path)

    total_counts_plot_path, overlap_lines_plot_path = get_global_plot_paths(plots_dir)
    if total_counts_plot_path.exists():
        existing_paths.append(total_counts_plot_path)
    if overlap_lines_plot_path.exists():
        existing_paths.append(overlap_lines_plot_path)
    if summary_path.exists():
        existing_paths.append(summary_path)

    if existing_paths:
        preview = "\n".join(f"  - {path}" for path in existing_paths[:10])
        suffix = "\n  - ..." if len(existing_paths) > 10 else ""
        raise FileExistsError(
            "Output files already exist. Use --overwrite to replace them:\n"
            f"{preview}{suffix}"
        )


def process_peak_file(
    input_path: Path,
    output_dir: Path,
    plots_dir: Path,
    eqtl_intervals: pd.DataFrame,
) -> dict[str, object]:
    """Remove complete peaks overlapping exact eQTL positions."""

    cell_type, output_path = make_output_path(input_path, output_dir)
    pie_chart_path = make_pie_chart_path(cell_type, plots_dir)
    LOGGER.info("Processing %s", input_path)

    original_peaks = read_bed_preserving_values(input_path)
    n_input = len(original_peaks)

    if original_peaks.empty:
        filtered_original = original_peaks
    else:
        peak_intervals = build_bioframe_peak_table(original_peaks)
        filtered_intervals = bioframe.setdiff(
            peak_intervals,
            eqtl_intervals,
            cols1=(BIOFRAME_CHROM_COLUMN, BIOFRAME_START_COLUMN, BIOFRAME_END_COLUMN),
            cols2=(BIOFRAME_CHROM_COLUMN, BIOFRAME_START_COLUMN, BIOFRAME_END_COLUMN),
        )
        retained_row_ids = (
            filtered_intervals[BIOFRAME_ROW_ID_COLUMN].astype("int64").sort_values().to_numpy()
        )
        filtered_original = original_peaks.iloc[retained_row_ids].copy()

    n_output = len(filtered_original)
    n_removed = n_input - n_output

    if n_input != n_removed + n_output:
        raise RuntimeError(f"Peak-count integrity check failed for {input_path}.")

    write_bed_atomic(filtered_original, output_path)
    create_pie_chart(
        pie_chart_path,
        cell_type=cell_type,
        n_eqtl_overlap=n_removed,
        n_without_eqtl=n_output,
    )

    fraction_removed = float(n_removed / n_input) if n_input > 0 else 0.0
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
        "pie_chart_file": str(pie_chart_path),
        "n_input_significant_peaks": n_input,
        "n_removed_eqtl_overlapping_peaks": n_removed,
        "n_significant_peaks_without_eqtl": n_output,
        "fraction_removed": fraction_removed,
    }


def run(overwrite: bool = False) -> pd.DataFrame:
    """Run exact-eQTL filtering for every significant-peak file."""

    input_paths = discover_input_files(UNFILTERED_SIG_PEAKS_DIR)
    FILTERED_SIG_PEAKS_DIR.mkdir(parents=True, exist_ok=True)
    SIGNIFICANT_PEAKS_EQTL_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    validate_output_paths(
        input_paths=input_paths,
        output_dir=FILTERED_SIG_PEAKS_DIR,
        plots_dir=SIGNIFICANT_PEAKS_EQTL_PLOTS_DIR,
        summary_path=SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH,
        overwrite=overwrite,
    )

    eqtl_intervals, n_unique_eqtl_positions = load_eqtl_intervals(
        LIVER_EQTL_SIGNIFICANT_PAIRS_PATH
    )

    records = [
        process_peak_file(
            input_path=input_path,
            output_dir=FILTERED_SIG_PEAKS_DIR,
            plots_dir=SIGNIFICANT_PEAKS_EQTL_PLOTS_DIR,
            eqtl_intervals=eqtl_intervals,
        )
        for input_path in input_paths
    ]

    summary = pd.DataFrame.from_records(records, columns=SUMMARY_COLUMNS).sort_values(
        "cell_type",
        kind="stable",
    )
    write_summary_atomic(summary, SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH)

    total_input = int(summary["n_input_significant_peaks"].sum())
    total_removed = int(summary["n_removed_eqtl_overlapping_peaks"].sum())
    total_output = int(summary["n_significant_peaks_without_eqtl"].sum())

    if total_input != total_removed + total_output:
        raise RuntimeError("Dataset-wide peak-count integrity check failed.")

    total_counts_plot_path, overlap_lines_plot_path = get_global_plot_paths(
        SIGNIFICANT_PEAKS_EQTL_PLOTS_DIR
    )
    create_total_counts_bar_plot(
        summary,
        n_unique_eqtl_positions=n_unique_eqtl_positions,
        output_path=total_counts_plot_path,
    )
    create_total_counts_with_overlap_lines_plot(
        summary,
        n_unique_eqtl_positions=n_unique_eqtl_positions,
        output_path=overlap_lines_plot_path,
    )

    LOGGER.info(
        "Completed %d cell types: retained %d/%d peaks; removed %d (%.2f%%). "
        "Summary: %s | Global plots: %s ; %s",
        len(summary),
        total_output,
        total_input,
        total_removed,
        100.0 * total_removed / total_input if total_input else 0.0,
        SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH,
        total_counts_plot_path,
        overlap_lines_plot_path,
    )
    return summary


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Remove significant peaks that overlap exact GTEx liver eQTL "
            "positions and create pie charts plus two summary bar plots."
        )
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing filtered BED files, plots, and summary.",
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
