"""Match snATAC matrix columns to the authors' cell annotations."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO

import pandas as pd

from constants import *

SNATAC_DIR = ATAC_SEQ_DIR

OUTPUT_DIR = MATCHING_ANNOTATED_COLUMNS_DIR


def load_annotations(path: Path) -> pd.DataFrame:
    """Load and validate the standardized Seurat annotations."""
    annotations = pd.read_csv(
        path,
        dtype={
            "sample": str,
            "raw_barcode": str,
            "seurat_cell_id": str,
            "cell_type_original": str,
            "cell_type": str,
            "condition": str,
        },
        low_memory=False,
    )

    required_columns = {
        "sample",
        "raw_barcode",
        "seurat_cell_id",
        "cell_type_original",
        "cell_type",
        "seurat_cluster",
        "condition",
        "is__cell_barcode",
    }

    missing_columns = required_columns - set(annotations.columns)

    if missing_columns:
        raise ValueError(
            "Annotation table is missing columns: "
            f"{sorted(missing_columns)}"
        )

    annotations["seurat_cluster"] = pd.to_numeric(
        annotations["seurat_cluster"],
        errors="raise",
    ).astype(int)

    annotations["is__cell_barcode"] = pd.to_numeric(
        annotations["is__cell_barcode"],
        errors="raise",
    ).astype(int)

    if annotations["seurat_cell_id"].duplicated().any():
        raise ValueError("seurat_cell_id is not globally unique.")

    if annotations.duplicated(["sample", "raw_barcode"]).any():
        raise ValueError(
            "The combination of sample and raw_barcode is not unique."
        )

    return annotations


def load_sample_mapping(path: Path) -> pd.DataFrame:
    """Load the one-to-one mapping between Seurat samples and GSM folders."""
    mapping = pd.read_csv(
        path,
        usecols=["seurat_sample", "sample_directory"],
        dtype=str,
    ).drop_duplicates()

    if mapping["seurat_sample"].duplicated().any():
        raise ValueError(
            "A Seurat sample maps to more than one GSM directory."
        )

    if mapping["sample_directory"].duplicated().any():
        raise ValueError(
            "A GSM directory maps to more than one Seurat sample."
        )

    return mapping.sort_values("sample_directory").reset_index(drop=True)


def find_single_file(directory: Path, pattern: str) -> Path:
    """Return the only file in a directory matching a glob pattern."""
    matches = sorted(directory.glob(pattern))

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one '{pattern}' file in {directory}; "
            f"found {len(matches)}."
        )

    return matches[0]


def load_barcodes(path: Path) -> pd.DataFrame:
    """Read matrix barcodes while preserving their original column order."""
    barcodes = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["raw_barcode"],
        dtype=str,
        compression="gzip",
    )

    if barcodes["raw_barcode"].duplicated().any():
        raise ValueError(f"Duplicate barcodes found in {path}.")

    barcodes.insert(
        loc=0,
        column="matrix_column_index_zero_based",
        value=range(len(barcodes)),
    )

    barcodes.insert(
        loc=1,
        column="matrix_column_index_one_based",
        value=range(1, len(barcodes) + 1),
    )

    return barcodes


def open_text_file(path: Path) -> TextIO:
    """Open a plain-text or gzip-compressed text file."""
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")

    return path.open(mode="rt", encoding="utf-8")


def read_matrix_market_shape(path: Path) -> tuple[int, int, int]:
    """Read the dimensions of a Matrix Market file without loading it."""
    with open_text_file(path) as handle:
        header = handle.readline().strip()

        if not header.startswith("%%MatrixMarket"):
            raise ValueError(
                f"{path} does not appear to be a Matrix Market file."
            )

        for line in handle:
            stripped = line.strip()

            if not stripped or stripped.startswith("%"):
                continue

            values = stripped.split()

            if len(values) != 3:
                raise ValueError(
                    f"Invalid Matrix Market dimensions line in {path}: "
                    f"{stripped}"
                )

            n_rows, n_columns, n_nonzero = map(int, values)
            return n_rows, n_columns, n_nonzero

    raise ValueError(f"No dimensions line found in {path}.")


def match_sample_columns(
    annotations: pd.DataFrame,
    seurat_sample: str,
    sample_directory_name: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Match one matrix's ordered barcode columns to Seurat annotations."""
    sample_directory = SNATAC_DIR / sample_directory_name

    if not sample_directory.is_dir():
        raise FileNotFoundError(
            f"Sample directory does not exist: {sample_directory}"
        )

    barcode_path = find_single_file(
        sample_directory,
        "*_barcodes.tsv.gz",
    )
    matrix_path = find_single_file(
        sample_directory,
        "*_matrix.mtx.gz",
    )

    barcodes = load_barcodes(barcode_path)

    n_matrix_rows, n_matrix_columns, n_nonzero = (
        read_matrix_market_shape(matrix_path)
    )

    if n_matrix_columns != len(barcodes):
        raise ValueError(
            f"{sample_directory_name}: matrix has "
            f"{n_matrix_columns:,} columns but barcode file has "
            f"{len(barcodes):,} rows."
        )

    sample_annotations = annotations.loc[
        annotations["sample"] == seurat_sample,
        [
            "raw_barcode",
            "seurat_cell_id",
            "cell_type_original",
            "cell_type",
            "seurat_cluster",
            "condition",
            "is__cell_barcode",
        ],
    ].copy()

    if sample_annotations.empty:
        raise ValueError(
            f"No annotations found for Seurat sample {seurat_sample}."
        )

    # This is an inner join: only matrix columns with an annotation survive.
    matched = barcodes.merge(
        sample_annotations,
        on="raw_barcode",
        how="inner",
        validate="one_to_one",
        sort=False,
    )

    matched = matched.sort_values(
        "matrix_column_index_zero_based"
    ).reset_index(drop=True)

    matched.insert(0, "sample_directory", sample_directory_name)
    matched.insert(1, "seurat_sample", seurat_sample)

    if matched["matrix_column_index_zero_based"].duplicated().any():
        raise ValueError(
            f"{sample_directory_name}: a matrix column matched twice."
        )

    if matched["seurat_cell_id"].duplicated().any():
        raise ValueError(
            f"{sample_directory_name}: a Seurat cell matched twice."
        )

    n_flagged_annotations = int(
        (sample_annotations["is__cell_barcode"] == 1).sum()
    )

    n_matched_flagged = int(
        (matched["is__cell_barcode"] == 1).sum()
    )

    n_matched_unflagged = int(
        (matched["is__cell_barcode"] == 0).sum()
    )

    summary = {
        "sample_directory": sample_directory_name,
        "seurat_sample": seurat_sample,
        "matrix_rows": n_matrix_rows,
        "matrix_columns": n_matrix_columns,
        "matrix_nonzero_entries": n_nonzero,
        "annotations_total": len(sample_annotations),
        "annotations_is_cell_barcode_1": n_flagged_annotations,
        "matched_columns": len(matched),
        "matched_is_cell_barcode_1": n_matched_flagged,
        "matched_is_cell_barcode_0": n_matched_unflagged,
        "unmatched_matrix_columns": n_matrix_columns - len(matched),
        "fraction_matrix_columns_matched": (
            len(matched) / n_matrix_columns
            if n_matrix_columns > 0
            else 0.0
        ),
        "fraction_annotations_matched": (
            len(matched) / len(sample_annotations)
            if len(sample_annotations) > 0
            else 0.0
        ),
    }

    return matched, summary


def main() -> None:
    """Match columns for all 12 snATAC matrices."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    annotations = load_annotations(CELL_ANNOTATIONS_PATH)
    sample_mapping = load_sample_mapping(SAMPLE_MAPPING_PATH)

    matched_tables: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []

    for row in sample_mapping.itertuples(index=False):
        matched, summary = match_sample_columns(
            annotations=annotations,
            seurat_sample=row.seurat_sample,
            sample_directory_name=row.sample_directory,
        )

        sample_output_directory = OUTPUT_DIR / row.sample_directory
        sample_output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        sample_output_path = (
            sample_output_directory
            / "matched_columns.csv.gz"
        )

        matched.to_csv(
            sample_output_path,
            index=False,
            compression="gzip",
        )

        matched_tables.append(matched)
        summaries.append(summary)

        print(
            f"{row.sample_directory}: "
            f"{len(matched):,} matched columns"
        )

    all_matched = pd.concat(
        matched_tables,
        ignore_index=True,
    )

    summary_table = pd.DataFrame.from_records(summaries)

    cell_type_counts = (
        all_matched.groupby(
            ["sample_directory", "seurat_sample", "cell_type"],
            observed=True,
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )

    all_matched.to_csv(
        OUTPUT_DIR / "all_matched_columns.csv.gz",
        index=False,
        compression="gzip",
    )

    summary_table.to_csv(
        OUTPUT_DIR / "match_summary.csv",
        index=False,
    )

    cell_type_counts.to_csv(
        OUTPUT_DIR / "matched_cell_type_counts.csv",
        index=False,
    )

    print("\nMatching summary:")
    print(
        summary_table[
            [
                "sample_directory",
                "seurat_sample",
                "matrix_columns",
                "annotations_total",
                "annotations_is_cell_barcode_1",
                "matched_columns",
                "matched_is_cell_barcode_1",
                "matched_is_cell_barcode_0",
            ]
        ].to_string(index=False)
    )

    print("\nTotal matched columns:")
    print(f"{len(all_matched):,}")

    print("\nMatched cells by cell type:")
    print(
        all_matched["cell_type"]
        .value_counts()
        .to_string()
    )

    expected_total = 55_014

    if len(all_matched) != expected_total:
        print(
            "\nWarning: based on the previous overlap results, "
            f"we expected {expected_total:,} matched columns, "
            f"but found {len(all_matched):,}."
        )


if __name__ == "__main__":
    main()