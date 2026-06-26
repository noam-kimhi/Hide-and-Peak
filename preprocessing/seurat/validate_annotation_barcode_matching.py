"""Validate Seurat annotations against the downloaded snATAC barcode files."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.optimize import linear_sum_assignment

from constants import *


def load_and_clean_annotations(path: Path) -> pd.DataFrame:
    """Load and standardize the per-nucleus Seurat metadata."""
    annotations = pd.read_csv(path, low_memory=False)

    required_columns = {
        "sample",
        "cells",
        "seurat_cell_id",
        "seurat_cluster",
        "cell_type_original",
        "active_ident",
    }
    missing_columns = required_columns - set(annotations.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required annotation columns: {sorted(missing_columns)}"
        )

    if not annotations["cell_type_original"].equals(annotations["active_ident"]):
        raise ValueError("cell_type_original and active_ident are not identical.")

    if not annotations["seurat_cell_id"].is_unique:
        raise ValueError("seurat_cell_id is not unique.")

    if not annotations["cells"].equals(annotations["seurat_cell_id"]):
        raise ValueError("cells and seurat_cell_id are not identical.")

    annotations["raw_barcode"] = (
        annotations["seurat_cell_id"]
        .astype(str)
        .str.extract(r"([ACGT]+-\d+)$", expand=False)
    )

    if annotations["raw_barcode"].isna().any():
        n_missing = int(annotations["raw_barcode"].isna().sum())
        raise ValueError(
            f"Could not extract a barcode from {n_missing} cell identifiers."
        )

    expected_cell_id = (
        annotations["sample"].astype(str)
        + "_"
        + annotations["raw_barcode"]
    )

    if not expected_cell_id.equals(annotations["seurat_cell_id"]):
        raise ValueError(
            "Some Seurat cell identifiers do not equal sample + raw barcode."
        )

    if annotations.duplicated(["sample", "raw_barcode"]).any():
        raise ValueError(
            "The combination of sample and raw_barcode is not unique."
        )

    annotations["cell_type"] = annotations[
        "cell_type_original"
    ].map(CELL_TYPE_STANDARDIZATION)

    if annotations["cell_type"].isna().any():
        unknown_values = sorted(
            annotations.loc[
                annotations["cell_type"].isna(),
                "cell_type_original",
            ].unique()
        )
        raise ValueError(
            f"Unrecognized cell-type labels: {unknown_values}"
        )

    annotations["condition"] = annotations["sample"].map(
        lambda sample: (
            "MASH"
            if str(sample).startswith("NASH")
            else "Normal"
        )
    )

    annotations["seurat_cluster"] = pd.to_numeric(
        annotations["seurat_cluster"],
        errors="raise",
    ).astype(int)

    return annotations


def find_barcode_files(snatac_dir: Path) -> dict[str, Path]:
    """Find the barcode file belonging to every downloaded GSM directory."""
    barcode_files: dict[str, Path] = {}

    for sample_dir in sorted(snatac_dir.glob("GSM*")):
        if not sample_dir.is_dir():
            continue

        matches = list(sample_dir.glob("*_barcodes.tsv.gz"))

        if len(matches) != 1:
            raise ValueError(
                f"Expected one barcode file in {sample_dir}, "
                f"found {len(matches)}."
            )

        barcode_files[sample_dir.name] = matches[0]

    if not barcode_files:
        raise FileNotFoundError(
            f"No GSM barcode files found beneath {snatac_dir}."
        )

    return barcode_files


def load_barcode_set(path: Path) -> set[str]:
    """Read a one-column 10x barcode file into a set."""
    barcodes = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["raw_barcode"],
        dtype=str,
        compression="gzip",
    )["raw_barcode"]

    if barcodes.duplicated().any():
        raise ValueError(f"Duplicate barcodes found in {path}.")

    return set(barcodes)


def calculate_overlap_table(
    annotations: pd.DataFrame,
    barcode_files: dict[str, Path],
) -> pd.DataFrame:
    """Calculate barcode overlap for every Seurat sample–GSM pair."""
    annotation_sets = {
        sample: set(group["raw_barcode"])
        for sample, group in annotations.groupby("sample")
    }

    downloaded_sets = {
        directory: load_barcode_set(path)
        for directory, path in barcode_files.items()
    }

    records: list[dict[str, object]] = []

    for seurat_sample, annotation_barcodes in annotation_sets.items():
        for sample_directory, downloaded_barcodes in downloaded_sets.items():
            overlap = len(
                annotation_barcodes & downloaded_barcodes
            )

            records.append(
                {
                    "seurat_sample": seurat_sample,
                    "sample_directory": sample_directory,
                    "annotation_cells": len(annotation_barcodes),
                    "downloaded_barcodes": len(downloaded_barcodes),
                    "overlap": overlap,
                    "fraction_annotations_matched": (
                        overlap / len(annotation_barcodes)
                    ),
                    "fraction_downloaded_matched": (
                        overlap / len(downloaded_barcodes)
                    ),
                }
            )

    return pd.DataFrame.from_records(records)


def infer_one_to_one_sample_mapping(
    overlap_table: pd.DataFrame,
) -> pd.DataFrame:
    """Infer the optimal one-to-one mapping by maximizing barcode overlap."""
    score_matrix = overlap_table.pivot(
        index="seurat_sample",
        columns="sample_directory",
        values="overlap",
    )

    row_indices, column_indices = linear_sum_assignment(
        -score_matrix.to_numpy()
    )

    mapping = pd.DataFrame(
        {
            "seurat_sample": score_matrix.index[row_indices],
            "sample_directory": score_matrix.columns[column_indices],
        }
    )

    mapping = mapping.merge(
        overlap_table,
        on=["seurat_sample", "sample_directory"],
        how="left",
        validate="one_to_one",
    )

    return mapping.sort_values("sample_directory").reset_index(drop=True)


def main() -> None:
    """Run all annotation and sample-matching validation steps."""
    annotations = load_and_clean_annotations(CELL_ANNOTATIONS_PATH)
    barcode_files = find_barcode_files(ATAC_SEQ_DIR)

    overlap_table = calculate_overlap_table(
        annotations=annotations,
        barcode_files=barcode_files,
    )

    mapping = infer_one_to_one_sample_mapping(overlap_table)

    annotation_output_path = (
        ANNOTATION_DIR / "GSE281367_cell_annotations.csv.gz"
    )
    overlap_output_path = ANNOTATION_DIR / "sample_barcode_overlaps.csv"
    mapping_output_path = ANNOTATION_DIR / "sample_directory_mapping.csv"

    annotations.to_csv(
        annotation_output_path,
        index=False,
        compression="gzip",
    )
    overlap_table.to_csv(overlap_output_path, index=False)
    mapping.to_csv(mapping_output_path, index=False)

    print(f"Number of annotated nuclei: {len(annotations):,}")
    print(f"Unique Seurat samples: {annotations['sample'].nunique()}")
    print(f"Downloaded GSM directories: {len(barcode_files)}")

    print("\nInferred mapping:")
    print(
        mapping[
            [
                "seurat_sample",
                "sample_directory",
                "annotation_cells",
                "downloaded_barcodes",
                "overlap",
                "fraction_annotations_matched",
                "fraction_downloaded_matched"
            ]
        ].to_string(index=False)
    )

    low_match = mapping[
        "fraction_annotations_matched"
    ] < 0.95

    if low_match.any():
        print(
            "At least one inferred sample mapping matched fewer than "
            "95% of the annotated cells."
        )

    print("\nCell-type counts:")
    print(
        annotations["cell_type"]
        .value_counts()
        .to_string()
    )

    print("\nCluster-to-cell-type table:")
    print(
        pd.crosstab(
            annotations["seurat_cluster"],
            annotations["cell_type"],
        ).to_string()
    )

    print(f"\nWrote annotations to: {annotation_output_path}")
    print(f"Wrote overlap table to: {overlap_output_path}")
    print(f"Wrote sample mapping to: {mapping_output_path}")


if __name__ == "__main__":
    main()