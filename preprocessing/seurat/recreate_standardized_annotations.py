"""Recreate the standardized GSE281367 per-cell annotation table."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from constants import *

SOURCE_METADATA_PATH = SEURAT_METADATA_PATH

OUTPUT_DIR = ANNOTATION_DIR

OUTPUT_ANNOTATIONS_PATH = CELL_ANNOTATIONS_PATH


def load_source_metadata(path: Path) -> pd.DataFrame:
    """Load the metadata exported directly from the Seurat object."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Could not find the extracted Seurat metadata: {path}"
        )

    metadata = pd.read_csv(
        path,
        low_memory=False,
        dtype={
            "sample": str,
            "cells": str,
            "CellType.2": str,
            "seurat_cell_id": str,
            "active_ident": str,
        },
    )

    required_columns = {
        "sample",
        "cells",
        "seurat_clusters",
        "CellType.2",
        "seurat_cell_id",
        "active_ident",
        "is__cell_barcode",
    }

    missing_columns = required_columns - set(metadata.columns)

    if missing_columns:
        raise ValueError(
            "The source metadata is missing columns: "
            f"{sorted(missing_columns)}"
        )

    return metadata


def validate_original_metadata(metadata: pd.DataFrame) -> None:
    """Validate key relationships in the exported Seurat metadata."""
    if metadata["seurat_cell_id"].duplicated().any():
        raise ValueError("seurat_cell_id is not globally unique.")

    if not metadata["cells"].equals(metadata["seurat_cell_id"]):
        raise ValueError(
            "The cells and seurat_cell_id columns are not identical."
        )

    if not metadata["CellType.2"].equals(metadata["active_ident"]):
        raise ValueError(
            "CellType.2 and active_ident are not identical."
        )


def standardize_annotations(metadata: pd.DataFrame) -> pd.DataFrame:
    """Create the derived columns required by downstream processing."""
    annotations = metadata.copy()

    annotations["raw_barcode"] = (
        annotations["seurat_cell_id"]
        .str.extract(r"([ACGT]+-\d+)$", expand=False)
    )

    if annotations["raw_barcode"].isna().any():
        n_missing = int(annotations["raw_barcode"].isna().sum())
        raise ValueError(
            f"Could not extract raw barcodes from {n_missing} cells."
        )

    expected_cell_ids = (
        annotations["sample"]
        + "_"
        + annotations["raw_barcode"]
    )

    invalid_ids = (
        expected_cell_ids
        != annotations["seurat_cell_id"]
    )

    if invalid_ids.any():
        raise ValueError(
            f"{int(invalid_ids.sum())} cell identifiers do not equal "
            "sample + '_' + raw_barcode."
        )

    annotations = annotations.rename(
        columns={
            "CellType.2": "cell_type_original",
            "seurat_clusters": "seurat_cluster",
        }
    )

    annotations["cell_type"] = annotations[
        "cell_type_original"
    ].map(CELL_TYPE_STANDARDIZATION)

    if annotations["cell_type"].isna().any():
        unrecognized = sorted(
            annotations.loc[
                annotations["cell_type"].isna(),
                "cell_type_original",
            ]
            .dropna()
            .unique()
        )

        raise ValueError(
            "Unrecognized original cell-type labels: "
            f"{unrecognized}"
        )

    annotations["condition"] = annotations["sample"].map(
        lambda sample: (
            "MASH"
            if sample.startswith("NASH")
            else "Normal"
        )
    )

    annotations["seurat_cluster"] = pd.to_numeric(
        annotations["seurat_cluster"],
        errors="raise",
    ).astype(int)

    annotations["is__cell_barcode"] = pd.to_numeric(
        annotations["is__cell_barcode"],
        errors="raise",
    ).astype(int)

    if annotations.duplicated(
        ["sample", "raw_barcode"]
    ).any():
        raise ValueError(
            "The combination of sample and raw_barcode is not unique."
        )

    return annotations


def validate_standardized_annotations(
    annotations: pd.DataFrame,
) -> None:
    """Check the expected structure of the recreated annotation table."""
    expected_nuclei = 69_601
    expected_samples = 12
    expected_clusters = 15

    if len(annotations) != expected_nuclei:
        raise ValueError(
            f"Expected {expected_nuclei:,} nuclei, "
            f"found {len(annotations):,}."
        )

    if annotations["sample"].nunique() != expected_samples:
        raise ValueError(
            f"Expected {expected_samples} samples, found "
            f"{annotations['sample'].nunique()}."
        )

    if annotations["seurat_cluster"].nunique() != expected_clusters:
        raise ValueError(
            f"Expected {expected_clusters} clusters, found "
            f"{annotations['seurat_cluster'].nunique()}."
        )

    observed_original_types = set(
        annotations["cell_type_original"].unique()
    )

    expected_original_types = set(
        CELL_TYPE_STANDARDIZATION
    )

    if observed_original_types != expected_original_types:
        raise ValueError(
            "Unexpected original cell-type labels. "
            f"Observed: {sorted(observed_original_types)}"
        )


def main() -> None:
    """Recreate and save the standardized annotation table."""
    metadata = load_source_metadata(SOURCE_METADATA_PATH)
    validate_original_metadata(metadata)

    annotations = standardize_annotations(metadata)
    validate_standardized_annotations(annotations)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    annotations.to_csv(
        OUTPUT_ANNOTATIONS_PATH,
        index=False,
        compression="gzip",
    )

    print(f"Wrote: {OUTPUT_ANNOTATIONS_PATH}")
    print(f"Annotated nuclei: {len(annotations):,}")
    print(f"Samples: {annotations['sample'].nunique()}")
    print(f"Clusters: {annotations['seurat_cluster'].nunique()}")

    print("\nCell-type counts:")
    print(
        annotations["cell_type"]
        .value_counts()
        .to_string()
    )

    print("\nOutput columns:")
    for column in annotations.columns:
        print(f"  {column}")


if __name__ == "__main__":
    main()