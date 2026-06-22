"""Split each snATAC peak-by-cell matrix by annotated cell type.

This step:

1. Loads the matched-column table for each sample.
2. Validates every matched barcode against the original barcode file.
3. Loads the original sparse peak-by-cell matrix.
4. Selects columns belonging to each annotated cell type.
5. Writes one sparse Matrix Market matrix per sample and cell type.

This step does not binarize values and does not filter peaks.
Therefore, every output matrix has exactly the same rows as its
corresponding original sample matrix.
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

from constants import (
    ATAC_SEQ_DIRS,
    BARCODES_SUFFIX,
    CELL_TYPE_MATRICES_DIR,
    CELL_TYPE_MATRIX_SPLIT_SUMMARY_PATH,
    CELL_TYPE_STANDARDIZATION,
    MATCHING_ANNOTATED_COLUMNS_DIR,
    MATRIX_SUFFIX,
    PEAKS_SUFFIX,
)


MATCHED_COLUMNS_FILENAME = "matched_columns.csv.gz"

ORIGINAL_COLUMN_INDEX = "matrix_column_index_zero_based"
OUTPUT_COLUMN_INDEX = "cell_type_column_index_zero_based"


def find_single_file(directory: Path, suffix: str) -> Path:
    """Return the unique file in a directory ending with the given suffix."""
    matches = sorted(directory.glob(f"*{suffix}"))

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one file ending with {suffix!r} in "
            f"{directory}, but found {len(matches)}."
        )

    return matches[0]


def load_barcodes(path: Path) -> pd.Series:
    """Load ordered barcodes from a compressed one-column TSV file."""
    barcodes = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["raw_barcode"],
        dtype=str,
        compression="gzip",
    )["raw_barcode"]

    if barcodes.isna().any():
        raise ValueError(f"Missing barcodes found in {path}.")

    if barcodes.duplicated().any():
        duplicated = barcodes[barcodes.duplicated()].head().tolist()
        raise ValueError(
            f"Duplicate barcodes found in {path}. "
            f"Examples: {duplicated}"
        )

    return barcodes


def load_matched_columns(path: Path) -> pd.DataFrame:
    """Load and validate the matched-column table for one sample."""
    matched = pd.read_csv(
        path,
        dtype={
            "sample_directory": str,
            "seurat_sample": str,
            "raw_barcode": str,
            "seurat_cell_id": str,
            "cell_type_original": str,
            "cell_type": str,
            "condition": str,
        },
        low_memory=False,
    )

    required_columns = {
        "sample_directory",
        "seurat_sample",
        ORIGINAL_COLUMN_INDEX,
        "raw_barcode",
        "cell_type",
        "seurat_cluster",
    }

    missing_columns = required_columns - set(matched.columns)

    if missing_columns:
        raise ValueError(
            f"{path} is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    matched[ORIGINAL_COLUMN_INDEX] = pd.to_numeric(
        matched[ORIGINAL_COLUMN_INDEX],
        errors="raise",
    ).astype(np.int64)

    matched["seurat_cluster"] = pd.to_numeric(
        matched["seurat_cluster"],
        errors="raise",
    ).astype(np.int64)

    if matched[ORIGINAL_COLUMN_INDEX].duplicated().any():
        raise ValueError(
            f"A matrix column occurs more than once in {path}."
        )

    if matched["raw_barcode"].duplicated().any():
        raise ValueError(
            f"A raw barcode occurs more than once in {path}."
        )

    valid_cell_types = set(CELL_TYPE_STANDARDIZATION.values())
    observed_cell_types = set(matched["cell_type"].dropna().unique())

    invalid_cell_types = observed_cell_types - valid_cell_types

    if invalid_cell_types:
        raise ValueError(
            f"Unrecognized standardized cell types in {path}: "
            f"{sorted(invalid_cell_types)}"
        )

    return matched.sort_values(
        ORIGINAL_COLUMN_INDEX
    ).reset_index(drop=True)


def validate_matched_columns(
    matched: pd.DataFrame,
    source_barcodes: pd.Series,
    sample_name: str,
) -> None:
    """Verify that matched indices point to the expected source barcodes."""
    if matched.empty:
        return

    indices = matched[ORIGINAL_COLUMN_INDEX].to_numpy(dtype=np.int64)

    if indices.min() < 0:
        raise ValueError(
            f"{sample_name}: a matched column index is negative."
        )

    if indices.max() >= len(source_barcodes):
        raise ValueError(
            f"{sample_name}: matched column index {indices.max()} "
            f"exceeds the largest valid index "
            f"{len(source_barcodes) - 1}."
        )

    barcodes_at_indices = (
        source_barcodes.iloc[indices]
        .reset_index(drop=True)
    )

    expected_barcodes = matched["raw_barcode"].reset_index(drop=True)

    mismatch_mask = barcodes_at_indices != expected_barcodes

    if mismatch_mask.any():
        first_mismatch = int(np.flatnonzero(mismatch_mask.to_numpy())[0])

        raise ValueError(
            f"{sample_name}: matched barcode validation failed at "
            f"matched-table row {first_mismatch}. "
            f"Source barcode is "
            f"{barcodes_at_indices.iloc[first_mismatch]!r}, but the "
            f"annotation table contains "
            f"{expected_barcodes.iloc[first_mismatch]!r}."
        )


def load_sparse_matrix(
    matrix_path: Path,
    expected_columns: int,
) -> sparse.csc_matrix:
    """Load a Matrix Market file and convert it to CSC for column slicing."""
    print(f"  Loading matrix: {matrix_path.name}")

    matrix = mmread(matrix_path)

    if not sparse.issparse(matrix):
        raise TypeError(
            f"Expected a sparse matrix in {matrix_path}, "
            f"but received {type(matrix).__name__}."
        )

    if matrix.shape[1] != expected_columns:
        raise ValueError(
            f"{matrix_path} contains {matrix.shape[1]:,} columns, "
            f"but its barcode file contains {expected_columns:,}."
        )

    # CSC is the appropriate sparse representation for selecting columns.
    matrix = matrix.tocsc()

    if np.issubdtype(matrix.dtype, np.floating):
        rounded = np.rint(matrix.data)

        if not np.allclose(matrix.data, rounded):
            raise ValueError(
                f"{matrix_path} contains non-integer accessibility values."
            )

        matrix.data = rounded.astype(np.int64)

    return matrix


def write_matrix_market_gzip(
    matrix: sparse.spmatrix,
    output_path: Path,
) -> None:
    """Write a sparse integer matrix as a gzip-compressed Matrix Market file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(
        output_path,
        mode="wb",
        compresslevel=6,
    ) as handle:
        mmwrite(
            handle,
            matrix,
            field="integer",
            symmetry="general",
        )


def write_barcodes(
    barcodes: pd.Series,
    output_path: Path,
) -> None:
    """Write ordered barcodes as a compressed one-column TSV file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    barcodes.to_csv(
        output_path,
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )


def create_output_prefix(
    sample_name: str,
    cell_type: str,
) -> str:
    """Create a filesystem-safe prefix for one sample–cell-type output."""
    safe_cell_type = (
        cell_type.strip()
        .replace(" ", "_")
        .replace("/", "_")
    )

    return f"{sample_name}_{safe_cell_type}"


def split_one_sample(
    sample_directory: Path,
) -> list[dict[str, Any]]:
    """Split one original sample matrix into cell-type-specific matrices."""
    sample_name = sample_directory.name

    matched_path = (
        MATCHING_ANNOTATED_COLUMNS_DIR
        / sample_name
        / MATCHED_COLUMNS_FILENAME
    )

    if not matched_path.is_file():
        raise FileNotFoundError(
            f"Matched-column table was not found for {sample_name}: "
            f"{matched_path}"
        )

    barcode_path = find_single_file(
        sample_directory,
        BARCODES_SUFFIX,
    )
    matrix_path = find_single_file(
        sample_directory,
        MATRIX_SUFFIX,
    )
    peaks_path = find_single_file(
        sample_directory,
        PEAKS_SUFFIX,
    )

    source_barcodes = load_barcodes(barcode_path)
    matched = load_matched_columns(matched_path)

    observed_sample_names = matched[
        "sample_directory"
    ].dropna().unique()

    if len(observed_sample_names) > 1:
        raise ValueError(
            f"{matched_path} contains multiple sample directories: "
            f"{observed_sample_names.tolist()}"
        )

    if (
        len(observed_sample_names) == 1
        and observed_sample_names[0] != sample_name
    ):
        raise ValueError(
            f"{matched_path} belongs to "
            f"{observed_sample_names[0]!r}, not {sample_name!r}."
        )

    validate_matched_columns(
        matched=matched,
        source_barcodes=source_barcodes,
        sample_name=sample_name,
    )

    if matched.empty:
        print(f"{sample_name}: no matched columns; skipping matrix load.")

        return [
            {
                "sample_directory": sample_name,
                "seurat_sample": None,
                "cell_type": None,
                "status": "skipped_no_matched_columns",
                "original_matrix_rows": None,
                "original_matrix_columns": len(source_barcodes),
                "total_matched_columns": 0,
                "output_matrix_rows": None,
                "output_matrix_columns": 0,
                "output_matrix_nonzero_entries": 0,
                "matrix_output_path": None,
                "barcodes_output_path": None,
                "peaks_output_path": None,
                "metadata_output_path": None,
            }
        ]

    matrix = load_sparse_matrix(
        matrix_path=matrix_path,
        expected_columns=len(source_barcodes),
    )

    original_n_rows, original_n_columns = matrix.shape
    sample_output_directory = CELL_TYPE_MATRICES_DIR / sample_name
    sample_output_directory.mkdir(parents=True, exist_ok=True)

    summary_records: list[dict[str, Any]] = []
    written_columns = 0

    grouped = matched.groupby(
        "cell_type",
        sort=True,
        observed=True,
    )

    for cell_type, cell_metadata in grouped:
        cell_metadata = (
            cell_metadata.sort_values(ORIGINAL_COLUMN_INDEX)
            .reset_index(drop=True)
            .copy()
        )

        source_column_indices = cell_metadata[
            ORIGINAL_COLUMN_INDEX
        ].to_numpy(dtype=np.int64)

        cell_type_matrix = matrix[:, source_column_indices]

        if cell_type_matrix.shape[0] != original_n_rows:
            raise RuntimeError(
                f"{sample_name}/{cell_type}: output row count changed."
            )

        if cell_type_matrix.shape[1] != len(cell_metadata):
            raise RuntimeError(
                f"{sample_name}/{cell_type}: output column count does "
                f"not equal the annotation count."
            )

        cell_metadata.insert(
            0,
            OUTPUT_COLUMN_INDEX,
            np.arange(len(cell_metadata), dtype=np.int64),
        )

        output_prefix = create_output_prefix(
            sample_name=sample_name,
            cell_type=cell_type,
        )

        cell_type_output_directory = (
            sample_output_directory / cell_type
        )
        cell_type_output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        matrix_output_path = (
            cell_type_output_directory
            / f"{output_prefix}{MATRIX_SUFFIX}"
        )
        barcodes_output_path = (
            cell_type_output_directory
            / f"{output_prefix}{BARCODES_SUFFIX}"
        )
        peaks_output_path = (
            cell_type_output_directory
            / f"{output_prefix}{PEAKS_SUFFIX}"
        )
        metadata_output_path = (
            cell_type_output_directory
            / f"{output_prefix}_cell_metadata.csv.gz"
        )

        print(
            f"  {cell_type}: "
            f"{cell_type_matrix.shape[0]:,} peaks × "
            f"{cell_type_matrix.shape[1]:,} cells; "
            f"{cell_type_matrix.nnz:,} nonzero entries"
        )

        write_matrix_market_gzip(
            matrix=cell_type_matrix,
            output_path=matrix_output_path,
        )

        write_barcodes(
            barcodes=cell_metadata["raw_barcode"],
            output_path=barcodes_output_path,
        )

        # Peak rows are unchanged at this stage.
        shutil.copy2(
            peaks_path,
            peaks_output_path,
        )

        cell_metadata.to_csv(
            metadata_output_path,
            index=False,
            compression="gzip",
        )

        written_columns += cell_type_matrix.shape[1]

        seurat_samples = cell_metadata[
            "seurat_sample"
        ].dropna().unique()

        seurat_sample = (
            seurat_samples[0]
            if len(seurat_samples) == 1
            else None
        )

        summary_records.append(
            {
                "sample_directory": sample_name,
                "seurat_sample": seurat_sample,
                "cell_type": cell_type,
                "status": "written",
                "original_matrix_rows": original_n_rows,
                "original_matrix_columns": original_n_columns,
                "total_matched_columns": len(matched),
                "output_matrix_rows": cell_type_matrix.shape[0],
                "output_matrix_columns": cell_type_matrix.shape[1],
                "output_matrix_nonzero_entries": cell_type_matrix.nnz,
                "matrix_output_path": str(matrix_output_path.name),
                "barcodes_output_path": str(barcodes_output_path.name),
                "peaks_output_path": str(peaks_output_path.name),
                "metadata_output_path": str(metadata_output_path.name),
            }
        )

        del cell_type_matrix

    if written_columns != len(matched):
        raise RuntimeError(
            f"{sample_name}: wrote {written_columns:,} columns across "
            f"cell types, but expected {len(matched):,} matched columns."
        )

    return summary_records


def validate_split_summary(summary: pd.DataFrame) -> None:
    """Validate global invariants of the completed split."""
    written = summary.loc[summary["status"] == "written"].copy()

    if written.empty:
        raise RuntimeError("No cell-type matrices were written.")

    if not (
        written["original_matrix_rows"]
        == written["output_matrix_rows"]
    ).all():
        raise RuntimeError(
            "At least one split matrix does not preserve the original "
            "number of peak rows."
        )

    output_columns_by_sample = (
        written.groupby("sample_directory", observed=True)[
            "output_matrix_columns"
        ]
        .sum()
    )

    expected_columns_by_sample = (
        written.groupby("sample_directory", observed=True)[
            "total_matched_columns"
        ]
        .first()
    )

    if not output_columns_by_sample.equals(expected_columns_by_sample):
        comparison = pd.concat(
            {
                "written": output_columns_by_sample,
                "expected": expected_columns_by_sample,
            },
            axis=1,
        )

        raise RuntimeError(
            "Split column totals do not equal matched-column totals:\n"
            f"{comparison}"
        )


def main() -> None:
    """Split all available snATAC matrices by annotated cell type."""
    CELL_TYPE_MATRICES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    all_summary_records: list[dict[str, Any]] = []

    for sample_directory in ATAC_SEQ_DIRS:
        print(f"\nProcessing {sample_directory.name}")

        sample_records = split_one_sample(sample_directory)
        all_summary_records.extend(sample_records)

    summary = pd.DataFrame.from_records(all_summary_records)

    validate_split_summary(summary)

    summary.to_csv(
        CELL_TYPE_MATRIX_SPLIT_SUMMARY_PATH,
        index=False,
    )

    written = summary.loc[summary["status"] == "written"]

    print("\nCompleted matrix splitting.")
    print(
        f"Cell-type matrices written: {len(written):,}"
    )
    print(
        "Total output columns: "
        f"{int(written['output_matrix_columns'].sum()):,}"
    )
    print(
        "Expected matched columns: 55,014"
    )
    print(
        f"Summary written to: "
        f"{CELL_TYPE_MATRIX_SPLIT_SUMMARY_PATH}"
    )

    print("\nOutput cell counts:")
    print(
        written.groupby("cell_type", observed=True)[
            "output_matrix_columns"
        ]
        .sum()
        .sort_values(ascending=False)
        .to_string()
    )


if __name__ == "__main__":
    main()