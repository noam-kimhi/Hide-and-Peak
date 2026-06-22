"""Binarize cell-type-specific snATAC-seq matrices.

Every positive accessibility count is replaced by 1:

    0 -> 0
    x -> 1, for every x > 0

The matrix dimensions, row order, column order, and sparsity pattern are
preserved. Original count matrices are not modified.
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
    BARCODES_SUFFIX,
    BINARIZED_CELL_TYPE_MATRICES_DIR,
    CELL_TYPE_MATRICES_DIR,
    MATRIX_BINARIZATION_SUMMARY_PATH,
    MATRIX_SUFFIX,
    PEAKS_SUFFIX,
)


CELL_METADATA_SUFFIX = "_cell_metadata.csv.gz"


def read_sparse_matrix(path: Path) -> sparse.csr_matrix:
    """Read a gzip-compressed Matrix Market file as a CSR matrix."""
    with gzip.open(path, mode="rb") as handle:
        matrix = mmread(handle)

    if not sparse.issparse(matrix):
        raise TypeError(
            f"Expected a sparse matrix in {path}, "
            f"but found {type(matrix).__name__}."
        )

    matrix = matrix.tocsr()

    raw_nnz = matrix.nnz

    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()

    if matrix.nnz != raw_nnz:
        raise ValueError(
            f"{path} contained duplicate coordinates or explicit zero "
            f"entries. Stored entries changed from {raw_nnz:,} to "
            f"{matrix.nnz:,} during sparse-matrix cleanup."
        )

    if matrix.nnz > 0:
        if not np.isfinite(matrix.data).all():
            raise ValueError(
                f"{path} contains non-finite matrix values."
            )

        if (matrix.data < 0).any():
            raise ValueError(
                f"{path} contains negative accessibility values."
            )

        if (matrix.data == 0).any():
            raise ValueError(
                f"{path} contains explicitly stored zero values."
            )

    return matrix


def binarize_sparse_matrix(
    matrix: sparse.csr_matrix,
) -> sparse.csr_matrix:
    """Return a binary matrix with the same sparsity pattern."""
    binary = matrix.copy()

    binary.data = np.ones(
        binary.nnz,
        dtype=np.int8,
    )

    return binary


def write_sparse_matrix(
    matrix: sparse.spmatrix,
    path: Path,
) -> None:
    """Write a sparse integer matrix as gzip-compressed Matrix Market."""
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with gzip.open(
        path,
        mode="wb",
        compresslevel=6,
    ) as handle:
        mmwrite(
            handle,
            matrix,
            field="integer",
            symmetry="general",
        )


def get_output_directory(
    source_matrix_path: Path,
) -> Path:
    """Return the mirrored output directory for a source matrix."""
    relative_parent = source_matrix_path.parent.relative_to(
        CELL_TYPE_MATRICES_DIR
    )

    return (
        BINARIZED_CELL_TYPE_MATRICES_DIR
        / relative_parent
    )


def get_file_prefix(matrix_path: Path) -> str:
    """Remove the matrix suffix from a matrix filename."""
    if not matrix_path.name.endswith(MATRIX_SUFFIX):
        raise ValueError(
            f"Unexpected matrix filename: {matrix_path.name}"
        )

    return matrix_path.name[
        :-len(MATRIX_SUFFIX)
    ]


def copy_associated_files(
    source_directory: Path,
    output_directory: Path,
    prefix: str,
) -> dict[str, Path]:
    """Copy barcode, peak, and cell-metadata files beside the binary matrix."""
    associated_suffixes = {
        "barcodes_output_path": BARCODES_SUFFIX,
        "peaks_output_path": PEAKS_SUFFIX,
        "metadata_output_path": CELL_METADATA_SUFFIX,
    }

    copied_paths: dict[str, Path] = {}

    for summary_field, suffix in associated_suffixes.items():
        source_path = source_directory / f"{prefix}{suffix}"

        if not source_path.is_file():
            raise FileNotFoundError(
                f"Required associated file was not found: {source_path}"
            )

        output_path = output_directory / source_path.name

        shutil.copy2(
            source_path,
            output_path,
        )

        copied_paths[summary_field] = output_path

    return copied_paths


def validate_binary_matrix(
    original: sparse.csr_matrix,
    binary: sparse.csr_matrix,
    source_path: Path,
) -> None:
    """Validate that binarization changed only stored values."""
    if original.shape != binary.shape:
        raise RuntimeError(
            f"{source_path}: matrix dimensions changed from "
            f"{original.shape} to {binary.shape}."
        )

    if original.nnz != binary.nnz:
        raise RuntimeError(
            f"{source_path}: nonzero count changed from "
            f"{original.nnz:,} to {binary.nnz:,}."
        )

    if binary.nnz > 0 and not np.all(binary.data == 1):
        unique_values = np.unique(binary.data)

        raise RuntimeError(
            f"{source_path}: binarized matrix contains values other "
            f"than 1: {unique_values.tolist()}"
        )

    if not np.array_equal(
        original.indptr,
        binary.indptr,
    ):
        raise RuntimeError(
            f"{source_path}: sparse row structure changed."
        )

    if not np.array_equal(
        original.indices,
        binary.indices,
    ):
        raise RuntimeError(
            f"{source_path}: sparse column indices changed."
        )


def binarize_one_matrix(
    source_matrix_path: Path,
) -> dict[str, Any]:
    """Binarize one cell-type-specific matrix and write its outputs."""
    relative_path = source_matrix_path.relative_to(
        CELL_TYPE_MATRICES_DIR
    )

    if len(relative_path.parts) < 3:
        raise ValueError(
            f"Unexpected matrix location: {source_matrix_path}"
        )

    sample_directory = relative_path.parts[0]
    cell_type = relative_path.parts[1]

    prefix = get_file_prefix(source_matrix_path)

    output_directory = get_output_directory(
        source_matrix_path
    )
    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_matrix_path = (
        output_directory
        / source_matrix_path.name
    )

    print(
        f"Processing {sample_directory} / {cell_type}"
    )

    matrix = read_sparse_matrix(
        source_matrix_path
    )

    input_nnz = matrix.nnz

    if input_nnz > 0:
        input_min_value: int | float | None = (
            matrix.data.min().item()
        )
        input_max_value: int | float | None = (
            matrix.data.max().item()
        )

        entries_equal_one = int(
            np.count_nonzero(matrix.data == 1)
        )
        entries_greater_than_one = int(
            np.count_nonzero(matrix.data > 1)
        )
    else:
        input_min_value = None
        input_max_value = None
        entries_equal_one = 0
        entries_greater_than_one = 0

    binary = binarize_sparse_matrix(
        matrix
    )

    validate_binary_matrix(
        original=matrix,
        binary=binary,
        source_path=source_matrix_path,
    )

    write_sparse_matrix(
        matrix=binary,
        path=output_matrix_path,
    )

    copied_paths = copy_associated_files(
        source_directory=source_matrix_path.parent,
        output_directory=output_directory,
        prefix=prefix,
    )

    summary: dict[str, Any] = {
        "sample_directory": sample_directory,
        "cell_type": cell_type,
        "status": "written",
        "input_matrix_rows": matrix.shape[0],
        "input_matrix_columns": matrix.shape[1],
        "input_nonzero_entries": input_nnz,
        "input_min_nonzero_value": input_min_value,
        "input_max_value": input_max_value,
        "input_entries_equal_one": entries_equal_one,
        "input_entries_greater_than_one": (
            entries_greater_than_one
        ),
        "output_matrix_rows": binary.shape[0],
        "output_matrix_columns": binary.shape[1],
        "output_nonzero_entries": binary.nnz,
        "dimensions_preserved": (
            matrix.shape == binary.shape
        ),
        "nonzero_entries_preserved": (
            matrix.nnz == binary.nnz
        ),
        "output_nonzero_value": (
            1 if binary.nnz > 0 else None
        ),
        "input_matrix_path": str(
            source_matrix_path.name
        ),
        "output_matrix_path": str(
            output_matrix_path.name
        ),
        **{
            key: str(value)
            for key, value in copied_paths.items()
        },
    }

    print(
        f"  Shape: {binary.shape[0]:,} × "
        f"{binary.shape[1]:,}"
    )
    print(
        f"  Nonzero entries: {binary.nnz:,}"
    )
    print(
        f"  Entries changed to 1: "
        f"{entries_greater_than_one:,}"
    )

    return summary


def find_input_matrices() -> list[Path]:
    """Find every cell-type-specific count matrix."""
    matrices = sorted(
        CELL_TYPE_MATRICES_DIR.rglob(
            f"*{MATRIX_SUFFIX}"
        )
    )

    if not matrices:
        raise FileNotFoundError(
            "No split matrices were found beneath "
            f"{CELL_TYPE_MATRICES_DIR}."
        )

    return matrices


def validate_global_summary(
    summary: pd.DataFrame,
    expected_matrix_count: int,
) -> None:
    """Validate all completed binarization outputs."""
    if len(summary) != expected_matrix_count:
        raise RuntimeError(
            f"Expected {expected_matrix_count:,} output matrices, "
            f"but summary contains {len(summary):,} rows."
        )

    if not summary["dimensions_preserved"].all():
        raise RuntimeError(
            "At least one matrix changed dimensions."
        )

    if not summary["nonzero_entries_preserved"].all():
        raise RuntimeError(
            "At least one matrix changed its number of nonzero entries."
        )

    if not (
        summary["input_matrix_rows"]
        == summary["output_matrix_rows"]
    ).all():
        raise RuntimeError(
            "At least one output has a different row count."
        )

    if not (
        summary["input_matrix_columns"]
        == summary["output_matrix_columns"]
    ).all():
        raise RuntimeError(
            "At least one output has a different column count."
        )

    if not (
        summary["input_nonzero_entries"]
        == summary["output_nonzero_entries"]
    ).all():
        raise RuntimeError(
            "At least one output has a different nonzero count."
        )


def main() -> None:
    """Binarize all split cell-type-specific matrices."""
    input_matrices = find_input_matrices()

    print(
        f"Found {len(input_matrices):,} matrices to binarize."
    )

    summary_records: list[dict[str, Any]] = []

    for matrix_path in input_matrices:
        summary_records.append(
            binarize_one_matrix(matrix_path)
        )

    summary = pd.DataFrame.from_records(
        summary_records
    )

    validate_global_summary(
        summary=summary,
        expected_matrix_count=len(input_matrices),
    )

    summary.to_csv(
        MATRIX_BINARIZATION_SUMMARY_PATH,
        index=False,
    )

    total_input_nnz = int(
        summary["input_nonzero_entries"].sum()
    )
    total_output_nnz = int(
        summary["output_nonzero_entries"].sum()
    )
    total_changed = int(
        summary["input_entries_greater_than_one"].sum()
    )

    print("\nBinarization completed successfully.")
    print(
        f"Matrices written: {len(summary):,}"
    )
    print(
        f"Total input nonzero entries: "
        f"{total_input_nnz:,}"
    )
    print(
        f"Total output nonzero entries: "
        f"{total_output_nnz:,}"
    )
    print(
        f"Entries changed from >1 to 1: "
        f"{total_changed:,}"
    )
    print(
        f"Summary written to: "
        f"{MATRIX_BINARIZATION_SUMMARY_PATH}"
    )


if __name__ == "__main__":
    main()