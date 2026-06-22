"""Report the number of unique values in every cell-type matrix.

The script analyzes both:

1. Original cell-type-specific count matrices.
2. Binarized cell-type-specific matrices.

Because the matrices are sparse, zero values are generally implicit and do not
appear in ``matrix.data``. The report therefore includes:

- number of unique stored nonzero values;
- whether zero is present implicitly;
- total number of unique matrix values, including zero;
- the actual unique nonzero values.

The output is written to:

    results/preprocessing/matrix_splitting/binarization/matrix_unique_values.csv
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from constants import (
    BINARIZED_CELL_TYPE_MATRICES_DIR,
    CELL_TYPE_MATRICES_DIR,
    MATRIX_SUFFIX,
    PREPROCESSING_RES_DIR,
)


OUTPUT_PATH = (
    PREPROCESSING_RES_DIR
    / "matrix_splitting"
    / "binarization"
    / "matrix_unique_values.csv"
)


def find_matrices(root_directory: Path) -> list[Path]:
    """Find all Matrix Market files beneath a root directory."""
    if not root_directory.is_dir():
        raise FileNotFoundError(
            f"Matrix root directory does not exist: {root_directory}"
        )

    return sorted(
        root_directory.rglob(f"*{MATRIX_SUFFIX}")
    )


def read_sparse_matrix(path: Path) -> sparse.coo_matrix:
    """Read a gzip-compressed Matrix Market file as a sparse COO matrix."""
    with gzip.open(path, mode="rb") as handle:
        matrix = mmread(handle)

    if not sparse.issparse(matrix):
        raise TypeError(
            f"Expected a sparse matrix in {path}, "
            f"but found {type(matrix).__name__}."
        )

    matrix = matrix.tocoo()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()

    if matrix.nnz > 0 and not np.isfinite(matrix.data).all():
        raise ValueError(
            f"Matrix contains non-finite values: {path}"
        )

    return matrix


def extract_matrix_identity(
    matrix_path: Path,
    root_directory: Path,
) -> tuple[str, str]:
    """Extract sample and cell type from a matrix path.

    Expected structure:

        <root>/<sample>/<cell_type>/<matrix_file>
    """
    relative_path = matrix_path.relative_to(
        root_directory
    )

    if len(relative_path.parts) < 3:
        raise ValueError(
            f"Unexpected matrix path structure: {matrix_path}"
        )

    sample_directory = relative_path.parts[0]
    cell_type = relative_path.parts[1]

    return sample_directory, cell_type


def calculate_unique_value_statistics(
    matrix: sparse.coo_matrix,
) -> dict[str, Any]:
    """Calculate unique-value statistics for one sparse matrix."""
    n_rows, n_columns = matrix.shape
    n_total_entries = n_rows * n_columns

    if matrix.nnz == 0:
        unique_nonzero_values = np.array([], dtype=np.int64)
    else:
        rounded_values = np.rint(matrix.data)

        if not np.allclose(matrix.data, rounded_values):
            raise ValueError(
                "Matrix contains non-integer values."
            )

        integer_values = rounded_values.astype(
            np.int64,
            copy=False,
        )

        if (integer_values <= 0).any():
            raise ValueError(
                "Matrix contains non-positive stored values."
            )

        unique_nonzero_values = np.unique(
            integer_values
        )

    has_implicit_zero = matrix.nnz < n_total_entries

    n_unique_nonzero_values = len(
        unique_nonzero_values
    )

    n_unique_values_including_zero = (
        n_unique_nonzero_values
        + int(has_implicit_zero)
    )

    unique_values_text = ";".join(
        str(int(value))
        for value in unique_nonzero_values
    )

    return {
        "matrix_rows": n_rows,
        "matrix_columns": n_columns,
        "total_entries": n_total_entries,
        "nonzero_entries": matrix.nnz,
        "has_zero": has_implicit_zero,
        "n_unique_nonzero_values": (
            n_unique_nonzero_values
        ),
        "n_unique_values_including_zero": (
            n_unique_values_including_zero
        ),
        "minimum_nonzero_value": (
            int(unique_nonzero_values.min())
            if n_unique_nonzero_values > 0
            else None
        ),
        "maximum_nonzero_value": (
            int(unique_nonzero_values.max())
            if n_unique_nonzero_values > 0
            else None
        ),
        "unique_nonzero_values": unique_values_text,
    }


def analyze_matrix(
    matrix_path: Path,
    root_directory: Path,
    matrix_type: str,
) -> dict[str, Any]:
    """Analyze one matrix and return a report row."""
    sample_directory, cell_type = (
        extract_matrix_identity(
            matrix_path=matrix_path,
            root_directory=root_directory,
        )
    )

    matrix = read_sparse_matrix(
        matrix_path
    )

    statistics = calculate_unique_value_statistics(
        matrix
    )

    return {
        "matrix_type": matrix_type,
        "sample_directory": sample_directory,
        "cell_type": cell_type,
        "matrix_filename": matrix_path.name,
        **statistics,
        "matrix_path": str(matrix_path),
    }


def analyze_matrix_collection(
    root_directory: Path,
    matrix_type: str,
) -> list[dict[str, Any]]:
    """Analyze every matrix in one matrix collection."""
    matrices = find_matrices(
        root_directory
    )

    print(
        f"Found {len(matrices):,} {matrix_type} matrices."
    )

    records: list[dict[str, Any]] = []

    for matrix_path in matrices:
        sample_directory, cell_type = (
            extract_matrix_identity(
                matrix_path=matrix_path,
                root_directory=root_directory,
            )
        )

        print(
            f"Processing {matrix_type}: "
            f"{sample_directory} / {cell_type}"
        )

        records.append(
            analyze_matrix(
                matrix_path=matrix_path,
                root_directory=root_directory,
                matrix_type=matrix_type,
            )
        )

    return records


def validate_binarized_matrices(
    report: pd.DataFrame,
) -> None:
    """Verify that all binarized matrices contain only zero and one."""
    binarized = report.loc[
        report["matrix_type"] == "binarized"
    ]

    invalid = binarized.loc[
        (binarized["n_unique_nonzero_values"] != 1)
        | (binarized["minimum_nonzero_value"] != 1)
        | (binarized["maximum_nonzero_value"] != 1)
    ]

    if not invalid.empty:
        raise RuntimeError(
            "At least one binarized matrix contains stored values "
            "other than 1:\n"
            + invalid[
                [
                    "sample_directory",
                    "cell_type",
                    "n_unique_nonzero_values",
                    "minimum_nonzero_value",
                    "maximum_nonzero_value",
                ]
            ].to_string(index=False)
        )


def main() -> None:
    """Analyze original and binarized matrices and write the report."""
    records: list[dict[str, Any]] = []

    records.extend(
        analyze_matrix_collection(
            root_directory=CELL_TYPE_MATRICES_DIR,
            matrix_type="original",
        )
    )

    records.extend(
        analyze_matrix_collection(
            root_directory=(
                BINARIZED_CELL_TYPE_MATRICES_DIR
            ),
            matrix_type="binarized",
        )
    )

    report = pd.DataFrame.from_records(
        records
    )

    report = report.sort_values(
        [
            "sample_directory",
            "cell_type",
            "matrix_type",
        ]
    ).reset_index(drop=True)

    validate_binarized_matrices(
        report
    )

    OUTPUT_PATH.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    report.to_csv(
        OUTPUT_PATH,
        index=False,
    )

    print("\nUnique-value summary:")
    print(
        report[
            [
                "matrix_type",
                "sample_directory",
                "cell_type",
                "n_unique_nonzero_values",
                "n_unique_values_including_zero",
                "minimum_nonzero_value",
                "maximum_nonzero_value",
            ]
        ].to_string(index=False)
    )

    print(
        f"\nAnalyzed {len(report):,} matrices."
    )
    print(
        f"Report written to: {OUTPUT_PATH}"
    )


if __name__ == "__main__":
    main()