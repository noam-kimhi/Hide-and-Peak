"""Plot original and binarized matrix-value distributions by cell type.

For every cell-type-specific matrix inside each replicate, this script creates
one PNG containing two side-by-side histograms:

- Left: original non-binarized matrix values.
- Right: binarized matrix values.

The output structure is:

    results/preprocessing/matrix_splitting/binarization/
    ├── GSM8619363_MASH_rep1/
    │   ├── Cholangiocyte.png
    │   ├── Endothelial.png
    │   ├── Hepatocytes.png
    │   └── ...
    ├── GSM8619364_MASH_rep2/
    │   └── ...
    └── ...

Implicit zero entries are excluded because the matrices are sparse. Including
them would dominate both distributions and obscure the binarization effect.
"""

from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.io import mmread

from constants import (
    ATAC_SEQ_DIRS,
    BINARIZED_CELL_TYPE_MATRICES_DIR,
    CELL_TYPE_MATRICES_DIR,
    MATRIX_BINARIZATION_PLOTS_DIR,
    MATRIX_SUFFIX,
    PLOT_DPI,
    PLOTS_MAIN_COLOR,
)


def find_sample_matrices(
    root_directory: Path,
    sample_name: str,
) -> dict[Path, Path]:
    """Find all cell-type matrices belonging to one replicate.

    Args:
        root_directory: Root directory containing all replicate directories.
        sample_name: Name of the replicate directory.

    Returns:
        Mapping from each matrix path relative to the replicate directory to
        its absolute path.
    """
    sample_directory = root_directory / sample_name

    if not sample_directory.is_dir():
        return {}

    return {
        matrix_path.relative_to(sample_directory): matrix_path
        for matrix_path in sorted(
            sample_directory.rglob(f"*{MATRIX_SUFFIX}")
        )
    }


def read_sparse_matrix(path: Path) -> sparse.coo_matrix:
    """Read and validate a gzip-compressed Matrix Market sparse matrix.

    Args:
        path: Path to the compressed Matrix Market file.

    Returns:
        Validated sparse matrix in COO format.

    Raises:
        TypeError: If the file does not contain a sparse matrix.
        ValueError: If the matrix is empty or contains invalid values.
    """
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

    if matrix.nnz == 0:
        raise ValueError(
            f"Matrix contains no nonzero entries: {path}"
        )

    if not np.isfinite(matrix.data).all():
        raise ValueError(
            f"Matrix contains non-finite values: {path}"
        )

    return matrix


def integer_value_counts(
    matrix: sparse.spmatrix,
    path: Path,
) -> Counter[int]:
    """Count the stored integer values in a sparse matrix.

    Args:
        matrix: Sparse matrix whose stored values should be counted.
        path: Matrix path, used in validation error messages.

    Returns:
        Counter mapping each positive integer value to its frequency.

    Raises:
        ValueError: If values are non-integer or non-positive.
    """
    rounded_values = np.rint(matrix.data)

    if not np.allclose(matrix.data, rounded_values):
        raise ValueError(
            f"Matrix contains non-integer values: {path}"
        )

    integer_values = rounded_values.astype(
        np.int64,
        copy=False,
    )

    if (integer_values <= 0).any():
        raise ValueError(
            f"Matrix contains non-positive stored values: {path}"
        )

    values, counts = np.unique(
        integer_values,
        return_counts=True,
    )

    return Counter(
        {
            int(value): int(count)
            for value, count in zip(values, counts)
        }
    )


def validate_matrix_pair(
    original_matrix: sparse.spmatrix,
    binary_matrix: sparse.spmatrix,
    original_path: Path,
    binary_path: Path,
) -> None:
    """Validate an original and binarized matrix pair.

    Args:
        original_matrix: Original count matrix.
        binary_matrix: Corresponding binarized matrix.
        original_path: Path to the original matrix.
        binary_path: Path to the binarized matrix.

    Raises:
        RuntimeError: If matrix dimensions, nonzero counts, or binary values
            are inconsistent.
    """
    if original_matrix.shape != binary_matrix.shape:
        raise RuntimeError(
            f"Matrix dimensions differ:\n"
            f"  Original: {original_path} -> {original_matrix.shape}\n"
            f"  Binary:   {binary_path} -> {binary_matrix.shape}"
        )

    if original_matrix.nnz != binary_matrix.nnz:
        raise RuntimeError(
            f"Nonzero-entry counts differ:\n"
            f"  Original: {original_path} -> "
            f"{original_matrix.nnz:,}\n"
            f"  Binary:   {binary_path} -> "
            f"{binary_matrix.nnz:,}"
        )

    binary_values = np.unique(binary_matrix.data)

    if not np.array_equal(
        binary_values,
        np.array([1], dtype=binary_values.dtype),
    ):
        raise RuntimeError(
            f"Binarized matrix contains values other than 1:\n"
            f"  Matrix: {binary_path}\n"
            f"  Values: {binary_values.tolist()}"
        )


def load_matrix_pair_distributions(
    original_path: Path,
    binary_path: Path,
) -> tuple[Counter[int], Counter[int], tuple[int, int]]:
    """Load one original/binary matrix pair and calculate distributions.

    Args:
        original_path: Path to the original cell-type matrix.
        binary_path: Path to the corresponding binarized matrix.

    Returns:
        Tuple containing:
        - original nonzero-value counts,
        - binarized nonzero-value counts,
        - matrix shape.
    """
    original_matrix = read_sparse_matrix(original_path)
    binary_matrix = read_sparse_matrix(binary_path)

    validate_matrix_pair(
        original_matrix=original_matrix,
        binary_matrix=binary_matrix,
        original_path=original_path,
        binary_path=binary_path,
    )

    original_counts = integer_value_counts(
        matrix=original_matrix,
        path=original_path,
    )
    binary_counts = integer_value_counts(
        matrix=binary_matrix,
        path=binary_path,
    )

    if sum(original_counts.values()) != sum(binary_counts.values()):
        raise RuntimeError(
            "Original and binarized distributions contain different "
            "numbers of nonzero entries."
        )

    return (
        original_counts,
        binary_counts,
        original_matrix.shape,
    )


def plot_weighted_histogram(
    axis: plt.Axes,
    value_counts: Counter[int],
    title: str,
    color: str,
) -> None:
    """Plot an integer-valued histogram from aggregated value counts.

    Args:
        axis: Matplotlib axis on which to draw.
        value_counts: Mapping from matrix values to their frequencies.
        title: Plot title.
        color: Histogram color.
    """
    if not value_counts:
        raise ValueError(
            f"Cannot plot an empty distribution for {title!r}."
        )

    values = np.array(
        sorted(value_counts),
        dtype=np.int64,
    )
    counts = np.array(
        [
            value_counts[int(value)]
            for value in values
        ],
        dtype=np.int64,
    )

    bin_edges = np.arange(
        -0.5,
        values.max() + 1.5,
        1.0,
    )

    axis.hist(
        values,
        bins=bin_edges,
        weights=counts,
        color=color,
    )

    axis.set_title(title)
    axis.set_xlabel("Stored matrix value")
    axis.set_ylabel("Number of nonzero entries")
    axis.set_yscale("log")
    axis.grid(
        axis="y",
        alpha=0.25,
    )

    axis.set_xlim(
        -0.5,
        values.max() + 0.5,
    )

    axis.set_xticks(
        np.unique(
            np.concatenate(
                (
                    axis.get_xticks(),
                    np.array([0.0, 1.0]),
                )
            )
        )
    )


def create_cell_type_plot(
    sample_name: str,
    cell_type: str,
    original_counts: Counter[int],
    binary_counts: Counter[int],
    matrix_shape: tuple[int, int],
) -> Path:
    """Create a side-by-side comparison plot for one cell type.

    Args:
        sample_name: Replicate name.
        cell_type: Annotated cell type.
        original_counts: Original matrix-value distribution.
        binary_counts: Binarized matrix-value distribution.
        matrix_shape: Peak-by-cell matrix dimensions.

    Returns:
        Path of the generated PNG file.
    """
    sample_output_directory = (
        MATRIX_BINARIZATION_PLOTS_DIR
        / sample_name
    )
    sample_output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(14, 5),
    )

    plot_weighted_histogram(
        axis=axes[0],
        value_counts=original_counts,
        title="Original counts",
        color=PLOTS_MAIN_COLOR,
    )

    plot_weighted_histogram(
        axis=axes[1],
        value_counts=binary_counts,
        title="Binarized counts",
        color=PLOTS_MAIN_COLOR,
    )

    total_nonzero_entries = sum(
        original_counts.values()
    )

    entries_greater_than_one = sum(
        count
        for value, count in original_counts.items()
        if value > 1
    )

    fraction_changed = (
        entries_greater_than_one
        / total_nonzero_entries
        if total_nonzero_entries > 0
        else 0.0
    )

    n_peaks, n_cells = matrix_shape

    figure.suptitle(
        f"{sample_name} — {cell_type}"
    )

    figure.text(
        0.5,
        0.01,
        (
            f"{n_peaks:,} peaks × {n_cells:,} cells; "
            f"{total_nonzero_entries:,} nonzero entries; "
            f"{fraction_changed:.1%} changed from >1 to 1. "
            "Implicit zero entries are excluded."
        ),
        ha="center",
    )

    figure.tight_layout(
        rect=(0.0, 0.06, 1.0, 0.93)
    )

    output_path = (
        sample_output_directory
        / f"{cell_type}.png"
    )

    figure.savefig(
        output_path,
        dpi=PLOT_DPI,
        bbox_inches="tight",
    )

    plt.close(figure)

    return output_path


def get_cell_type(
    relative_matrix_path: Path,
) -> str:
    """Extract the cell type from a matrix path relative to its sample.

    The expected path structure is:

        <cell_type>/<matrix_filename>

    Args:
        relative_matrix_path: Matrix path relative to its sample directory.

    Returns:
        Cell-type directory name.

    Raises:
        ValueError: If the path does not have the expected structure.
    """
    if len(relative_matrix_path.parts) < 2:
        raise ValueError(
            "Expected a matrix path beneath a cell-type directory, "
            f"but found: {relative_matrix_path}"
        )

    return relative_matrix_path.parts[0]


def process_sample(
    sample_name: str,
) -> int:
    """Create one plot for every cell-type matrix in a replicate.

    Args:
        sample_name: Replicate directory name.

    Returns:
        Number of plots created.

    Raises:
        RuntimeError: If original and binarized matrix sets differ.
    """
    original_matrices = find_sample_matrices(
        root_directory=CELL_TYPE_MATRICES_DIR,
        sample_name=sample_name,
    )

    binary_matrices = find_sample_matrices(
        root_directory=BINARIZED_CELL_TYPE_MATRICES_DIR,
        sample_name=sample_name,
    )

    if not original_matrices and not binary_matrices:
        return 0

    if set(original_matrices) != set(binary_matrices):
        missing_binary = sorted(
            set(original_matrices)
            - set(binary_matrices)
        )

        missing_original = sorted(
            set(binary_matrices)
            - set(original_matrices)
        )

        raise RuntimeError(
            f"{sample_name}: original and binary matrix sets differ.\n"
            f"Missing binary matrices: {missing_binary}\n"
            f"Missing original matrices: {missing_original}"
        )

    n_created = 0
    seen_cell_types: set[str] = set()

    for relative_path in sorted(original_matrices):
        cell_type = get_cell_type(relative_path)

        if cell_type in seen_cell_types:
            raise RuntimeError(
                f"{sample_name}: found more than one matrix for "
                f"cell type {cell_type!r}."
            )

        seen_cell_types.add(cell_type)

        original_path = original_matrices[
            relative_path
        ]
        binary_path = binary_matrices[
            relative_path
        ]

        print(
            f"  Processing {cell_type}"
        )

        (
            original_counts,
            binary_counts,
            matrix_shape,
        ) = load_matrix_pair_distributions(
            original_path=original_path,
            binary_path=binary_path,
        )

        output_path = create_cell_type_plot(
            sample_name=sample_name,
            cell_type=cell_type,
            original_counts=original_counts,
            binary_counts=binary_counts,
            matrix_shape=matrix_shape,
        )

        print(
            f"    Wrote: {output_path}"
        )

        n_created += 1

    return n_created


def main() -> None:
    """Create one binarization-distribution plot per cell type."""
    MATRIX_BINARIZATION_PLOTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    total_created = 0
    skipped_samples: list[str] = []

    for sample_directory in ATAC_SEQ_DIRS:
        sample_name = sample_directory.name

        print(
            f"\nProcessing {sample_name}"
        )

        n_created = process_sample(
            sample_name=sample_name,
        )

        if n_created == 0:
            print(
                "  No split matrices were found; skipping."
            )
            skipped_samples.append(sample_name)
            continue

        total_created += n_created

        print(
            f"  Created {n_created} plots."
        )

    print(
        f"\nCreated {total_created} cell-type plots in "
        f"{MATRIX_BINARIZATION_PLOTS_DIR}."
    )

    if skipped_samples:
        print(
            "Skipped samples with no split matrices: "
            + ", ".join(skipped_samples)
        )


if __name__ == "__main__":
    main()