"""Apply local hybrid peak filtering to binarized cell-type matrices.

For every available sample-by-cell-type matrix, this script applies the chosen
local filtering policy:

    1. A matrix is usable only when it contains at least
       ``MIN_PEAK_FILTER_GROUP_SIZE`` cells.

    2. For a usable matrix with ``n_cells`` columns, a peak is retained when:

           support_cells >= max(
               MIN_PEAK_FILTER_CELL_SUPPORT,
               ceil(MIN_PEAK_FILTER_CELL_FRACTION * n_cells),
           )

       where ``support_cells`` is the number of cells in which the binary
       matrix contains 1 for that peak.

The script preserves the original binarized matrices and writes filtered
outputs beneath ``FILTERED_CELL_TYPE_MATRICES_DIR``. Matrix rows and peak BED
rows are filtered with exactly the same Boolean mask. Cell columns, barcodes,
and cell metadata are preserved unchanged.

The resulting peaks remain sample-specific. Cross-sample consensus alignment
and replicate-level filtering must be performed in a later processing stage.
"""

from __future__ import annotations

import gzip
import hashlib
import math
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

from constants import (
    ATAC_SEQ_DIRS,
    BARCODES_SUFFIX,
    BINARIZED_CELL_TYPE_MATRICES_DIR,
    CELL_TYPE_STANDARDIZATION,
    FILTERED_CELL_TYPE_MATRICES_DIR,
    MATRIX_SUFFIX,
    MIN_PEAK_FILTER_CELL_FRACTION,
    MIN_PEAK_FILTER_CELL_SUPPORT,
    MIN_PEAK_FILTER_GROUP_SIZE,
    PEAK_FILTERING_SUMMARY_PATH,
    PEAKS_SUFFIX,
)


CELL_METADATA_SUFFIX: Final[str] = "_cell_metadata.csv.gz"
PEAK_FILTERING_REPORT_SUFFIX: Final[str] = "_peak_filtering.csv.gz"

PREFERRED_CELL_TYPE_ORDER: Final[tuple[str, ...]] = (
    "Hepatocytes",
    "Endothelial",
    "Cholangiocyte",
    "Kupffer",
    "Stellate",
    "T_NK_B",
    "Unknown",
)


def validate_filtering_parameters(
    minimum_group_size: int,
    minimum_cell_support: int,
    minimum_cell_fraction: float,
) -> None:
    """
    Validate the configured local peak-filtering parameters.

    :param minimum_group_size: Minimum number of cells required for a matrix
        to be eligible.
    :param minimum_cell_support: Minimum absolute number of accessible cells
        required for a peak.
    :param minimum_cell_fraction: Minimum fraction of accessible cells required
        for a peak.
    """
    if minimum_group_size < 1:
        raise ValueError(
            "MIN_PEAK_FILTER_GROUP_SIZE must be a positive integer."
        )

    if minimum_cell_support < 1:
        raise ValueError(
            "MIN_PEAK_FILTER_CELL_SUPPORT must be a positive integer."
        )

    if not 0.0 <= minimum_cell_fraction <= 1.0:
        raise ValueError(
            "MIN_PEAK_FILTER_CELL_FRACTION must be in the interval [0, 1]."
        )


def get_cell_type_order() -> tuple[str, ...]:
    """
    Determine and validate the standardized cell-type order.

    :return: Ordered standardized cell-type names.
    """
    standardized_cell_types = tuple(
        dict.fromkeys(CELL_TYPE_STANDARDIZATION.values())
    )

    if set(standardized_cell_types) != set(PREFERRED_CELL_TYPE_ORDER):
        raise ValueError(
            "CELL_TYPE_STANDARDIZATION contains an unexpected set of "
            f"standardized cell types: {standardized_cell_types}"
        )

    return PREFERRED_CELL_TYPE_ORDER


def infer_condition(sample_name: str) -> str:
    """
    Infer the biological condition from a sample-directory name.

    :param sample_name: Sample-directory name.
    :return: Either ``MASH`` or ``Normal``.
    """
    lowered_name = sample_name.lower()

    if "mash" in lowered_name:
        return "MASH"

    if "normal" in lowered_name:
        return "Normal"

    raise ValueError(
        f"Could not infer condition from sample name: {sample_name}"
    )


def prepare_output_directory() -> None:
    """
    Recreate the derived filtering-output directory to prevent stale outputs.
    """
    input_root = BINARIZED_CELL_TYPE_MATRICES_DIR.resolve()
    output_root = FILTERED_CELL_TYPE_MATRICES_DIR.resolve()

    if input_root == output_root:
        raise RuntimeError(
            "The filtered output directory must differ from the binarized "
            "input directory."
        )

    if FILTERED_CELL_TYPE_MATRICES_DIR.exists():
        shutil.rmtree(FILTERED_CELL_TYPE_MATRICES_DIR)

    FILTERED_CELL_TYPE_MATRICES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


def index_input_matrices(
    cell_type_order: Sequence[str],
) -> dict[tuple[str, str], Path]:
    """
    Index all available binarized matrices by sample and cell type.

    :param cell_type_order: Expected standardized cell-type names.
    :return: Mapping from ``(sample_directory, cell_type)`` to matrix path.
    """
    known_samples = {
        sample_directory.name
        for sample_directory in ATAC_SEQ_DIRS
    }
    known_cell_types = set(cell_type_order)

    matrix_paths = sorted(
        BINARIZED_CELL_TYPE_MATRICES_DIR.rglob(
            f"*{MATRIX_SUFFIX}"
        )
    )

    if not matrix_paths:
        raise FileNotFoundError(
            "No binarized cell-type matrices were found beneath "
            f"{BINARIZED_CELL_TYPE_MATRICES_DIR}."
        )

    indexed_matrices: dict[tuple[str, str], Path] = {}

    for matrix_path in matrix_paths:
        relative_path = matrix_path.relative_to(
            BINARIZED_CELL_TYPE_MATRICES_DIR
        )

        if len(relative_path.parts) != 3:
            raise ValueError(
                "Expected matrix paths with structure "
                "<sample>/<cell_type>/<matrix_file>, but found: "
                f"{matrix_path}"
            )

        sample_name = relative_path.parts[0]
        cell_type = relative_path.parts[1]

        if sample_name not in known_samples:
            raise ValueError(
                f"Unexpected sample directory: {sample_name}"
            )

        if cell_type not in known_cell_types:
            raise ValueError(
                f"Unexpected cell type beneath {sample_name}: {cell_type}"
            )

        key = (sample_name, cell_type)

        if key in indexed_matrices:
            raise RuntimeError(
                f"More than one matrix was found for "
                f"{sample_name} / {cell_type}."
            )

        indexed_matrices[key] = matrix_path

    return indexed_matrices


def get_file_prefix(matrix_path: Path) -> str:
    """
    Remove the configured matrix suffix from a matrix filename.

    :param matrix_path: Matrix path.
    :return: Filename prefix shared by associated files.
    """
    if not matrix_path.name.endswith(MATRIX_SUFFIX):
        raise ValueError(
            f"Unexpected matrix filename: {matrix_path.name}"
        )

    return matrix_path.name[:-len(MATRIX_SUFFIX)]


def get_associated_input_paths(
    matrix_path: Path,
) -> dict[str, Path]:
    """
    Construct and validate associated input-file paths for one matrix.

    :param matrix_path: Binarized matrix path.
    :return: Mapping containing matrix, peaks, barcodes, and metadata paths.
    """
    prefix = get_file_prefix(matrix_path)
    source_directory = matrix_path.parent

    paths = {
        "matrix": matrix_path,
        "peaks": source_directory / f"{prefix}{PEAKS_SUFFIX}",
        "barcodes": source_directory / f"{prefix}{BARCODES_SUFFIX}",
        "metadata": source_directory / f"{prefix}{CELL_METADATA_SUFFIX}",
    }

    missing_paths = [
        path
        for path in paths.values()
        if not path.is_file()
    ]

    if missing_paths:
        formatted_paths = "\n".join(
            f"  - {path}"
            for path in missing_paths
        )
        raise FileNotFoundError(
            "Required associated input files were not found:\n"
            f"{formatted_paths}"
        )

    return paths


def get_output_paths(
    sample_name: str,
    cell_type: str,
    matrix_path: Path,
) -> dict[str, Path]:
    """
    Construct output paths for one eligible sample-by-cell-type group.

    :param sample_name: Sample-directory name.
    :param cell_type: Standardized cell type.
    :param matrix_path: Input matrix path.
    :return: Mapping containing all output paths.
    """
    prefix = get_file_prefix(matrix_path)
    output_directory = (
        FILTERED_CELL_TYPE_MATRICES_DIR
        / sample_name
        / cell_type
    )

    return {
        "directory": output_directory,
        "matrix": output_directory / matrix_path.name,
        "peaks": output_directory / f"{prefix}{PEAKS_SUFFIX}",
        "barcodes": output_directory / f"{prefix}{BARCODES_SUFFIX}",
        "metadata": output_directory / f"{prefix}{CELL_METADATA_SUFFIX}",
        "report": (
            output_directory
            / f"{prefix}{PEAK_FILTERING_REPORT_SUFFIX}"
        ),
    }


def read_binary_matrix(matrix_path: Path) -> sparse.csr_matrix:
    """
    Read and validate one gzip-compressed binary Matrix Market file.

    :param matrix_path: Path to the binarized matrix.
    :return: Validated binary matrix in CSR format.
    """
    with gzip.open(matrix_path, mode="rb") as handle:
        matrix = mmread(handle)

    if not sparse.issparse(matrix):
        raise TypeError(
            f"Expected a sparse matrix in {matrix_path}, "
            f"but found {type(matrix).__name__}."
        )

    matrix = matrix.tocsr()
    raw_nnz = matrix.nnz

    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()

    if matrix.nnz != raw_nnz:
        raise ValueError(
            f"{matrix_path} contained duplicate coordinates or explicit zero "
            f"entries. Stored entries changed from {raw_nnz:,} to "
            f"{matrix.nnz:,} during cleanup."
        )

    if matrix.nnz > 0:
        if not np.isfinite(matrix.data).all():
            raise ValueError(
                f"{matrix_path} contains non-finite values."
            )

        if not np.all(matrix.data == 1):
            unique_values = np.unique(matrix.data)
            raise ValueError(
                f"{matrix_path} is not binary. Stored values are "
                f"{unique_values.tolist()}."
            )

    return matrix


def read_peak_lines(
    peaks_path: Path,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a compressed BED file while preserving every original line.

    :param peaks_path: Path to the gzip-compressed peaks BED file.
    :return: Original BED lines and arrays containing chromosome, start, and
        end coordinates.
    """
    lines: list[str] = []
    chromosomes: list[str] = []
    starts: list[int] = []
    ends: list[int] = []

    with gzip.open(
        peaks_path,
        mode="rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        for row_index, line in enumerate(handle):
            normalized_line = (
                line
                if line.endswith("\n")
                else f"{line}\n"
            )
            stripped_line = normalized_line.rstrip("\r\n")

            if not stripped_line:
                raise ValueError(
                    f"{peaks_path} contains an empty line at row {row_index}."
                )

            fields = stripped_line.split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"{peaks_path} row {row_index} contains fewer than three "
                    "BED columns."
                )

            chromosome = fields[0]

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{peaks_path} row {row_index} contains non-integer "
                    "coordinates."
                ) from error

            if start < 0 or end <= start:
                raise ValueError(
                    f"{peaks_path} row {row_index} has invalid interval "
                    f"{chromosome}:{start}-{end}."
                )

            lines.append(normalized_line)
            chromosomes.append(chromosome)
            starts.append(start)
            ends.append(end)

    return (
        lines,
        np.asarray(chromosomes, dtype=object),
        np.asarray(starts, dtype=np.int64),
        np.asarray(ends, dtype=np.int64),
    )


def count_barcode_rows(barcodes_path: Path) -> int:
    """
    Count and validate barcode rows in a compressed TSV file.

    :param barcodes_path: Path to the gzip-compressed barcode file.
    :return: Number of non-empty barcode rows.
    """
    row_count = 0

    with gzip.open(
        barcodes_path,
        mode="rt",
        encoding="utf-8",
    ) as handle:
        for row_index, line in enumerate(handle):
            if not line.strip():
                raise ValueError(
                    f"{barcodes_path} contains an empty barcode at row "
                    f"{row_index}."
                )

            row_count += 1

    return row_count


def count_metadata_rows(metadata_path: Path) -> int:
    """
    Count cell-metadata records in a compressed CSV file.

    :param metadata_path: Path to the gzip-compressed cell-metadata CSV.
    :return: Number of metadata records excluding the header.
    """
    metadata = pd.read_csv(
        metadata_path,
        compression="gzip",
    )

    return len(metadata)


def calculate_peak_support(
    matrix: sparse.csr_matrix,
    matrix_path: Path,
) -> np.ndarray:
    """
    Calculate the number of accessible cells for every peak row.

    :param matrix: Binary peak-by-cell matrix.
    :param matrix_path: Matrix path used in validation messages.
    :return: Integer support count for every peak row.
    """
    support = np.asarray(
        matrix.sum(axis=1)
    ).ravel().astype(
        np.int64,
        copy=False,
    )

    if support.shape[0] != matrix.shape[0]:
        raise RuntimeError(
            f"{matrix_path}: calculated {support.shape[0]:,} support values "
            f"for {matrix.shape[0]:,} rows."
        )

    if support.size > 0:
        if support.min() < 0 or support.max() > matrix.shape[1]:
            raise RuntimeError(
                f"{matrix_path}: support values are outside the valid range "
                f"[0, {matrix.shape[1]}]."
            )

    if int(support.sum()) != int(matrix.nnz):
        raise RuntimeError(
            f"{matrix_path}: support sum {int(support.sum()):,} does not "
            f"equal matrix nnz {matrix.nnz:,}."
        )

    return support


def calculate_effective_support_threshold(
    n_cells: int,
    minimum_cell_support: int,
    minimum_cell_fraction: float,
) -> tuple[int, int]:
    """
    Calculate the fractional and effective support thresholds for one matrix.

    :param n_cells: Number of cells in the matrix.
    :param minimum_cell_support: Configured absolute support threshold.
    :param minimum_cell_fraction: Configured fractional support threshold.
    :return: Fractional ceiling and final effective support threshold.
    """
    fractional_support_ceiling = math.ceil(
        minimum_cell_fraction * n_cells
    )
    effective_support_threshold = max(
        minimum_cell_support,
        fractional_support_ceiling,
    )

    return (
        fractional_support_ceiling,
        effective_support_threshold,
    )


def validate_input_alignment(
    matrix: sparse.csr_matrix,
    peak_lines: Sequence[str],
    barcode_count: int,
    metadata_count: int,
    matrix_path: Path,
) -> None:
    """
    Validate row and column alignment across all input files.

    :param matrix: Input binary matrix.
    :param peak_lines: BED rows corresponding to matrix rows.
    :param barcode_count: Number of barcode-file rows.
    :param metadata_count: Number of cell-metadata rows.
    :param matrix_path: Matrix path used in validation messages.
    """
    if matrix.shape[0] != len(peak_lines):
        raise RuntimeError(
            f"{matrix_path}: matrix has {matrix.shape[0]:,} rows but the "
            f"peaks file has {len(peak_lines):,} rows."
        )

    if matrix.shape[1] != barcode_count:
        raise RuntimeError(
            f"{matrix_path}: matrix has {matrix.shape[1]:,} columns but the "
            f"barcode file has {barcode_count:,} rows."
        )

    if matrix.shape[1] != metadata_count:
        raise RuntimeError(
            f"{matrix_path}: matrix has {matrix.shape[1]:,} columns but the "
            f"cell-metadata file has {metadata_count:,} rows."
        )


def filter_matrix_rows(
    matrix: sparse.csr_matrix,
    retained_mask: np.ndarray,
) -> sparse.csr_matrix:
    """
    Filter matrix rows while preserving cell columns and retained-row order.

    :param matrix: Input binary matrix.
    :param retained_mask: Boolean row-retention mask.
    :return: Filtered binary matrix in CSR format.
    """
    if retained_mask.dtype != np.bool_:
        raise TypeError(
            "The retained-row mask must have Boolean dtype."
        )

    if retained_mask.shape != (matrix.shape[0],):
        raise ValueError(
            f"Expected a mask of shape {(matrix.shape[0],)}, "
            f"but found {retained_mask.shape}."
        )

    filtered_matrix = matrix[retained_mask, :].tocsr()
    filtered_matrix.sort_indices()

    return filtered_matrix


def write_sparse_matrix(
    matrix: sparse.spmatrix,
    output_path: Path,
) -> None:
    """
    Write a sparse binary matrix as gzip-compressed Matrix Market.

    :param matrix: Sparse matrix to write.
    :param output_path: Destination Matrix Market path.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

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


def write_filtered_peaks(
    peak_lines: Sequence[str],
    retained_mask: np.ndarray,
    output_path: Path,
) -> None:
    """
    Write retained BED rows in their original order.

    :param peak_lines: Original BED lines.
    :param retained_mask: Boolean row-retention mask.
    :param output_path: Destination compressed BED path.
    """
    if len(peak_lines) != retained_mask.shape[0]:
        raise ValueError(
            "Peak-line count and retained-mask length differ."
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with gzip.open(
        output_path,
        mode="wt",
        encoding="utf-8",
        newline="",
        compresslevel=6,
    ) as handle:
        for line, retained in zip(peak_lines, retained_mask):
            if retained:
                handle.write(line)


def build_peak_filtering_report(
    chromosomes: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    support: np.ndarray,
    n_cells: int,
    fractional_support_ceiling: int,
    effective_support_threshold: int,
    retained_mask: np.ndarray,
) -> pd.DataFrame:
    """
    Build a row-level peak-filtering provenance table.

    :param chromosomes: Chromosome value for every original peak row.
    :param starts: Start coordinate for every original peak row.
    :param ends: End coordinate for every original peak row.
    :param support: Accessible-cell count for every original peak row.
    :param n_cells: Number of cells in the matrix.
    :param fractional_support_ceiling: Cell count implied by the fractional
        threshold.
    :param effective_support_threshold: Final support threshold applied.
    :param retained_mask: Boolean row-retention mask.
    :return: Per-peak filtering report.
    """
    n_peaks = support.shape[0]

    filtered_row_values = np.full(
        n_peaks,
        np.nan,
        dtype=float,
    )
    filtered_row_values[retained_mask] = np.arange(
        int(retained_mask.sum()),
        dtype=np.int64,
    )

    report = pd.DataFrame(
        {
            "original_peak_row": np.arange(
                n_peaks,
                dtype=np.int64,
            ),
            "chromosome": chromosomes,
            "start": starts,
            "end": ends,
            "support_cells": support,
            "support_fraction": support / n_cells,
            "minimum_cell_support": (
                MIN_PEAK_FILTER_CELL_SUPPORT
            ),
            "minimum_cell_fraction": (
                MIN_PEAK_FILTER_CELL_FRACTION
            ),
            "fractional_support_ceiling": (
                fractional_support_ceiling
            ),
            "effective_min_support_cells": (
                effective_support_threshold
            ),
            "retained": retained_mask,
            "filtered_peak_row": pd.array(
                filtered_row_values,
                dtype="Int64",
            ),
        }
    )

    return report


def write_peak_filtering_report(
    report: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Write a compressed per-peak filtering report.

    :param report: Per-peak filtering provenance table.
    :param output_path: Destination compressed CSV path.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    report.to_csv(
        output_path,
        index=False,
        compression={
            "method": "gzip",
            "compresslevel": 6,
        },
    )


def calculate_sha256(path: Path) -> str:
    """
    Calculate the SHA-256 checksum of a file.

    :param path: File whose checksum should be calculated.
    :return: Lowercase hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)

    return digest.hexdigest()


def copy_unchanged_file(
    source_path: Path,
    output_path: Path,
) -> None:
    """
    Copy a file and verify that its bytes are unchanged.

    :param source_path: Source file path.
    :param output_path: Destination file path.
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    shutil.copy2(
        source_path,
        output_path,
    )

    source_digest = calculate_sha256(source_path)
    output_digest = calculate_sha256(output_path)

    if source_digest != output_digest:
        raise RuntimeError(
            f"Copied file checksum differs from source:\n"
            f"  Source: {source_path}\n"
            f"  Output: {output_path}"
        )


def validate_written_matrix(
    expected_matrix: sparse.csr_matrix,
    output_path: Path,
) -> None:
    """
    Re-read and validate a written filtered matrix.

    :param expected_matrix: In-memory filtered matrix.
    :param output_path: Written Matrix Market path.
    """
    written_matrix = read_binary_matrix(output_path)

    if written_matrix.shape != expected_matrix.shape:
        raise RuntimeError(
            f"{output_path}: written shape {written_matrix.shape} differs "
            f"from expected shape {expected_matrix.shape}."
        )

    if written_matrix.nnz != expected_matrix.nnz:
        raise RuntimeError(
            f"{output_path}: written nnz {written_matrix.nnz:,} differs "
            f"from expected nnz {expected_matrix.nnz:,}."
        )

    expected_matrix = expected_matrix.copy()
    expected_matrix.sort_indices()

    if not np.array_equal(
        written_matrix.indptr,
        expected_matrix.indptr,
    ):
        raise RuntimeError(
            f"{output_path}: written CSR row structure differs from expected."
        )

    if not np.array_equal(
        written_matrix.indices,
        expected_matrix.indices,
    ):
        raise RuntimeError(
            f"{output_path}: written CSR column indices differ from expected."
        )

    if not np.array_equal(
        written_matrix.data,
        expected_matrix.data,
    ):
        raise RuntimeError(
            f"{output_path}: written matrix values differ from expected."
        )


def validate_written_peaks(
    original_peak_lines: Sequence[str],
    retained_mask: np.ndarray,
    output_path: Path,
) -> None:
    """
    Re-read and validate the written filtered BED file.

    :param original_peak_lines: Original BED lines.
    :param retained_mask: Boolean row-retention mask.
    :param output_path: Written compressed BED path.
    """
    expected_lines = [
        line
        for line, retained in zip(
            original_peak_lines,
            retained_mask,
        )
        if retained
    ]

    with gzip.open(
        output_path,
        mode="rt",
        encoding="utf-8",
        newline="",
    ) as handle:
        written_lines = [
            line
            if line.endswith("\n")
            else f"{line}\n"
            for line in handle
        ]

    if written_lines != expected_lines:
        raise RuntimeError(
            f"{output_path}: written peak rows or their order differ from "
            "the expected retained rows."
        )


def validate_written_report(
    output_path: Path,
    expected_row_count: int,
    expected_retained_count: int,
) -> None:
    """
    Re-read and validate a written peak-filtering report.

    :param output_path: Written compressed CSV path.
    :param expected_row_count: Expected number of original peak rows.
    :param expected_retained_count: Expected number of retained peak rows.
    """
    report = pd.read_csv(
        output_path,
        compression="gzip",
    )

    if len(report) != expected_row_count:
        raise RuntimeError(
            f"{output_path}: report has {len(report):,} rows, expected "
            f"{expected_row_count:,}."
        )

    if int(report["retained"].sum()) != expected_retained_count:
        raise RuntimeError(
            f"{output_path}: report marks "
            f"{int(report['retained'].sum()):,} retained rows, expected "
            f"{expected_retained_count:,}."
        )

    retained_filtered_rows = report.loc[
        report["retained"],
        "filtered_peak_row",
    ].to_numpy(dtype=np.int64)

    expected_filtered_rows = np.arange(
        expected_retained_count,
        dtype=np.int64,
    )

    if not np.array_equal(
        retained_filtered_rows,
        expected_filtered_rows,
    ):
        raise RuntimeError(
            f"{output_path}: filtered row indices are not consecutive and "
            "order-preserving."
        )


def build_missing_summary_record(
    sample_name: str,
    condition: str,
    cell_type: str,
) -> dict[str, Any]:
    """
    Build a summary record for an expected matrix that is unavailable.

    :param sample_name: Sample-directory name.
    :param condition: Biological condition.
    :param cell_type: Standardized cell type.
    :return: Missing-matrix summary record.
    """
    return {
        "sample_directory": sample_name,
        "condition": condition,
        "cell_type": cell_type,
        "status": "missing_matrix",
        "n_cells": None,
        "n_input_peaks": None,
        "minimum_group_size": MIN_PEAK_FILTER_GROUP_SIZE,
        "minimum_cell_support": MIN_PEAK_FILTER_CELL_SUPPORT,
        "minimum_cell_fraction": MIN_PEAK_FILTER_CELL_FRACTION,
        "fractional_support_ceiling": None,
        "effective_min_support_cells": None,
        "n_retained_peaks": None,
        "n_removed_peaks": None,
        "fraction_retained": None,
        "input_nonzero_entries": None,
        "output_nonzero_entries": None,
        "input_matrix_path": None,
        "output_matrix_path": None,
        "input_peaks_path": None,
        "output_peaks_path": None,
        "output_barcodes_path": None,
        "output_metadata_path": None,
        "output_peak_filtering_report_path": None,
    }


def process_available_matrix(
    sample_name: str,
    condition: str,
    cell_type: str,
    matrix_path: Path,
) -> dict[str, Any]:
    """
    Validate, filter, write, and summarize one available matrix.

    :param sample_name: Sample-directory name.
    :param condition: Biological condition.
    :param cell_type: Standardized cell type.
    :param matrix_path: Binarized input matrix path.
    :return: One matrix-level filtering summary record.
    """
    input_paths = get_associated_input_paths(matrix_path)
    matrix = read_binary_matrix(input_paths["matrix"])

    peak_lines, chromosomes, starts, ends = read_peak_lines(
        input_paths["peaks"]
    )
    barcode_count = count_barcode_rows(
        input_paths["barcodes"]
    )
    metadata_count = count_metadata_rows(
        input_paths["metadata"]
    )

    validate_input_alignment(
        matrix=matrix,
        peak_lines=peak_lines,
        barcode_count=barcode_count,
        metadata_count=metadata_count,
        matrix_path=matrix_path,
    )

    n_input_peaks, n_cells = matrix.shape

    (
        fractional_support_ceiling,
        effective_support_threshold,
    ) = calculate_effective_support_threshold(
        n_cells=n_cells,
        minimum_cell_support=MIN_PEAK_FILTER_CELL_SUPPORT,
        minimum_cell_fraction=MIN_PEAK_FILTER_CELL_FRACTION,
    )

    common_summary: dict[str, Any] = {
        "sample_directory": sample_name,
        "condition": condition,
        "cell_type": cell_type,
        "n_cells": n_cells,
        "n_input_peaks": n_input_peaks,
        "minimum_group_size": MIN_PEAK_FILTER_GROUP_SIZE,
        "minimum_cell_support": MIN_PEAK_FILTER_CELL_SUPPORT,
        "minimum_cell_fraction": MIN_PEAK_FILTER_CELL_FRACTION,
        "fractional_support_ceiling": fractional_support_ceiling,
        "effective_min_support_cells": effective_support_threshold,
        "input_nonzero_entries": matrix.nnz,
        "input_matrix_path": str(input_paths["matrix"]),
        "input_peaks_path": str(input_paths["peaks"]),
    }

    if n_cells < MIN_PEAK_FILTER_GROUP_SIZE:
        print(
            f"  Excluded: {n_cells:,} cells < "
            f"{MIN_PEAK_FILTER_GROUP_SIZE:,}"
        )

        return {
            **common_summary,
            "status": "excluded_group_too_small",
            "n_retained_peaks": None,
            "n_removed_peaks": None,
            "fraction_retained": None,
            "output_nonzero_entries": None,
            "output_matrix_path": None,
            "output_peaks_path": None,
            "output_barcodes_path": None,
            "output_metadata_path": None,
            "output_peak_filtering_report_path": None,
        }

    support = calculate_peak_support(
        matrix=matrix,
        matrix_path=matrix_path,
    )
    retained_mask = (
        support >= effective_support_threshold
    )
    filtered_matrix = filter_matrix_rows(
        matrix=matrix,
        retained_mask=retained_mask,
    )

    n_retained_peaks = int(retained_mask.sum())
    n_removed_peaks = n_input_peaks - n_retained_peaks
    fraction_retained = (
        n_retained_peaks / n_input_peaks
        if n_input_peaks > 0
        else 0.0
    )

    if n_retained_peaks == 0:
        raise RuntimeError(
            f"{matrix_path}: the configured rule removed every peak."
        )

    if filtered_matrix.shape != (
        n_retained_peaks,
        n_cells,
    ):
        raise RuntimeError(
            f"{matrix_path}: filtered matrix shape "
            f"{filtered_matrix.shape} differs from expected "
            f"{(n_retained_peaks, n_cells)}."
        )

    if filtered_matrix.nnz != int(
        support[retained_mask].sum()
    ):
        raise RuntimeError(
            f"{matrix_path}: filtered matrix nnz does not equal retained "
            "support sum."
        )

    output_paths = get_output_paths(
        sample_name=sample_name,
        cell_type=cell_type,
        matrix_path=matrix_path,
    )
    output_paths["directory"].mkdir(
        parents=True,
        exist_ok=True,
    )

    write_sparse_matrix(
        matrix=filtered_matrix,
        output_path=output_paths["matrix"],
    )
    write_filtered_peaks(
        peak_lines=peak_lines,
        retained_mask=retained_mask,
        output_path=output_paths["peaks"],
    )

    copy_unchanged_file(
        source_path=input_paths["barcodes"],
        output_path=output_paths["barcodes"],
    )
    copy_unchanged_file(
        source_path=input_paths["metadata"],
        output_path=output_paths["metadata"],
    )

    report = build_peak_filtering_report(
        chromosomes=chromosomes,
        starts=starts,
        ends=ends,
        support=support,
        n_cells=n_cells,
        fractional_support_ceiling=fractional_support_ceiling,
        effective_support_threshold=effective_support_threshold,
        retained_mask=retained_mask,
    )
    write_peak_filtering_report(
        report=report,
        output_path=output_paths["report"],
    )

    validate_written_matrix(
        expected_matrix=filtered_matrix,
        output_path=output_paths["matrix"],
    )
    validate_written_peaks(
        original_peak_lines=peak_lines,
        retained_mask=retained_mask,
        output_path=output_paths["peaks"],
    )
    validate_written_report(
        output_path=output_paths["report"],
        expected_row_count=n_input_peaks,
        expected_retained_count=n_retained_peaks,
    )

    written_barcode_count = count_barcode_rows(
        output_paths["barcodes"]
    )
    written_metadata_count = count_metadata_rows(
        output_paths["metadata"]
    )

    if written_barcode_count != n_cells:
        raise RuntimeError(
            f"{output_paths['barcodes']}: written barcode count differs "
            f"from matrix column count {n_cells:,}."
        )

    if written_metadata_count != n_cells:
        raise RuntimeError(
            f"{output_paths['metadata']}: written metadata count differs "
            f"from matrix column count {n_cells:,}."
        )

    print(
        f"  Retained {n_retained_peaks:,} / "
        f"{n_input_peaks:,} peaks "
        f"({fraction_retained:.1%}); "
        f"effective support >= {effective_support_threshold}"
    )

    return {
        **common_summary,
        "status": "written",
        "n_retained_peaks": n_retained_peaks,
        "n_removed_peaks": n_removed_peaks,
        "fraction_retained": fraction_retained,
        "output_nonzero_entries": filtered_matrix.nnz,
        "output_matrix_path": str(output_paths["matrix"]),
        "output_peaks_path": str(output_paths["peaks"]),
        "output_barcodes_path": str(output_paths["barcodes"]),
        "output_metadata_path": str(output_paths["metadata"]),
        "output_peak_filtering_report_path": str(
            output_paths["report"]
        ),
    }


def validate_global_summary(
    summary: pd.DataFrame,
    expected_group_count: int,
) -> None:
    """
    Validate the completed global peak-filtering summary.

    :param summary: Matrix-level filtering summary.
    :param expected_group_count: Expected number of sample-by-cell-type groups.
    """
    if len(summary) != expected_group_count:
        raise RuntimeError(
            f"Summary contains {len(summary):,} rows, expected "
            f"{expected_group_count:,}."
        )

    allowed_statuses = {
        "written",
        "excluded_group_too_small",
        "missing_matrix",
    }
    observed_statuses = set(summary["status"])

    if not observed_statuses.issubset(allowed_statuses):
        raise RuntimeError(
            f"Unexpected summary statuses: "
            f"{sorted(observed_statuses - allowed_statuses)}"
        )

    written = summary.loc[
        summary["status"] == "written"
    ]

    if written.empty:
        raise RuntimeError(
            "No filtered matrices were written."
        )

    if (
        written["n_retained_peaks"]
        > written["n_input_peaks"]
    ).any():
        raise RuntimeError(
            "At least one matrix retained more peaks than it contained."
        )

    if (
        written["n_removed_peaks"]
        != (
            written["n_input_peaks"]
            - written["n_retained_peaks"]
        )
    ).any():
        raise RuntimeError(
            "At least one summary row has an inconsistent removed-peak count."
        )

    expected_fraction = (
        written["n_retained_peaks"]
        / written["n_input_peaks"]
    )

    if not np.allclose(
        written["fraction_retained"],
        expected_fraction,
    ):
        raise RuntimeError(
            "At least one summary row has an inconsistent retained fraction."
        )

    excluded = summary.loc[
        summary["status"]
        == "excluded_group_too_small"
    ]

    if (
        excluded["n_cells"]
        >= MIN_PEAK_FILTER_GROUP_SIZE
    ).any():
        raise RuntimeError(
            "At least one excluded group meets the configured group-size "
            "threshold."
        )


def print_global_summary(summary: pd.DataFrame) -> None:
    """
    Print overall and condition-specific filtering totals.

    :param summary: Matrix-level filtering summary.
    """
    status_counts = summary["status"].value_counts()

    print("\nPeak filtering completed successfully.")
    print(
        f"Written matrices: "
        f"{int(status_counts.get('written', 0)):,}"
    )
    print(
        f"Groups excluded as too small: "
        f"{int(status_counts.get('excluded_group_too_small', 0)):,}"
    )
    print(
        f"Missing matrices: "
        f"{int(status_counts.get('missing_matrix', 0)):,}"
    )

    written = summary.loc[
        summary["status"] == "written"
    ]

    for condition in ("MASH", "Normal"):
        condition_rows = written.loc[
            written["condition"] == condition
        ]

        total_input_peaks = int(
            condition_rows["n_input_peaks"].sum()
        )
        total_retained_peaks = int(
            condition_rows["n_retained_peaks"].sum()
        )
        pooled_fraction = (
            total_retained_peaks / total_input_peaks
            if total_input_peaks > 0
            else 0.0
        )

        print(
            f"{condition}: retained "
            f"{total_retained_peaks:,} / "
            f"{total_input_peaks:,} local peak rows "
            f"({pooled_fraction:.1%})."
        )

    print(
        f"Summary written to: {PEAK_FILTERING_SUMMARY_PATH}"
    )


def main() -> None:
    """
    Apply local hybrid peak filtering to all expected sample-cell-type groups.
    """
    validate_filtering_parameters(
        minimum_group_size=MIN_PEAK_FILTER_GROUP_SIZE,
        minimum_cell_support=MIN_PEAK_FILTER_CELL_SUPPORT,
        minimum_cell_fraction=MIN_PEAK_FILTER_CELL_FRACTION,
    )

    cell_type_order = get_cell_type_order()
    sample_order = tuple(
        sample_directory.name
        for sample_directory in ATAC_SEQ_DIRS
    )
    expected_group_count = (
        len(sample_order) * len(cell_type_order)
    )

    indexed_matrices = index_input_matrices(
        cell_type_order=cell_type_order,
    )

    prepare_output_directory()

    summary_records: list[dict[str, Any]] = []

    for sample_name in sample_order:
        condition = infer_condition(sample_name)

        for cell_type in cell_type_order:
            print(
                f"Processing {sample_name} / {cell_type}"
            )

            matrix_path = indexed_matrices.get(
                (sample_name, cell_type)
            )

            if matrix_path is None:
                print("  Missing matrix.")
                summary_records.append(
                    build_missing_summary_record(
                        sample_name=sample_name,
                        condition=condition,
                        cell_type=cell_type,
                    )
                )
                continue

            summary_records.append(
                process_available_matrix(
                    sample_name=sample_name,
                    condition=condition,
                    cell_type=cell_type,
                    matrix_path=matrix_path,
                )
            )

    summary = pd.DataFrame.from_records(
        summary_records
    )

    validate_global_summary(
        summary=summary,
        expected_group_count=expected_group_count,
    )

    summary.to_csv(
        PEAK_FILTERING_SUMMARY_PATH,
        index=False,
    )

    print_global_summary(
        summary=summary,
    )


if __name__ == "__main__":
    main()
