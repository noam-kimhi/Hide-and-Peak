#!/usr/bin/env python3
"""Create eQTL-overlap and non-eQTL peak sets for filtered cell-type matrices.

For every ``*_peaks.bed.gz`` file under the filtered cell-type matrix
directory, this script writes:

1. ``*_peaks_without_eqt.bed.gz``
   Complete original peaks containing no exact GTEx liver eQTL variant.
2. ``*_peaks_eqt_overlap.bed.gz``
   Complete original peaks containing at least one exact GTEx liver eQTL
   variant.
3. ``*_peaks_eqt_overlap_annotation.csv.gz``
   One row per original peak, including its original zero-based matrix row
   index and detailed eQTL overlap annotations.

The eQTL variant position is converted from the one-based GTEx variant ID to a
one-base, zero-based, half-open BED interval: ``[position - 1, position)``.

The original peak coordinates are never split or modified. The script uses
``bioframe.count_overlaps`` to classify peaks and ``bioframe.overlap`` to
collect the overlapping variant and gene identifiers.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Iterable, Sequence

import bioframe
import pandas as pd

from constants import (
    EQTL_PEAK_OVERLAP_SUMMARY_PATH,
    FILTERED_CELL_TYPE_MATRICES_DIR,
    LIVER_EQTL_SIGNIFICANT_PAIRS_PATH,
    PEAKS_EQT_OVERLAP_ANNOTATION_SUFFIX,
    PEAKS_EQT_OVERLAP_SUFFIX,
    PEAKS_SUFFIX,
    PEAKS_WITHOUT_EQT_SUFFIX,
)

LOGGER = logging.getLogger(__name__)

SCRIPT_VERSION = "2026-06-23.2"
BIOFRAME_PEAK_SUFFIX = "_peak"
BIOFRAME_EQTL_SUFFIX = "_eqt"
BIOFRAME_PEAK_INDEX_COLUMN = f"peak_index{BIOFRAME_PEAK_SUFFIX}"
BIOFRAME_VARIANT_ID_COLUMN = f"variant_id{BIOFRAME_EQTL_SUFFIX}"
BIOFRAME_GENE_IDS_COLUMN = f"gene_ids{BIOFRAME_EQTL_SUFFIX}"

BED_COORDINATE_COLUMNS: tuple[str, str, str] = ("chrom", "start", "end")
EQTL_REQUIRED_COLUMNS: tuple[str, str] = ("variant_id", "gene_id")
EXPECTED_EQTL_BUILD = "b38"
REPLICATE_PATTERN = re.compile(r"(rep\d+)", flags=re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Split every filtered peak BED file into exact liver-eQTL-overlap "
            "and non-overlap peak sets, and write a full annotation table."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=FILTERED_CELL_TYPE_MATRICES_DIR,
        help=(
            "Root containing sample/cell-type directories. "
            f"Default: {FILTERED_CELL_TYPE_MATRICES_DIR}"
        ),
    )
    parser.add_argument(
        "--eqtl-path",
        type=Path,
        default=LIVER_EQTL_SIGNIFICANT_PAIRS_PATH,
        help=(
            "GTEx liver significant variant-gene-pair table. "
            f"Default: {LIVER_EQTL_SIGNIFICANT_PAIRS_PATH}"
        ),
    )
    parser.add_argument(
        "--summary-path",
        type=Path,
        default=EQTL_PEAK_OVERLAP_SUMMARY_PATH,
        help=(
            "Dataset-wide output summary CSV. "
            f"Default: {EQTL_PEAK_OVERLAP_SUMMARY_PATH}"
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing generated outputs.",
    )
    return parser.parse_args()


def normalize_chromosome_name(value: object) -> str:
    """Return a chromosome name in ``chr*`` form for interval matching.

    The returned name is used only for Bioframe operations. Original chromosome
    labels are retained in all output BED and annotation files.
    """
    chrom = str(value).strip()
    if not chrom:
        raise ValueError("Encountered an empty chromosome name.")

    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    if chrom.upper() == "MT":
        chrom = "M"
    elif chrom.upper() in {"X", "Y", "M"}:
        chrom = chrom.upper()
    elif chrom.isdigit():
        chrom = str(int(chrom))

    return f"chr{chrom}"


def join_unique(values: Iterable[object], separator: str = ";") -> str:
    """Join unique, non-empty values in deterministic lexical order."""
    unique_values = {
        str(value).strip()
        for value in values
        if pd.notna(value) and str(value).strip()
    }
    return separator.join(sorted(unique_values))


def join_unique_tokens(values: Iterable[object], separator: str = ";") -> str:
    """Flatten separator-delimited values and join their unique tokens."""
    tokens: set[str] = set()
    for value in values:
        if pd.isna(value):
            continue
        tokens.update(
            token.strip()
            for token in str(value).split(separator)
            if token.strip()
        )
    return separator.join(sorted(tokens))


def load_eqtl_variants(eqtl_path: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    """Load and aggregate significant liver eQTL variants.

    One row is returned per unique ``variant_id``. A variant associated with
    several genes remains one interval, while all of its associated genes are
    retained in the ``gene_ids`` column.
    """
    if not eqtl_path.is_file():
        raise FileNotFoundError(f"eQTL file does not exist: {eqtl_path}")

    LOGGER.info("Reading GTEx liver eQTL pairs from %s", eqtl_path)
    pairs = pd.read_csv(
        eqtl_path,
        sep="\t",
        compression="infer",
        usecols=list(EQTL_REQUIRED_COLUMNS),
        dtype={"variant_id": "string", "gene_id": "string"},
    )

    missing_columns = set(EQTL_REQUIRED_COLUMNS) - set(pairs.columns)
    if missing_columns:
        raise ValueError(
            "The eQTL table is missing required columns: "
            f"{sorted(missing_columns)}"
        )
    if pairs.empty:
        raise ValueError(f"The eQTL table contains no rows: {eqtl_path}")
    if pairs[list(EQTL_REQUIRED_COLUMNS)].isna().any().any():
        missing_counts = (
            pairs[list(EQTL_REQUIRED_COLUMNS)]
            .isna()
            .sum()
            .to_dict()
        )
        raise ValueError(
            "The eQTL table contains missing required values: "
            f"{missing_counts}"
        )

    pairs = pairs.drop_duplicates(subset=list(EQTL_REQUIRED_COLUMNS))
    n_variant_gene_pairs = len(pairs)

    genes_by_variant = (
        pairs.groupby("variant_id", sort=False, observed=True)["gene_id"]
        .agg(join_unique)
        .rename("gene_ids")
        .reset_index()
    )

    variants = genes_by_variant[["variant_id"]].copy()
    parts = variants["variant_id"].str.rsplit("_", n=4, expand=True)
    if parts.shape[1] != 5:
        raise ValueError(
            "Could not parse GTEx variant IDs as "
            "'chrom_position_ref_alt_build'."
        )

    parts.columns = ["chrom", "position", "ref", "alt", "build"]
    malformed = parts.isna().any(axis=1)
    if malformed.any():
        examples = variants.loc[malformed, "variant_id"].head(5).tolist()
        raise ValueError(
            "Malformed GTEx variant IDs were found. Examples: "
            f"{examples}"
        )

    positions = pd.to_numeric(parts["position"], errors="coerce")
    invalid_position = positions.isna() | (positions < 1)
    if invalid_position.any():
        examples = variants.loc[invalid_position, "variant_id"].head(5).tolist()
        raise ValueError(
            "GTEx variant IDs contain invalid one-based positions. Examples: "
            f"{examples}"
        )

    builds = parts["build"].str.lower()
    unexpected_builds = sorted(set(builds) - {EXPECTED_EQTL_BUILD})
    if unexpected_builds:
        raise ValueError(
            "Expected only GRCh38/b38 GTEx variants, but found builds: "
            f"{unexpected_builds}"
        )

    variants = pd.concat(
        [
            variants,
            parts[["chrom", "ref", "alt", "build"]],
            positions.astype("int64").rename("position"),
        ],
        axis=1,
    )
    variants = variants.merge(
        genes_by_variant,
        on="variant_id",
        how="left",
        validate="one_to_one",
    )

    variants["chrom"] = variants["chrom"].map(normalize_chromosome_name)
    variants["start"] = variants["position"] - 1
    variants["end"] = variants["position"]

    variants = variants[
        [
            "chrom",
            "start",
            "end",
            "variant_id",
            "gene_ids",
            "position",
            "ref",
            "alt",
            "build",
        ]
    ].sort_values(
        ["chrom", "start", "end", "variant_id"],
        kind="stable",
        ignore_index=True,
    )

    statistics = {
        "n_variant_gene_pairs": n_variant_gene_pairs,
        "n_unique_variants": int(variants["variant_id"].nunique()),
        "n_unique_genes": int(pairs["gene_id"].nunique()),
    }

    LOGGER.info(
        "Loaded %d unique variants across %d variant-gene pairs and %d genes",
        statistics["n_unique_variants"],
        statistics["n_variant_gene_pairs"],
        statistics["n_unique_genes"],
    )
    return variants, statistics


def read_peak_bed(peak_path: Path) -> pd.DataFrame:
    """Read a headerless BED-like peak file while preserving extra columns."""
    try:
        peaks = pd.read_csv(
            peak_path,
            sep="\t",
            compression="infer",
            header=None,
            comment="#",
            dtype={0: "string"},
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=list(BED_COORDINATE_COLUMNS))

    if peaks.shape[1] < 3:
        raise ValueError(
            f"Peak file must have at least three BED columns: {peak_path}"
        )

    extra_columns = [
        f"bed_column_{column_number}"
        for column_number in range(4, peaks.shape[1] + 1)
    ]
    peaks.columns = [*BED_COORDINATE_COLUMNS, *extra_columns]

    peaks["chrom"] = peaks["chrom"].astype("string")
    peaks["start"] = pd.to_numeric(peaks["start"], errors="coerce")
    peaks["end"] = pd.to_numeric(peaks["end"], errors="coerce")

    invalid = (
        peaks[list(BED_COORDINATE_COLUMNS)].isna().any(axis=1)
        | (peaks["start"] < 0)
        | (peaks["end"] <= peaks["start"])
        | (peaks["start"] % 1 != 0)
        | (peaks["end"] % 1 != 0)
    )
    if invalid.any():
        bad_rows = peaks.loc[invalid].head(5).to_dict(orient="records")
        raise ValueError(
            f"Invalid BED intervals in {peak_path}. Examples: {bad_rows}"
        )

    peaks["start"] = peaks["start"].astype("int64")
    peaks["end"] = peaks["end"].astype("int64")
    return peaks


def make_overlap_peak_frame(peaks: pd.DataFrame) -> pd.DataFrame:
    """Create a normalized BedFrame for Bioframe without changing output data."""
    overlap_peaks = peaks.copy()
    overlap_peaks.insert(0, "peak_index", range(len(overlap_peaks)))
    overlap_peaks["chrom"] = overlap_peaks["chrom"].map(
        normalize_chromosome_name
    )
    return overlap_peaks


def build_output_paths(peak_path: Path) -> tuple[Path, Path, Path]:
    """Construct the three output paths for one input peak BED."""
    if not peak_path.name.endswith(PEAKS_SUFFIX):
        raise ValueError(
            f"Input filename does not end with {PEAKS_SUFFIX}: {peak_path}"
        )

    prefix = peak_path.name[: -len(PEAKS_SUFFIX)]
    without_eqt_path = peak_path.with_name(
        f"{prefix}{PEAKS_WITHOUT_EQT_SUFFIX}"
    )
    eqt_overlap_path = peak_path.with_name(
        f"{prefix}{PEAKS_EQT_OVERLAP_SUFFIX}"
    )
    annotation_path = peak_path.with_name(
        f"{prefix}{PEAKS_EQT_OVERLAP_ANNOTATION_SUFFIX}"
    )
    return without_eqt_path, eqt_overlap_path, annotation_path


def write_dataframe_atomic(
    dataframe: pd.DataFrame,
    output_path: Path,
    *,
    header: bool,
    index: bool,
) -> None:
    """Write a DataFrame atomically, with gzip selected from the final suffix."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f".{output_path.name}.tmp")

    compression = "gzip" if output_path.name.endswith(".gz") else None
    dataframe.to_csv(
        temporary_path,
        sep="\t" if ".bed" in output_path.name else ",",
        header=header,
        index=index,
        compression=compression,
    )
    temporary_path.replace(output_path)


def normalize_bioframe_overlap_columns(
    overlap_pairs: pd.DataFrame,
) -> pd.DataFrame:
    """Rename Bioframe's suffixed overlap columns to stable internal names.

    With ``suffixes=("_peak", "_eqt")``, current Bioframe versions return
    ``peak_index_peak``, ``variant_id_eqt`` and ``gene_ids_eqt``. This helper
    also accepts unsuffixed names for compatibility, but always returns the
    stable names ``peak_index``, ``variant_id`` and ``gene_ids``.
    """
    aliases: dict[str, tuple[str, ...]] = {
        "peak_index": (BIOFRAME_PEAK_INDEX_COLUMN, "peak_index"),
        "variant_id": (BIOFRAME_VARIANT_ID_COLUMN, "variant_id"),
        "gene_ids": (BIOFRAME_GENE_IDS_COLUMN, "gene_ids"),
    }

    rename_map: dict[str, str] = {}
    for target_name, candidates in aliases.items():
        source_name = next(
            (name for name in candidates if name in overlap_pairs.columns),
            None,
        )
        if source_name is None:
            raise RuntimeError(
                "Unexpected columns returned by bioframe.overlap; could not "
                f"resolve {target_name!r}. Expected one of {list(candidates)}. "
                f"Available columns: {list(overlap_pairs.columns)}"
            )
        if source_name != target_name:
            rename_map[source_name] = target_name

    normalized = overlap_pairs.rename(columns=rename_map).copy()
    required = {"peak_index", "variant_id", "gene_ids"}
    missing = required - set(normalized.columns)
    if missing:
        raise RuntimeError(
            "Failed to normalize bioframe.overlap columns; missing "
            f"{sorted(missing)}. Available columns: {list(normalized.columns)}"
        )
    return normalized


def aggregate_overlap_details(overlap_pairs: pd.DataFrame) -> pd.DataFrame:
    """Aggregate Bioframe peak-variant pairs to one row per original peak."""
    if overlap_pairs.empty:
        return pd.DataFrame(
            columns=[
                "peak_index",
                "n_overlapping_eqt_variants",
                "eqt_variant_ids",
                "eqt_gene_ids",
            ]
        )

    normalized = normalize_bioframe_overlap_columns(overlap_pairs)
    return (
        normalized.groupby("peak_index", sort=True, observed=True)
        .agg(
            n_overlapping_eqt_variants=("variant_id", "nunique"),
            eqt_variant_ids=("variant_id", join_unique),
            eqt_gene_ids=("gene_ids", join_unique_tokens),
        )
        .reset_index()
    )

def infer_sample_metadata(
    peak_path: Path,
    input_dir: Path,
) -> dict[str, str]:
    """Infer sample, condition, replicate and cell type from the directory tree."""
    relative_path = peak_path.relative_to(input_dir)
    parts = relative_path.parts

    sample = parts[0] if len(parts) >= 3 else peak_path.parent.parent.name
    cell_type = parts[-2] if len(parts) >= 2 else peak_path.parent.name

    sample_lower = sample.lower()
    if "mash" in sample_lower or "masld" in sample_lower:
        condition = "MASLD"
    elif "normal" in sample_lower or "healthy" in sample_lower:
        condition = "Normal"
    else:
        condition = "Unknown"

    replicate_match = REPLICATE_PATTERN.search(sample)
    replicate = replicate_match.group(1) if replicate_match else ""

    return {
        "sample": sample,
        "condition": condition,
        "replicate": replicate,
        "cell_type": cell_type,
    }


def process_peak_file(
    peak_path: Path,
    input_dir: Path,
    eqtl_variants: pd.DataFrame,
    *,
    overwrite: bool,
) -> dict[str, object]:
    """Create all three outputs for one peak BED and return summary metrics."""
    without_eqt_path, eqt_overlap_path, annotation_path = build_output_paths(
        peak_path
    )
    output_paths = (without_eqt_path, eqt_overlap_path, annotation_path)

    existing_outputs = [path for path in output_paths if path.exists()]
    if existing_outputs and not overwrite:
        raise FileExistsError(
            "Generated outputs already exist. Re-run with --overwrite to "
            f"replace them: {existing_outputs}"
        )

    LOGGER.info("Processing %s", peak_path)
    peaks = read_peak_bed(peak_path)
    peak_columns = list(peaks.columns)

    if peaks.empty:
        empty_annotation = peaks.copy()
        empty_annotation.insert(0, "peak_index", pd.Series(dtype="int64"))
        empty_annotation["overlaps_eqt"] = pd.Series(dtype="bool")
        empty_annotation["n_overlapping_eqt_variants"] = pd.Series(
            dtype="int64"
        )
        empty_annotation["eqt_variant_ids"] = pd.Series(dtype="string")
        empty_annotation["eqt_gene_ids"] = pd.Series(dtype="string")

        write_dataframe_atomic(
            peaks, without_eqt_path, header=False, index=False
        )
        write_dataframe_atomic(
            peaks, eqt_overlap_path, header=False, index=False
        )
        write_dataframe_atomic(
            empty_annotation, annotation_path, header=True, index=False
        )

        metadata = infer_sample_metadata(peak_path, input_dir)
        return {
            **metadata,
            "input_peak_path": str(peak_path.relative_to(input_dir)),
            "n_input_peaks": 0,
            "n_peaks_without_eqt": 0,
            "n_peaks_with_eqt_overlap": 0,
            "fraction_peaks_with_eqt_overlap": 0.0,
            "n_distinct_overlapping_eqt_variants": 0,
            "n_distinct_overlapping_eqt_genes": 0,
        }

    overlap_peaks = make_overlap_peak_frame(peaks)
    peak_chromosomes = set(overlap_peaks["chrom"])
    relevant_eqtl_variants = eqtl_variants[
        eqtl_variants["chrom"].isin(peak_chromosomes)
    ].copy()

    shared_chromosomes = peak_chromosomes & set(relevant_eqtl_variants["chrom"])
    if not shared_chromosomes:
        raise ValueError(
            "No chromosome names are shared between peaks and eQTL variants "
            f"after normalization for {peak_path}."
        )

    counted = bioframe.count_overlaps(
        overlap_peaks.copy(),
        relevant_eqtl_variants,
        return_input=True,
    )
    if "count" not in counted.columns:
        raise RuntimeError(
            "bioframe.count_overlaps did not return the expected 'count' column."
        )
    if len(counted) != len(peaks):
        raise RuntimeError(
            "bioframe.count_overlaps changed the number of peak rows for "
            f"{peak_path}: expected {len(peaks)}, observed {len(counted)}."
        )

    counts = (
        counted[["peak_index", "count"]]
        .sort_values("peak_index", kind="stable")
        .reset_index(drop=True)
    )
    expected_indices = pd.Series(range(len(peaks)), name="peak_index")
    if not counts["peak_index"].reset_index(drop=True).equals(expected_indices):
        raise RuntimeError(
            "Peak order/index was not preserved by bioframe.count_overlaps "
            f"for {peak_path}."
        )

    overlap_pairs = bioframe.overlap(
        overlap_peaks,
        relevant_eqtl_variants,
        how="inner",
        return_input=True,
        return_index=False,
        suffixes=(BIOFRAME_PEAK_SUFFIX, BIOFRAME_EQTL_SUFFIX),
    )
    overlap_details = aggregate_overlap_details(overlap_pairs)

    annotation = peaks.copy()
    annotation.insert(0, "peak_index", range(len(annotation)))
    annotation["n_overlapping_eqt_variants"] = counts["count"].astype("int64")
    annotation["overlaps_eqt"] = (
        annotation["n_overlapping_eqt_variants"] > 0
    )
    annotation = annotation.merge(
        overlap_details,
        on="peak_index",
        how="left",
        suffixes=("", "_details"),
        validate="one_to_one",
        sort=False,
    )

    detail_count_column = "n_overlapping_eqt_variants_details"
    has_detail = annotation[detail_count_column].notna()
    if not (
        annotation.loc[has_detail, "n_overlapping_eqt_variants"]
        .astype("int64")
        .equals(
            annotation.loc[has_detail, detail_count_column].astype("int64")
        )
    ):
        raise RuntimeError(
            "Overlap counts from bioframe.count_overlaps and bioframe.overlap "
            f"disagree for {peak_path}."
        )

    annotation = annotation.drop(columns=[detail_count_column])
    annotation["eqt_variant_ids"] = (
        annotation["eqt_variant_ids"].fillna("").astype("string")
    )
    annotation["eqt_gene_ids"] = (
        annotation["eqt_gene_ids"].fillna("").astype("string")
    )

    annotation_columns = [
        "peak_index",
        *peak_columns,
        "overlaps_eqt",
        "n_overlapping_eqt_variants",
        "eqt_variant_ids",
        "eqt_gene_ids",
    ]
    annotation = annotation[annotation_columns]

    no_eqt_mask = ~annotation["overlaps_eqt"]
    eqt_mask = annotation["overlaps_eqt"]

    peaks_without_eqt = peaks.loc[no_eqt_mask.to_numpy()].copy()
    peaks_with_eqt = peaks.loc[eqt_mask.to_numpy()].copy()

    if len(peaks_without_eqt) + len(peaks_with_eqt) != len(peaks):
        raise RuntimeError(
            f"Peak partition is incomplete for {peak_path}."
        )

    write_dataframe_atomic(
        peaks_without_eqt,
        without_eqt_path,
        header=False,
        index=False,
    )
    write_dataframe_atomic(
        peaks_with_eqt,
        eqt_overlap_path,
        header=False,
        index=False,
    )
    write_dataframe_atomic(
        annotation,
        annotation_path,
        header=True,
        index=False,
    )

    overlapping_variant_ids: set[str] = set()
    for variant_ids in overlap_details["eqt_variant_ids"].dropna():
        overlapping_variant_ids.update(
            variant_id.strip()
            for variant_id in str(variant_ids).split(";")
            if variant_id.strip()
        )

    overlapping_gene_ids: set[str] = set()
    for gene_ids in overlap_details["eqt_gene_ids"].dropna():
        overlapping_gene_ids.update(
            gene_id.strip()
            for gene_id in str(gene_ids).split(";")
            if gene_id.strip()
        )

    n_input_peaks = len(peaks)
    n_with_eqt = len(peaks_with_eqt)
    metadata = infer_sample_metadata(peak_path, input_dir)
    summary = {
        **metadata,
        "input_peak_path": str(peak_path.relative_to(input_dir)),
        "n_input_peaks": n_input_peaks,
        "n_peaks_without_eqt": len(peaks_without_eqt),
        "n_peaks_with_eqt_overlap": n_with_eqt,
        "fraction_peaks_with_eqt_overlap": n_with_eqt / n_input_peaks,
        "n_distinct_overlapping_eqt_variants": len(
            overlapping_variant_ids
        ),
        "n_distinct_overlapping_eqt_genes": len(overlapping_gene_ids),
    }

    LOGGER.info(
        "Finished %s: %d / %d peaks (%.2f%%) overlap exact eQTL variants",
        peak_path.name,
        n_with_eqt,
        n_input_peaks,
        100.0 * summary["fraction_peaks_with_eqt_overlap"],
    )
    return summary


def discover_peak_files(input_dir: Path) -> list[Path]:
    """Return all original filtered peak BED files recursively."""
    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Filtered cell-type matrix directory does not exist: {input_dir}"
        )

    peak_paths = sorted(
        path
        for path in input_dir.rglob(f"*{PEAKS_SUFFIX}")
        if path.is_file() and path.name.endswith(PEAKS_SUFFIX)
    )
    if not peak_paths:
        raise FileNotFoundError(
            f"No '*{PEAKS_SUFFIX}' files were found under {input_dir}"
        )
    return peak_paths


def validate_output_preflight(
    peak_paths: Sequence[Path],
    summary_path: Path,
    *,
    overwrite: bool,
) -> None:
    """Fail before processing if outputs already exist without permission."""
    if overwrite:
        return

    existing: list[Path] = []
    for peak_path in peak_paths:
        existing.extend(
            path for path in build_output_paths(peak_path) if path.exists()
        )
    if summary_path.exists():
        existing.append(summary_path)

    if existing:
        preview = "\n".join(f"  - {path}" for path in existing[:10])
        remaining = len(existing) - min(len(existing), 10)
        if remaining:
            preview += f"\n  ... and {remaining} more"
        raise FileExistsError(
            "One or more generated outputs already exist. Use --overwrite "
            f"to replace them:\n{preview}"
        )


def run(
    input_dir: Path,
    eqtl_path: Path,
    summary_path: Path,
    *,
    overwrite: bool,
) -> pd.DataFrame:
    """Run the complete exact-eQTL peak-overlap workflow."""
    input_dir = input_dir.resolve()
    eqtl_path = eqtl_path.resolve()
    summary_path = summary_path.resolve()

    peak_paths = discover_peak_files(input_dir)
    validate_output_preflight(
        peak_paths,
        summary_path,
        overwrite=overwrite,
    )

    eqtl_variants, eqtl_statistics = load_eqtl_variants(eqtl_path)

    summaries = [
        process_peak_file(
            peak_path,
            input_dir,
            eqtl_variants,
            overwrite=overwrite,
        )
        for peak_path in peak_paths
    ]
    summary = pd.DataFrame(summaries)
    summary.insert(
        0,
        "eqtl_source_path",
        str(eqtl_path),
    )
    summary.insert(
        1,
        "n_eqtl_variant_gene_pairs_in_source",
        eqtl_statistics["n_variant_gene_pairs"],
    )
    summary.insert(
        2,
        "n_unique_eqtl_variants_in_source",
        eqtl_statistics["n_unique_variants"],
    )
    summary.insert(
        3,
        "n_unique_eqtl_genes_in_source",
        eqtl_statistics["n_unique_genes"],
    )

    write_dataframe_atomic(
        summary,
        summary_path,
        header=True,
        index=False,
    )
    LOGGER.info(
        "Wrote dataset-wide summary for %d peak files to %s",
        len(summary),
        summary_path,
    )
    return summary


def main() -> int:
    """CLI entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    args = parse_args()
    LOGGER.info("create_peaks_without_eqt.py version %s", SCRIPT_VERSION)

    try:
        run(
            input_dir=args.input_dir,
            eqtl_path=args.eqtl_path,
            summary_path=args.summary_path,
            overwrite=args.overwrite,
        )
    except Exception:
        LOGGER.exception("eQTL peak-overlap processing failed")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
