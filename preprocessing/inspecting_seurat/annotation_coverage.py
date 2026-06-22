"""Per-sample cell-type annotation coverage analysis.

Reads the exported Seurat cell metadata and, for each sample, reports what
fraction of cells carry a meaningful cell-type label versus an unannotated
placeholder value (NA / "Unknow" / "Unknown" / empty string).

Outputs
-------
- Printed summary table to stdout.
- ``annotation_coverage_by_sample.csv`` written to ANNOTATION_DIR.
"""

import sys
import pandas as pd
from constants import ANNOTATIONS_PATH, ANNOTATION_DIR

# Values that represent the absence of a real cell-type label.
UNANNOTATED_VALUES: frozenset[str] = frozenset({"Unknow", "Unknown", ""})

# Columns to evaluate; the first is the primary annotation, the second is
# used as a cross-check (it mirrors the active Seurat identity).
ANNOTATION_COLUMNS: list[str] = ["CellType.2", "active_ident"]


def is_annotated(series: pd.Series) -> pd.Series:
    """Return a boolean Series: True where the cell has a real cell-type label."""
    return series.notna() & ~series.isin(UNANNOTATED_VALUES)


def coverage_by_sample(
    metadata: pd.DataFrame,
    annotation_col: str,
) -> pd.DataFrame:
    """Compute per-sample annotation coverage for one annotation column."""
    annotated = is_annotated(metadata[annotation_col])

    grouped = (
        metadata.assign(_annotated=annotated)
        .groupby("sample", sort=True)["_annotated"]
        .agg(n_total="count", n_annotated="sum")
        .reset_index()
    )

    grouped["n_unannotated"] = grouped["n_total"] - grouped["n_annotated"]
    grouped["pct_annotated"] = (
        grouped["n_annotated"] / grouped["n_total"] * 100
    ).round(1)
    grouped.insert(0, "annotation_column", annotation_col)

    return grouped


def overall_coverage(df: pd.DataFrame) -> str:
    """Return a single-line overall coverage summary string."""
    total = df["n_total"].sum()
    annotated = df["n_annotated"].sum()
    pct = annotated / total * 100 if total > 0 else 0.0
    return (
        f"Overall: {annotated:,} / {total:,} cells annotated "
        f"({pct:.1f}%)"
    )


def main() -> None:
    if not ANNOTATIONS_PATH.is_file():
        print(
            f"ERROR: Cell metadata file not found:\n  {ANNOTATIONS_PATH}\n"
            "Run preprocessing/inspecting_seurat/inspect.py first to export it.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Reading cell metadata from:\n  {ANNOTATIONS_PATH}\n")
    metadata = pd.read_csv(ANNOTATIONS_PATH, low_memory=False)
    print(f"Loaded {len(metadata):,} cells.\n")

    all_frames: list[pd.DataFrame] = []

    for col in ANNOTATION_COLUMNS:
        if col not in metadata.columns:
            print(f"WARNING: Column '{col}' not found in metadata — skipping.")
            continue

        df = coverage_by_sample(metadata, col)
        all_frames.append(df)

        print(f"=== Annotation coverage — {col} ===")
        print(
            df.to_string(
                index=False,
                columns=[
                    "sample",
                    "n_total",
                    "n_annotated",
                    "n_unannotated",
                    "pct_annotated",
                ],
            )
        )
        print()
        print(overall_coverage(df))

        # Flag samples that look under-annotated (< 50 %).
        low_coverage = df[df["pct_annotated"] < 50.0]
        if not low_coverage.empty:
            print(
                f"\nSamples with < 50 % annotated cells ({col}):"
            )
            for _, row in low_coverage.iterrows():
                print(
                    f"  {row['sample']}: "
                    f"{row['n_annotated']:,} / {row['n_total']:,} "
                    f"({row['pct_annotated']}%)"
                )
        print()

    if not all_frames:
        print("No annotation columns found in metadata.", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(all_frames, ignore_index=True)

    ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)
    output_path = ANNOTATION_DIR / "annotation_coverage_by_sample.csv"
    combined.to_csv(output_path, index=False)
    print(f"Coverage table written to:\n  {output_path}")


if __name__ == "__main__":
    main()
