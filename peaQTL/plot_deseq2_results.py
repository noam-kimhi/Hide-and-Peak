#!/usr/bin/env python3
"""Create diagnostic and presentation plots for DESeq2 peak results."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, PercentFormatter

import constants as C


REQUIRED_RESULT_COLUMNS = {
    "peak_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}

NUMERIC_RESULT_COLUMNS = (
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)

REQUIRED_SUMMARY_COLUMNS = {
    "cell_type",
    "n_peaks_total",
    "n_padj_available",
    "n_pvalue_lt_0p05",
    "n_padj_lt_0p1",
    "n_padj_lt_0p05",
    "n_padj_lt_0p01",
    "default_n_significant",
    "default_n_significant_positive_log2fc",
    "default_n_significant_negative_log2fc",
}

REQUIRED_GRID_COLUMNS = {
    "cell_type",
    "padj_threshold",
    "abs_log2fc_threshold",
    "n_significant",
}

CHI_SQUARE_1DF_MEDIAN = 0.4549364231195727

BACKGROUND_COLOR = "#bdbdbd"
SMALL_EFFECT_COLOR = "#f4a261"
POSITIVE_COLOR = "#d95f02"
NEGATIVE_COLOR = "#1f78b4"
PADJ_AVAILABLE_COLOR = "#4c78a8"
PADJ_MISSING_COLOR = "#d9d9d9"


def configure_plotting() -> None:
    """Set presentation-friendly Matplotlib defaults."""
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.frameon": False,
            "font.size": 9,
            "savefig.facecolor": "white",
        }
    )


def validate_columns(
    dataframe: pd.DataFrame,
    required: set[str],
    path: Path,
) -> None:
    """Validate that a table contains all required columns."""
    missing = required.difference(dataframe.columns)

    if missing:
        raise ValueError(
            f"{path} is missing required columns: "
            f"{', '.join(sorted(missing))}"
        )


def slugify(value: str) -> str:
    """Return a filesystem-safe name."""
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return slug.strip("_") or "cell_type"


def cell_type_from_path(path: Path) -> str:
    """Extract the cell type from a DESeq2 result filename."""
    if not path.name.endswith(C.DESEQ2_RESULTS_SUFFIX):
        raise ValueError(
            f"Unexpected filename {path.name!r}; expected suffix "
            f"{C.DESEQ2_RESULTS_SUFFIX!r}."
        )

    return path.name.removesuffix(C.DESEQ2_RESULTS_SUFFIX)


def discover_result_files() -> list[Path]:
    """Find all cell-type DESeq2 result files."""
    if not C.DESEQ2_RESULTS_DIR.exists():
        raise FileNotFoundError(
            f"Result directory does not exist: "
            f"{C.DESEQ2_RESULTS_DIR}"
        )

    files = sorted(
        path
        for path in C.DESEQ2_RESULTS_DIR.glob(
            f"*{C.DESEQ2_RESULTS_SUFFIX}"
        )
        if path.is_file()
    )

    if not files:
        raise FileNotFoundError(
            f"No files ending with {C.DESEQ2_RESULTS_SUFFIX!r} were "
            f"found in {C.DESEQ2_RESULTS_DIR}."
        )

    return files


def load_result_file(path: Path) -> pd.DataFrame:
    """Load and clean one DESeq2 result file."""
    dataframe = pd.read_csv(path)
    validate_columns(dataframe, REQUIRED_RESULT_COLUMNS, path)

    dataframe = dataframe.copy()
    dataframe["peak_id"] = dataframe["peak_id"].astype("string")

    for column in NUMERIC_RESULT_COLUMNS:
        dataframe[column] = pd.to_numeric(
            dataframe[column],
            errors="coerce",
        )

    for column in ("pvalue", "padj"):
        valid = (
            np.isfinite(dataframe[column])
            & dataframe[column].between(0.0, 1.0)
        )
        dataframe.loc[~valid, column] = np.nan

    return dataframe


def load_summary() -> pd.DataFrame:
    """Load the previously generated cross-cell-type summary."""
    if not C.DESEQ2_SUMMARY_CSV.exists():
        raise FileNotFoundError(
            f"Missing summary file: {C.DESEQ2_SUMMARY_CSV}. "
            "Run the analysis script first."
        )

    summary = pd.read_csv(C.DESEQ2_SUMMARY_CSV)
    validate_columns(
        summary,
        REQUIRED_SUMMARY_COLUMNS,
        C.DESEQ2_SUMMARY_CSV,
    )

    return summary


def load_threshold_grid() -> pd.DataFrame:
    """Load the previously generated threshold-sensitivity table."""
    if not C.DESEQ2_THRESHOLD_GRID_CSV.exists():
        raise FileNotFoundError(
            f"Missing threshold-grid file: "
            f"{C.DESEQ2_THRESHOLD_GRID_CSV}. "
            "Run the analysis script first."
        )

    grid = pd.read_csv(C.DESEQ2_THRESHOLD_GRID_CSV)
    validate_columns(
        grid,
        REQUIRED_GRID_COLUMNS,
        C.DESEQ2_THRESHOLD_GRID_CSV,
    )

    return grid


def ordered_cell_types(summary: pd.DataFrame) -> list[str]:
    """Order cell types by their default number of discoveries."""
    ordered = summary.sort_values(
        ["default_n_significant", "cell_type"],
        ascending=[False, True],
    )

    return ordered["cell_type"].astype(str).tolist()


def save_figure(
    figure: Figure,
    path: Path,
) -> None:
    """Save and close one figure."""
    path.parent.mkdir(parents=True, exist_ok=True)

    figure.savefig(
        path,
        dpi=C.PLOT_DPI,
        bbox_inches="tight",
    )
    plt.close(figure)


def negative_log10(
    values: pd.Series | np.ndarray,
) -> np.ndarray:
    """Calculate -log10 safely, including exact zero p-values."""
    array = np.asarray(values, dtype=float)
    result = np.full(array.shape, np.nan, dtype=float)

    valid = (
        np.isfinite(array)
        & (array >= 0.0)
        & (array <= 1.0)
    )

    result[valid] = -np.log10(
        np.clip(
            array[valid],
            np.finfo(float).tiny,
            1.0,
        )
    )

    return result


def clip_extreme_values(
    values: np.ndarray,
    maximum_cap: float,
    minimum_cap: float = 2.0,
) -> tuple[np.ndarray, float, int]:
    """Clip only the extreme upper tail for readable plotting."""
    finite = values[np.isfinite(values)]

    if finite.size == 0:
        return values, minimum_cap, 0

    observed_maximum = float(np.max(finite))
    quantile_cap = float(np.quantile(finite, 0.995)) * 1.15

    cap = min(
        max(minimum_cap, quantile_cap),
        maximum_cap,
    )
    cap = min(cap, observed_maximum)

    clipped = np.minimum(values, cap)
    n_clipped = int(np.sum(values > cap))

    return clipped, cap, n_clipped


def fdr_mask(
    dataframe: pd.DataFrame,
) -> pd.Series:
    """Return peaks passing the configured adjusted-p-value threshold."""
    return (
        dataframe["padj"].notna()
        & (
            dataframe["padj"]
            < C.DESEQ2_DEFAULT_PADJ_THRESHOLD
        )
    )


def significant_mask(
    dataframe: pd.DataFrame,
) -> pd.Series:
    """Return peaks passing both FDR and effect-size thresholds."""
    return (
        fdr_mask(dataframe)
        & dataframe["log2FoldChange"].notna()
        & (
            dataframe["log2FoldChange"].abs()
            >= C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
        )
    )


def classify_peaks(
    dataframe: pd.DataFrame,
) -> pd.Series:
    """Assign each peak to a plotting category."""
    categories = pd.Series(
        "Not FDR-significant",
        index=dataframe.index,
        dtype="object",
    )

    passes_fdr = fdr_mask(dataframe)

    passes_effect = (
        dataframe["log2FoldChange"].abs()
        >= C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
    )

    categories.loc[
        passes_fdr & ~passes_effect
    ] = "FDR-significant, below effect threshold"

    categories.loc[
        passes_fdr
        & passes_effect
        & (dataframe["log2FoldChange"] > 0)
    ] = C.DESEQ2_POSITIVE_LFC_LABEL

    categories.loc[
        passes_fdr
        & passes_effect
        & (dataframe["log2FoldChange"] < 0)
    ] = C.DESEQ2_NEGATIVE_LFC_LABEL

    return categories


def scatter_categories(
    axis: Axes,
    dataframe: pd.DataFrame,
    x_column: str,
    y_column: str,
) -> None:
    """Plot categories with discoveries drawn above background points."""
    styles = (
        (
            "Not FDR-significant",
            BACKGROUND_COLOR,
            7,
            0.25,
        ),
        (
            "FDR-significant, below effect threshold",
            SMALL_EFFECT_COLOR,
            13,
            0.80,
        ),
        (
            C.DESEQ2_NEGATIVE_LFC_LABEL,
            NEGATIVE_COLOR,
            18,
            0.90,
        ),
        (
            C.DESEQ2_POSITIVE_LFC_LABEL,
            POSITIVE_COLOR,
            18,
            0.90,
        ),
    )

    for category, color, size, alpha in styles:
        subset = dataframe[
            dataframe["category"] == category
        ]

        if subset.empty:
            continue

        axis.scatter(
            subset[x_column],
            subset[y_column],
            s=size,
            alpha=alpha,
            color=color,
            edgecolors="none",
            rasterized=True,
            label=f"{category} ({len(subset):,})",
        )


def peaks_to_label(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """Choose significant peaks to label in a volcano plot."""
    significant = dataframe[
        significant_mask(dataframe)
    ].copy()

    significant = significant.sort_values(
        ["padj", "pvalue"],
        ascending=[True, True],
        na_position="last",
    )

    if len(significant) <= 15:
        return significant

    return significant.head(
        C.DESEQ2_TOP_PEAK_LABELS
    )


def plot_volcano(
    axis: Axes,
    dataframe: pd.DataFrame,
) -> None:
    """Draw a volcano plot using raw p-values and FDR-based colors."""
    plot_data = dataframe[
        dataframe["pvalue"].notna()
        & dataframe["log2FoldChange"].notna()
    ].copy()

    raw_y = negative_log10(
        plot_data["pvalue"]
    )

    (
        plot_data["negative_log10_pvalue"],
        y_cap,
        n_clipped,
    ) = clip_extreme_values(
        raw_y,
        C.DESEQ2_MAX_NEG_LOG10_FOR_PLOTS,
    )

    plot_data["category"] = classify_peaks(
        plot_data
    )

    scatter_categories(
        axis,
        plot_data,
        "log2FoldChange",
        "negative_log10_pvalue",
    )

    effect_threshold = (
        C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
    )

    axis.axvline(
        -effect_threshold,
        color="black",
        linestyle="--",
        alpha=0.55,
    )

    axis.axvline(
        effect_threshold,
        color="black",
        linestyle="--",
        alpha=0.55,
    )

    nominal_y = -math.log10(0.05)

    if nominal_y <= y_cap:
        axis.axhline(
            nominal_y,
            color="black",
            linestyle=":",
            alpha=0.45,
        )

        axis.text(
            0.99,
            nominal_y,
            " nominal p = 0.05",
            transform=axis.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=7,
        )

    for _, row in peaks_to_label(
        plot_data
    ).iterrows():
        axis.annotate(
            str(row["peak_id"]),
            (
                float(row["log2FoldChange"]),
                float(row["negative_log10_pvalue"]),
            ),
            xytext=(3, 4),
            textcoords="offset points",
            fontsize=6,
        )

    axis.set_xlabel("log2 fold change")
    axis.set_ylabel("-log10(raw p-value)")
    axis.set_title("Volcano plot")
    axis.grid(alpha=0.15)

    axis.legend(
        loc="best",
        fontsize=7,
        markerscale=1.2,
    )

    if n_clipped:
        axis.text(
            0.02,
            0.98,
            (
                f"{n_clipped} extreme point(s) "
                f"clipped at {y_cap:.1f}"
            ),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=7,
        )


def plot_ma(
    axis: Axes,
    dataframe: pd.DataFrame,
) -> None:
    """Draw an MA plot of abundance against fold change."""
    plot_data = dataframe[
        dataframe["baseMean"].notna()
        & dataframe["log2FoldChange"].notna()
        & (dataframe["baseMean"] >= 0)
    ].copy()

    plot_data["log10_base_mean_plus_one"] = (
        np.log10(plot_data["baseMean"] + 1.0)
    )

    plot_data["category"] = classify_peaks(
        plot_data
    )

    finite_effects = (
        plot_data["log2FoldChange"].abs()
    )
    finite_effects = finite_effects[
        np.isfinite(finite_effects)
    ]

    if finite_effects.empty:
        effect_limit = 1.0
    else:
        effect_limit = min(
            C.DESEQ2_MA_MAX_ABS_LOG2FC,
            max(
                1.0,
                float(
                    finite_effects.quantile(0.995)
                )
                * 1.10,
            ),
        )

    plot_data["display_log2fc"] = (
        plot_data["log2FoldChange"].clip(
            -effect_limit,
            effect_limit,
        )
    )

    scatter_categories(
        axis,
        plot_data,
        "log10_base_mean_plus_one",
        "display_log2fc",
    )

    effect_threshold = (
        C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD
    )

    axis.axhline(
        0,
        color="black",
        linewidth=1,
        alpha=0.65,
    )

    axis.axhline(
        effect_threshold,
        color="black",
        linestyle="--",
        alpha=0.45,
    )

    axis.axhline(
        -effect_threshold,
        color="black",
        linestyle="--",
        alpha=0.45,
    )

    n_clipped = int(
        (
            plot_data["log2FoldChange"].abs()
            > effect_limit
        ).sum()
    )

    axis.set_ylim(
        -effect_limit,
        effect_limit,
    )

    axis.set_xlabel("log10(baseMean + 1)")
    axis.set_ylabel("log2 fold change")
    axis.set_title("MA plot")
    axis.grid(alpha=0.15)

    if n_clipped:
        axis.text(
            0.02,
            0.98,
            (
                f"{n_clipped} extreme "
                f"effect(s) clipped"
            ),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=7,
        )


def plot_qq(
    axis: Axes,
    dataframe: pd.DataFrame,
) -> None:
    """Draw a raw-p-value QQ plot."""
    pvalues = (
        dataframe["pvalue"]
        .dropna()
        .to_numpy(dtype=float)
    )

    pvalues = pvalues[
        np.isfinite(pvalues)
        & (pvalues >= 0.0)
        & (pvalues <= 1.0)
    ]

    if pvalues.size == 0:
        axis.text(
            0.5,
            0.5,
            "No valid p-values",
            transform=axis.transAxes,
            ha="center",
            va="center",
        )
        axis.set_title("Raw-p-value QQ plot")
        return

    observed_pvalues = np.sort(
        np.clip(
            pvalues,
            np.finfo(float).tiny,
            1.0,
        )
    )

    expected_pvalues = (
        np.arange(
            1,
            observed_pvalues.size + 1,
        )
        - 0.5
    ) / observed_pvalues.size

    observed = -np.log10(
        observed_pvalues
    )

    expected = -np.log10(
        expected_pvalues
    )

    _, cap, n_clipped = clip_extreme_values(
        np.concatenate(
            (observed, expected)
        ),
        C.DESEQ2_MAX_NEG_LOG10_FOR_PLOTS,
    )

    axis.scatter(
        np.minimum(expected, cap),
        np.minimum(observed, cap),
        s=7,
        alpha=0.45,
        color="#4d4d4d",
        edgecolors="none",
        rasterized=True,
    )

    axis.plot(
        [0, cap],
        [0, cap],
        color="black",
        linestyle="--",
        linewidth=1,
    )

    statistic = pd.to_numeric(
        dataframe["stat"],
        errors="coerce",
    )

    statistic = statistic[
        np.isfinite(statistic)
    ]

    annotation_lines = [
        f"n = {pvalues.size:,}"
    ]

    if not statistic.empty:
        lambda_gc = (
            np.median(
                np.square(statistic)
            )
            / CHI_SQUARE_1DF_MEDIAN
        )

        annotation_lines.append(
            f"λGC = {lambda_gc:.3f}"
        )

    if n_clipped:
        annotation_lines.append(
            f"{n_clipped} value(s) clipped"
        )

    axis.text(
        0.04,
        0.96,
        "\n".join(annotation_lines),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )

    axis.set_xlim(
        0,
        cap * 1.03,
    )

    axis.set_ylim(
        0,
        cap * 1.03,
    )

    axis.set_xlabel(
        "Expected -log10(p-value)"
    )

    axis.set_ylabel(
        "Observed -log10(p-value)"
    )

    axis.set_title(
        "Raw-p-value QQ plot"
    )

    axis.grid(alpha=0.15)


def create_dashboard(
    path: Path,
    dataframe: pd.DataFrame,
) -> None:
    """Create a volcano, MA, and QQ dashboard for one cell type."""
    cell_type = cell_type_from_path(path)

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(18, 5.5),
        constrained_layout=True,
    )

    plot_volcano(
        axes[0],
        dataframe,
    )

    plot_ma(
        axes[1],
        dataframe,
    )

    plot_qq(
        axes[2],
        dataframe,
    )

    n_significant = int(
        significant_mask(dataframe).sum()
    )

    figure.suptitle(
        (
            f"{cell_type}: differential "
            f"chromatin accessibility\n"
            f"Highlighted: padj < "
            f"{C.DESEQ2_DEFAULT_PADJ_THRESHOLD:g}, "
            f"|log2FC| ≥ "
            f"{C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD:g}; "
            f"n = {n_significant:,}"
        ),
        fontsize=14,
        fontweight="bold",
    )

    output_path = (
        C.DESEQ2_PER_CELL_PLOTS_DIR
        / (
            f"{slugify(cell_type)}"
            f"_differential_accessibility_dashboard."
            f"{C.DESEQ2_PLOT_FORMAT}"
        )
    )

    save_figure(
        figure,
        output_path,
    )


def create_forest_plot(
    path: Path,
    dataframe: pd.DataFrame,
) -> bool:
    """Plot top significant peaks with 95% Wald intervals."""
    cell_type = cell_type_from_path(path)

    candidates = dataframe[
        significant_mask(dataframe)
    ].copy()

    candidates = candidates[
        candidates["lfcSE"].notna()
        & np.isfinite(candidates["lfcSE"])
        & (candidates["lfcSE"] >= 0)
    ]

    if candidates.empty:
        print(
            f"  Skipping forest plot for "
            f"{cell_type}: no significant "
            f"peak has a valid lfcSE."
        )
        return False

    total_significant = int(
        significant_mask(dataframe).sum()
    )

    candidates["abs_log2fc"] = (
        candidates["log2FoldChange"].abs()
    )

    candidates = (
        candidates.sort_values(
            [
                "padj",
                "pvalue",
                "abs_log2fc",
            ],
            ascending=[
                True,
                True,
                False,
            ],
            na_position="last",
        )
        .head(C.DESEQ2_TOP_FOREST_PEAKS)
        .copy()
    )

    candidates["ci95"] = (
        1.96 * candidates["lfcSE"]
    )

    candidates = (
        candidates.sort_values(
            "log2FoldChange"
        )
        .reset_index(drop=True)
    )

    figure_height = max(
        4.5,
        0.34 * len(candidates) + 1.8,
    )

    figure, axis = plt.subplots(
        figsize=(10, figure_height),
        constrained_layout=True,
    )

    for index, row in candidates.iterrows():
        color = (
            POSITIVE_COLOR
            if row["log2FoldChange"] >= 0
            else NEGATIVE_COLOR
        )

        axis.errorbar(
            float(row["log2FoldChange"]),
            index,
            xerr=float(row["ci95"]),
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=1.3,
            capsize=2.5,
            markersize=5,
        )

    axis.axvline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )

    axis.set_yticks(
        np.arange(len(candidates))
    )

    axis.set_yticklabels(
        candidates["peak_id"].astype(str)
    )

    axis.set_xlabel(
        "log2 fold change with 95% Wald interval"
    )

    axis.set_ylabel("Peak")

    axis.set_title(
        (
            f"{cell_type}: top significant "
            f"differential peaks\n"
            f"Showing {len(candidates):,} of "
            f"{total_significant:,} significant peaks"
        )
    )

    axis.grid(
        axis="x",
        alpha=0.2,
    )

    axis.text(
        0.01,
        1.01,
        f"← {C.DESEQ2_NEGATIVE_LFC_LABEL}",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        color=NEGATIVE_COLOR,
        fontsize=8,
    )

    axis.text(
        0.99,
        1.01,
        f"{C.DESEQ2_POSITIVE_LFC_LABEL} →",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        color=POSITIVE_COLOR,
        fontsize=8,
    )

    output_path = (
        C.DESEQ2_PER_CELL_PLOTS_DIR
        / (
            f"{slugify(cell_type)}"
            f"_top_significant_peaks_forest."
            f"{C.DESEQ2_PLOT_FORMAT}"
        )
    )

    save_figure(
        figure,
        output_path,
    )

    return True


def create_pvalue_histograms(
    result_files: Iterable[Path],
) -> None:
    """Plot raw-p-value distributions for all cell types."""
    files = list(result_files)

    n_columns = 3
    n_rows = math.ceil(
        len(files) / n_columns
    )

    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(
            5.2 * n_columns,
            3.8 * n_rows,
        ),
        squeeze=False,
        constrained_layout=True,
    )

    bins = np.linspace(
        0.0,
        1.0,
        C.DESEQ2_PVALUE_HISTOGRAM_BINS + 1,
    )

    for axis, path in zip(
        axes.flat,
        files,
    ):
        dataframe = load_result_file(path)
        cell_type = cell_type_from_path(path)

        pvalues = (
            dataframe["pvalue"]
            .dropna()
        )

        pvalues = pvalues[
            np.isfinite(pvalues)
            & pvalues.between(0.0, 1.0)
        ]

        axis.hist(
            pvalues,
            bins=bins,
            density=True,
            color=C.PLOTS_MAIN_COLOR,
            alpha=0.8,
            edgecolor="white",
            linewidth=0.4,
        )

        axis.axhline(
            1.0,
            color="black",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )

        axis.set_xlim(
            0.0,
            1.0,
        )

        axis.set_xlabel(
            "Raw p-value"
        )

        axis.set_ylabel(
            "Density"
        )

        axis.set_title(
            f"{cell_type}\nn = {len(pvalues):,}"
        )

        axis.grid(alpha=0.12)

    for axis in axes.flat[len(files):]:
        axis.set_visible(False)

    figure.suptitle(
        (
            "Raw p-value distributions by cell type\n"
            "The dashed line is the "
            "Uniform(0,1) null density"
        ),
        fontsize=14,
        fontweight="bold",
    )

    output_path = (
        C.DESEQ2_PLOTS_DIR
        / (
            "all_cell_types_pvalue_histograms."
            f"{C.DESEQ2_PLOT_FORMAT}"
        )
    )

    save_figure(
        figure,
        output_path,
    )


def create_nominal_vs_fdr_plot(
    summary: pd.DataFrame,
) -> None:
    """Compare nominal and FDR-adjusted peak counts."""
    cell_order = ordered_cell_types(
        summary
    )

    ordered = (
        summary.set_index("cell_type")
        .loc[cell_order]
        .reset_index()
    )

    criteria = (
        (
            "Raw p < 0.05",
            "n_pvalue_lt_0p05",
        ),
        (
            "padj < 0.10",
            "n_padj_lt_0p1",
        ),
        (
            "padj < 0.05",
            "n_padj_lt_0p05",
        ),
        (
            "padj < 0.01",
            "n_padj_lt_0p01",
        ),
    )

    x_positions = np.arange(
        len(ordered)
    )

    bar_width = 0.19

    figure, axis = plt.subplots(
        figsize=(
            max(
                10,
                1.4 * len(ordered),
            ),
            6,
        ),
        constrained_layout=True,
    )

    for criterion_index, (
        label,
        column,
    ) in enumerate(criteria):
        counts = ordered[
            column
        ].to_numpy(dtype=float)

        positions = (
            x_positions
            + (
                criterion_index
                - (len(criteria) - 1) / 2
            )
            * bar_width
        )

        bars = axis.bar(
            positions,
            np.where(
                counts > 0,
                counts,
                np.nan,
            ),
            width=bar_width,
            label=label,
        )

        for position, count, _ in zip(
            positions,
            counts,
            bars,
        ):
            if count > 0:
                axis.text(
                    position,
                    count * 1.12,
                    f"{int(count):,}",
                    ha="center",
                    va="bottom",
                    rotation=90,
                    fontsize=6,
                )
            else:
                axis.text(
                    position,
                    0.72,
                    "0",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                )

    axis.set_yscale("log")

    axis.set_ylim(
        0.6,
        max(
            10.0,
            float(
                ordered[
                    "n_pvalue_lt_0p05"
                ].max()
            )
            * 2.5,
        ),
    )

    axis.set_xticks(
        x_positions
    )

    axis.set_xticklabels(
        ordered["cell_type"],
        rotation=30,
        ha="right",
    )

    axis.set_ylabel(
        "Number of peaks (log scale)"
    )

    axis.set_title(
        (
            "Nominal evidence versus "
            "multiple-testing-adjusted discoveries"
        )
    )

    axis.legend(ncols=2)

    axis.grid(
        axis="y",
        which="both",
        alpha=0.2,
    )

    output_path = (
        C.DESEQ2_PLOTS_DIR
        / (
            "all_cell_types_nominal_vs_fdr_counts."
            f"{C.DESEQ2_PLOT_FORMAT}"
        )
    )

    save_figure(
        figure,
        output_path,
    )


def create_direction_plot(
    summary: pd.DataFrame,
) -> None:
    """Create a diverging count plot by fold-change direction."""
    cell_order = ordered_cell_types(
        summary
    )

    ordered = (
        summary.set_index("cell_type")
        .loc[cell_order]
        .reset_index()
        .iloc[::-1]
        .reset_index(drop=True)
    )

    positive = ordered[
        "default_n_significant_positive_log2fc"
    ].to_numpy(dtype=int)

    negative = ordered[
        "default_n_significant_negative_log2fc"
    ].to_numpy(dtype=int)

    y_positions = np.arange(
        len(ordered)
    )

    figure, axis = plt.subplots(
        figsize=(
            10,
            max(
                4.8,
                0.65 * len(ordered),
            ),
        ),
        constrained_layout=True,
    )

    axis.barh(
        y_positions,
        -negative,
        color=NEGATIVE_COLOR,
        label=C.DESEQ2_NEGATIVE_LFC_LABEL,
    )

    axis.barh(
        y_positions,
        positive,
        color=POSITIVE_COLOR,
        label=C.DESEQ2_POSITIVE_LFC_LABEL,
    )

    for (
        y_position,
        positive_count,
        negative_count,
    ) in zip(
        y_positions,
        positive,
        negative,
    ):
        if negative_count:
            axis.text(
                -negative_count - 0.8,
                y_position,
                str(negative_count),
                ha="right",
                va="center",
            )

        if positive_count:
            axis.text(
                positive_count + 0.8,
                y_position,
                str(positive_count),
                ha="left",
                va="center",
            )

        if (
            positive_count == 0
            and negative_count == 0
        ):
            axis.text(
                0,
                y_position,
                "0",
                ha="center",
                va="center",
            )

    maximum_positive = (
        int(np.max(positive))
        if positive.size
        else 0
    )

    maximum_negative = (
        int(np.max(negative))
        if negative.size
        else 0
    )

    maximum_count = max(
        1,
        maximum_positive,
        maximum_negative,
    )

    axis.set_xlim(
        -maximum_count * 1.35,
        maximum_count * 1.35,
    )

    axis.axvline(
        0,
        color="black",
        linewidth=1,
    )

    axis.set_yticks(
        y_positions
    )

    axis.set_yticklabels(
        ordered["cell_type"]
    )

    axis.xaxis.set_major_formatter(
        FuncFormatter(
            lambda value, _: str(
                abs(int(value))
            )
        )
    )

    axis.set_xlabel(
        "Number of significant peaks"
    )

    axis.set_title(
        (
            "Significant differential peaks "
            "by accessibility direction\n"
            f"padj < "
            f"{C.DESEQ2_DEFAULT_PADJ_THRESHOLD:g}, "
            f"|log2FC| ≥ "
            f"{C.DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD:g}"
        )
    )

    axis.legend(
        loc="lower right"
    )

    axis.grid(
        axis="x",
        alpha=0.2,
    )

    output_path = (
        C.DESEQ2_PLOTS_DIR
        / (
            "all_cell_types_"
            "significant_direction_counts."
            f"{C.DESEQ2_PLOT_FORMAT}"
        )
    )

    save_figure(
        figure,
        output_path,
    )


def create_independent_filtering_plot(
    summary: pd.DataFrame,
) -> None:
    """Plot fractions with available and missing adjusted p-values."""
    cell_order = ordered_cell_types(
        summary
    )

    ordered = (
        summary.set_index("cell_type")
        .loc[cell_order]
        .reset_index()
    )

    total = ordered[
        "n_peaks_total"
    ].to_numpy(dtype=float)

    available = ordered[
        "n_padj_available"
    ].to_numpy(dtype=float)

    missing = total - available

    available_fraction = np.divide(
        available,
        total,
        out=np.zeros_like(available),
        where=total > 0,
    )

    missing_fraction = np.divide(
        missing,
        total,
        out=np.zeros_like(missing),
        where=total > 0,
    )

    x_positions = np.arange(
        len(ordered)
    )

    figure, axis = plt.subplots(
        figsize=(
            max(
                10,
                1.35 * len(ordered),
            ),
            5.8,
        ),
        constrained_layout=True,
    )

    axis.bar(
        x_positions,
        available_fraction,
        color=PADJ_AVAILABLE_COLOR,
        label="padj available",
    )

    axis.bar(
        x_positions,
        missing_fraction,
        bottom=available_fraction,
        color=PADJ_MISSING_COLOR,
        label=(
            "padj missing after "
            "independent filtering"
        ),
    )

    for index, fraction in enumerate(
        missing_fraction
    ):
        if fraction >= 0.03:
            text_y = (
                available_fraction[index]
                + fraction / 2
            )
            vertical_alignment = "center"
        else:
            text_y = min(
                1.03,
                available_fraction[index] + 0.015,
            )
            vertical_alignment = "bottom"

        axis.text(
            index,
            text_y,
            f"{100 * fraction:.1f}%",
            ha="center",
            va=vertical_alignment,
            fontsize=8,
        )

    axis.set_xticks(
        x_positions
    )

    axis.set_xticklabels(
        ordered["cell_type"],
        rotation=30,
        ha="right",
    )

    axis.set_ylim(
        0,
        1.08,
    )

    axis.yaxis.set_major_formatter(
        PercentFormatter(1.0)
    )

    axis.set_ylabel(
        "Fraction of input peaks"
    )

    axis.set_title(
        (
            "Availability of adjusted p-values "
            "after DESeq2 independent filtering"
        )
    )

    axis.legend(
        loc="upper right"
    )

    axis.grid(
        axis="y",
        alpha=0.2,
    )

    output_path = (
        C.DESEQ2_PLOTS_DIR
        / (
            "all_cell_types_independent_filtering."
            f"{C.DESEQ2_PLOT_FORMAT}"
        )
    )

    save_figure(
        figure,
        output_path,
    )


def create_threshold_sensitivity_plot(
    grid: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    """Plot significant counts across selected FDR/effect thresholds."""
    cell_types = ordered_cell_types(
        summary
    )

    criteria = [
        (padj_threshold, effect_threshold)
        for padj_threshold
        in C.DESEQ2_PADJ_THRESHOLDS
        for effect_threshold
        in C.DESEQ2_THRESHOLD_PLOT_ABS_LOG2FC_THRESHOLDS
    ]

    counts = np.full(
        (
            len(cell_types),
            len(criteria),
        ),
        np.nan,
        dtype=float,
    )

    for row_index, cell_type in enumerate(
        cell_types
    ):
        cell_grid = grid[
            grid["cell_type"].astype(str)
            == cell_type
        ]

        for column_index, (
            padj_threshold,
            effect_threshold,
        ) in enumerate(criteria):
            match = cell_grid[
                np.isclose(
                    cell_grid["padj_threshold"],
                    padj_threshold,
                )
                & np.isclose(
                    cell_grid[
                        "abs_log2fc_threshold"
                    ],
                    effect_threshold,
                )
            ]

            if not match.empty:
                counts[
                    row_index,
                    column_index,
                ] = float(
                    match.iloc[0][
                        "n_significant"
                    ]
                )

    color_values = np.log10(
        counts + 1.0
    )

    labels = [
        (
            f"padj<{padj_threshold:g}\n"
            f"|LFC|≥{effect_threshold:g}"
        )
        for (
            padj_threshold,
            effect_threshold,
        ) in criteria
    ]

    figure, axis = plt.subplots(
        figsize=(
            max(
                11,
                1.05 * len(criteria),
            ),
            max(
                4.8,
                0.72 * len(cell_types),
            ),
        ),
        constrained_layout=True,
    )

    image = axis.imshow(
        color_values,
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
    )

    finite = color_values[
        np.isfinite(color_values)
    ]

    midpoint = (
        float(np.max(finite)) / 2
        if finite.size
        else 0.0
    )

    for row_index in range(
        counts.shape[0]
    ):
        for column_index in range(
            counts.shape[1]
        ):
            value = counts[
                row_index,
                column_index,
            ]

            if not np.isfinite(value):
                continue

            text_color = (
                "white"
                if color_values[
                    row_index,
                    column_index,
                ]
                > midpoint
                else "black"
            )

            axis.text(
                column_index,
                row_index,
                f"{int(value):,}",
                ha="center",
                va="center",
                fontsize=8,
                color=text_color,
            )

    axis.set_xticks(
        np.arange(len(labels))
    )

    axis.set_xticklabels(
        labels,
        rotation=40,
        ha="right",
    )

    axis.set_yticks(
        np.arange(len(cell_types))
    )

    axis.set_yticklabels(
        cell_types
    )

    axis.set_title(
        (
            "Sensitivity of significant-peak "
            "counts to threshold choice"
        )
    )

    colorbar = figure.colorbar(
        image,
        ax=axis,
    )

    colorbar.set_label(
        "log10(significant peak count + 1)"
    )

    output_path = (
        C.DESEQ2_PLOTS_DIR
        / (
            "all_cell_types_"
            "threshold_sensitivity_counts."
            f"{C.DESEQ2_PLOT_FORMAT}"
        )
    )

    save_figure(
        figure,
        output_path,
    )


def main() -> None:
    """Create all DESeq2 result plots."""
    configure_plotting()

    C.DESEQ2_PLOTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    C.DESEQ2_PER_CELL_PLOTS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    result_files = discover_result_files()
    summary = load_summary()
    threshold_grid = load_threshold_grid()

    preferred_order = ordered_cell_types(
        summary
    )

    path_by_cell_type = {
        cell_type_from_path(path): path
        for path in result_files
    }

    ordered_files = [
        path_by_cell_type[cell_type]
        for cell_type in preferred_order
        if cell_type in path_by_cell_type
    ]

    ordered_files.extend(
        path
        for path in result_files
        if cell_type_from_path(path)
        not in preferred_order
    )

    forest_plot_count = 0

    for path in ordered_files:
        print(
            f"Creating plots for {path.name}..."
        )

        dataframe = load_result_file(path)

        create_dashboard(
            path,
            dataframe,
        )

        if create_forest_plot(
            path,
            dataframe,
        ):
            forest_plot_count += 1

    print(
        "Creating cross-cell-type plots..."
    )

    create_pvalue_histograms(
        ordered_files
    )

    create_nominal_vs_fdr_plot(
        summary
    )

    create_direction_plot(
        summary
    )

    create_independent_filtering_plot(
        summary
    )

    create_threshold_sensitivity_plot(
        threshold_grid,
        summary,
    )

    print(
        f"Finished: {len(ordered_files)} dashboards "
        f"and {forest_plot_count} forest plots."
    )

    print(
        f"Plots saved under:\n"
        f"  {C.DESEQ2_PLOTS_DIR}"
    )


if __name__ == "__main__":
    main()