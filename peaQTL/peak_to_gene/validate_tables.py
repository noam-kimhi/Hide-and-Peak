"""Print validation summaries for generated pseudobulk count tables.

The script reads one pseudobulk CSV per standardized cell type, reports
per-sample totals and cross-sample correlations, and prints an aggregate summary
without writing new files.
"""
from __future__ import annotations

import sys
import textwrap
from pathlib import Path

# Ensure the project root (parent of this file's directory) is on sys.path
# so that constants.py can be imported regardless of working directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from constants import (
    CELL_TYPE_STANDARDIZATION,
    PSEUDOBULK_OUTPUT_DIR,
)

# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #

NORMAL_COLS = [f"Normal_rep{i}" for i in range(1, 7)]
MASH_COLS   = [f"MASH_rep{i}"   for i in range(1, 7)]


def _pearson_mean(mat: np.ndarray) -> tuple[float, float]:
    """Return mean and std of off-diagonal Pearson correlations."""
    corr = np.corrcoef(mat.T)
    n = corr.shape[0]
    upper = corr[np.triu_indices(n, k=1)]
    return float(np.mean(upper)), float(np.std(upper))


def _spearman_mean(mat: np.ndarray) -> tuple[float, float]:
    """Return mean and std of off-diagonal Spearman correlations."""
    n = mat.shape[1]
    vals: list[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            r, _ = spearmanr(mat[:, i], mat[:, j])
            vals.append(float(r))
    return float(np.mean(vals)), float(np.std(vals))


def _section(title: str) -> str:
    bar = "=" * len(title)
    return f"\n{bar}\n{title}\n{bar}"


def _subsection(title: str) -> str:
    return f"\n--- {title} ---"


# --------------------------------------------------------------------------- #
# Per-cell-type report                                                         #
# --------------------------------------------------------------------------- #

def report_cell_type(cell_type: str, df: pd.DataFrame) -> dict:
    """Print statistics for one cell type and return a summary dict."""
    n_peaks, n_samples = df.shape
    mat = df.values.astype(np.float64)

    lines: list[str] = [_section(f"Cell type: {cell_type}")]
    lines.append(f"Shape : {n_peaks:,} consensus peaks × {n_samples} samples")

    # ------------------------------------------------------------------ #
    # Per-sample statistics                                                #
    # ------------------------------------------------------------------ #
    lines.append(_subsection("Per-sample statistics"))
    header = f"  {'Sample':<18} {'Total events':>14} {'Open peaks':>11} {'Open %':>8} {'Mean (nz)':>10} {'Median (nz)':>12} {'Max':>8}"
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))

    for col in df.columns:
        v = df[col].values.astype(np.float64)
        total   = int(v.sum())
        n_open  = int((v > 0).sum())
        pct     = 100.0 * n_open / n_peaks if n_peaks else 0.0
        nz      = v[v > 0]
        mean_nz = float(nz.mean()) if len(nz) else 0.0
        med_nz  = float(np.median(nz)) if len(nz) else 0.0
        mx      = int(v.max())
        lines.append(
            f"  {col:<18} {total:>14,} {n_open:>11,} {pct:>7.1f}% "
            f"{mean_nz:>10.2f} {med_nz:>12.2f} {mx:>8,}"
        )

    # ------------------------------------------------------------------ #
    # Normal vs. MASH group comparison                                    #
    # ------------------------------------------------------------------ #
    lines.append(_subsection("Normal vs. MASH group comparison"))

    present_normal = [c for c in NORMAL_COLS if c in df.columns]
    present_mash   = [c for c in MASH_COLS   if c in df.columns]

    def _group_stats(cols: list[str]) -> dict:
        sub = df[cols].values.astype(np.float64)
        col_totals = sub.sum(axis=0)
        return {
            "n_samples"         : len(cols),
            "total_events"      : int(sub.sum()),
            "mean_events_sample": float(col_totals.mean()),
            "peaks_open_any"    : int((sub.sum(axis=1) > 0).sum()),
            "peaks_open_all"    : int((sub.min(axis=1) > 0).sum()),
        }

    g_normal = _group_stats(present_normal) if present_normal else {}
    g_mash   = _group_stats(present_mash)   if present_mash   else {}

    for label, g in (("Normal", g_normal), ("MASH", g_mash)):
        if not g:
            continue
        lines.append(
            f"  {label}: {g['total_events']:,} total events across "
            f"{g['n_samples']} samples  "
            f"(mean/sample: {g['mean_events_sample']:,.0f})  "
            f"peaks open in ≥1 sample: {g['peaks_open_any']:,}  "
            f"in all samples: {g['peaks_open_all']:,}"
        )

    # Exclusive peaks (open in ≥1 sample of one group, zero in the other)
    if present_normal and present_mash:
        normal_any = df[present_normal].sum(axis=1) > 0
        mash_any   = df[present_mash].sum(axis=1)   > 0
        excl_normal = int((normal_any & ~mash_any).sum())
        excl_mash   = int((mash_any & ~normal_any).sum())
        shared      = int((normal_any & mash_any).sum())
        lines.append(
            f"  Peaks open only in Normal: {excl_normal:,}  |  "
            f"only in MASH: {excl_mash:,}  |  shared: {shared:,}"
        )

    # ------------------------------------------------------------------ #
    # Cross-sample reproducibility                                        #
    # ------------------------------------------------------------------ #
    lines.append(_subsection("Cross-sample reproducibility"))

    # Use log1p-transformed counts for Pearson to reduce skew
    log_mat = np.log1p(mat)
    pr_mean, pr_std = _pearson_mean(log_mat)
    sp_mean, sp_std = _spearman_mean(mat)

    lines.append(f"  Pearson  r (log1p counts): {pr_mean:.4f} ± {pr_std:.4f}")
    lines.append(f"  Spearman r (raw counts)  : {sp_mean:.4f} ± {sp_std:.4f}")

    # Within-condition correlations
    for label, cols in (("Normal", present_normal), ("MASH", present_mash)):
        if len(cols) < 2:
            continue
        sub_log = np.log1p(df[cols].values.astype(np.float64))
        sub_raw = df[cols].values.astype(np.float64)
        pm, ps = _pearson_mean(sub_log)
        sm, ss = _spearman_mean(sub_raw)
        lines.append(
            f"  {label} within-group — "
            f"Pearson (log1p): {pm:.4f} ± {ps:.4f}  |  "
            f"Spearman: {sm:.4f} ± {ss:.4f}"
        )

    # ------------------------------------------------------------------ #
    # Peak coverage across all samples                                    #
    # ------------------------------------------------------------------ #
    lines.append(_subsection("Peak coverage across samples"))
    open_counts = (mat > 0).sum(axis=1)  # per peak: how many samples have it open
    for k in [12, 11, 10, 6, 3, 2, 1]:
        if k > n_samples:
            continue
        cnt = int((open_counts >= k).sum())
        lines.append(
            f"  Open in ≥{k:2d} of {n_samples} samples: {cnt:>8,}  "
            f"({100.0 * cnt / n_peaks:.1f}%)"
        )
    lines.append(
        f"  Never open (all-zero)         : "
        f"{int((open_counts == 0).sum()):>8,}  "
        f"({100.0 * (open_counts == 0).sum() / n_peaks:.1f}%)"
    )

    print("\n".join(lines))

    return {
        "cell_type"          : cell_type,
        "n_peaks"            : n_peaks,
        "n_samples"          : n_samples,
        "total_events"       : int(mat.sum()),
        "frac_open_any"      : float((open_counts > 0).mean()),
        "frac_open_all"      : float((open_counts == n_samples).mean()),
        "pearson_log1p_mean" : pr_mean,
        "spearman_mean"      : sp_mean,
    }


# --------------------------------------------------------------------------- #
# Cross-cell-type summary                                                      #
# --------------------------------------------------------------------------- #

def report_cross_cell_type(summaries: list[dict]) -> None:
    """Print a ranked cross-cell-type comparison table."""
    df = pd.DataFrame(summaries).sort_values("total_events", ascending=False)

    print(_section("Cross-cell-type summary (ranked by total events)"))
    header = (
        f"  {'Cell type':<22} {'Total events':>14} "
        f"{'Open %':>8} {'Open all %':>11} "
        f"{'Pearson r':>10} {'Spearman r':>11}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))

    for _, row in df.iterrows():
        print(
            f"  {row['cell_type']:<22} {int(row['total_events']):>14,} "
            f"{100 * row['frac_open_any']:>7.1f}% "
            f"{100 * row['frac_open_all']:>10.1f}% "
            f"{row['pearson_log1p_mean']:>10.4f} "
            f"{row['spearman_mean']:>11.4f}"
        )


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #

def main() -> None:
    cell_types = sorted(set(CELL_TYPE_STANDARDIZATION.values()))

    print(
        textwrap.dedent(f"""\
        Pseudo-bulk table validation
        ============================
        Output directory : {PSEUDOBULK_OUTPUT_DIR}
        Cell types found : {', '.join(cell_types)}
        """)
    )

    summaries: list[dict] = []

    for cell_type in cell_types:
        csv_path = PSEUDOBULK_OUTPUT_DIR / f"{cell_type}_pseudobulk.csv"

        if not csv_path.exists():
            print(f"[SKIP] {csv_path.name} not found.\n")
            continue

        df = pd.read_csv(csv_path, index_col="peak_id")
        summary = report_cell_type(cell_type, df)
        summaries.append(summary)

    if summaries:
        report_cross_cell_type(summaries)

    print("\nValidation complete.\n")


if __name__ == "__main__":
    main()
