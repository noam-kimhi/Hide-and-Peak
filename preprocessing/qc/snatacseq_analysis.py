import gzip
import os
import scipy.io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from constants import *

SAMPLE_DIR = MASH_REP_1_DIR
PREFIX = "GSM8619363_MASH_rep1"

# ── 1. Barcodes ──────────────────────────────────────────────────────────────
with gzip.open(os.path.join(SAMPLE_DIR, f"{PREFIX}_barcodes.tsv.gz"), 'rt') as f:
    barcodes = [line.strip() for line in f]

print(f"=== Barcodes ({len(barcodes)} cells) ===")
print("First 5:", barcodes[:5])
print()

# ── 2. Peaks (BED) ───────────────────────────────────────────────────────────
peaks = []
with gzip.open(os.path.join(SAMPLE_DIR, f"{PREFIX}_peaks.bed.gz"), 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        peaks.append(parts)

print(f"=== Peaks ({len(peaks)} peaks) ===")
print("Columns inferred: chrom, start, end (standard BED3)")
print("First 5 peaks:")
for p in peaks[:5]:
    print(" ", p)
peak_lengths = [int(p[2]) - int(p[1]) for p in peaks]
print(f"Peak length stats: min={min(peak_lengths)}, max={max(peak_lengths)}, "
      f"mean={np.mean(peak_lengths):.0f}, median={np.median(peak_lengths):.0f}")
print()

# ── 3. Sparse matrix (peaks x cells) ─────────────────────────────────────────
mat = scipy.io.mmread(os.path.join(SAMPLE_DIR, f"{PREFIX}_matrix.mtx.gz")).tocsr()
print(f"=== Matrix ===")
print(f"Shape: {mat.shape[0]} peaks x {mat.shape[1]} cells")
print(f"Non-zero entries: {mat.nnz:,}  (sparsity: {100 * (1 - mat.nnz / (mat.shape[0] * mat.shape[1])):.2f}%)")

fragments_per_cell = np.array(mat.sum(axis=0)).flatten()
peaks_per_cell     = np.array((mat > 0).sum(axis=0)).flatten()
cells_per_peak     = np.array((mat > 0).sum(axis=1)).flatten()

print(f"\nFragments per cell: min={fragments_per_cell.min():.0f}, "
      f"max={fragments_per_cell.max():.0f}, mean={fragments_per_cell.mean():.1f}, "
      f"median={np.median(fragments_per_cell):.0f}")
print(f"Peaks per cell:     min={peaks_per_cell.min():.0f}, "
      f"max={peaks_per_cell.max():.0f}, mean={peaks_per_cell.mean():.1f}, "
      f"median={np.median(peaks_per_cell):.0f}")
print(f"Cells per peak:     min={cells_per_peak.min():.0f}, "
      f"max={cells_per_peak.max():.0f}, mean={cells_per_peak.mean():.1f}, "
      f"median={np.median(cells_per_peak):.0f}")

unique_vals = np.unique(mat.data)  # .data contains only non-zero entries
print(f"\nUnique non-zero values in matrix: {sorted(unique_vals.tolist())}")
print(f"(0 is implicit for all unrecorded entries)")

print("\n=== Matrix sample (first 10 peaks x first 10 cells) ===")
sample = mat[:10, :10].toarray()
header = "peak\\cell  " + "  ".join(f"c{i:<4}" for i in range(sample.shape[1]))
print(header)
for i, row in enumerate(sample):
    print(f"peak {i:<5}" + "  ".join(f"{v:<6}" for v in row))

# ── 4. Visualizations ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 9))
fig.suptitle("GSM8619363 MASH rep1 – snATAC-seq QC", fontsize=14, fontweight='bold')

ax = axes[0, 0]
ax.hist(fragments_per_cell, bins=60, color='steelblue', edgecolor='none')
ax.set_xlabel("Total fragments per cell"); ax.set_ylabel("# cells")
ax.set_title("Fragments per cell")

ax = axes[0, 1]
ax.hist(peaks_per_cell, bins=60, color='darkorange', edgecolor='none')
ax.set_xlabel("# peaks accessible per cell"); ax.set_ylabel("# cells")
ax.set_title("Peaks per cell")

ax = axes[1, 0]
ax.hist(cells_per_peak, bins=60, color='seagreen', edgecolor='none')
ax.set_xlabel("# cells with peak accessible"); ax.set_ylabel("# peaks")
ax.set_title("Cells per peak")

ax = axes[1, 1]
ax.hist(peak_lengths, bins=60, color='mediumpurple', edgecolor='none')
ax.set_xlabel("Peak length (bp)"); ax.set_ylabel("# peaks")
ax.set_title("Peak length distribution")

ax = axes[1, 2]
vals, counts = np.unique(mat.data, return_counts=True)  # non-zero entries
n_zeros = mat.shape[0] * mat.shape[1] - mat.nnz
all_vals   = np.concatenate([[0], vals])
all_counts = np.concatenate([[n_zeros], counts])
ax.bar(all_vals, all_counts, color='tomato', edgecolor='black', linewidth=0.5)
ax.set_xticks([0, 50, 100, 150, 200, 250])
ax.set_xlabel("Matrix value"); ax.set_ylabel("Count (log scale)")
ax.set_yscale('log')
ax.set_title("Matrix value distribution (incl. zeros)")

axes[0, 2].axis('off')  # unused 3rd column top cell

plt.tight_layout()
out_path = ATAC_SEQ_PREPROCESSING_RES_DIR / "snATAC_QC.png"
plt.savefig(out_path, dpi=150)
print(f"\nPlot saved to {out_path}")
