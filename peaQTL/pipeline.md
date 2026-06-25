# Hackathon CompGenomics: MASLD snATAC-seq Pipeline

**Author:** Dvir Mateless, Or Forshmit and Noam Kimhi  
**Date:** June 2026  

---

## Project Objective
We integrate snATAC-seq data from MASLD (Metabolic Dysfunction-Associated Steatotic Liver Disease) and healthy liver samples with eQTL summary statistics, in order to identify and model Differential Activity Variants.

---

## Data Processing Pipeline

### Step 1: Consensus Peak Generation
Since peak boundaries vary across samples, we unify them into a standard coordinate system to allow cross-sample comparisons.
* **Blacklist Filtering:** Subtract ENCODE hg38 blacklist regions from the peaks to remove technical artifacts and genomic noise.
* **Merging:** Merge overlapping intervals (allowing a 50bp distance threshold) using the `bioframe` library to create a unified feature space.

### Step 2: Split to Cell Types Using Available Annotations
Process the raw `matrix.mtx.gz` files and split them for all the different cell-types. Now each folder (e.g., `/data/snatac_gse281367/GSM8619363_MASH_rep1`) contains sub-matrices, where each matrix contains cells and peaks only from a specific cell-type.

### Step 3: Data Processing on the Cell-Specific Sub-Matrices
Process the raw matrices files:
* **Binarization:** Flatten any count value $> 0$ to $1$. This represents a pure open/closed chromatin state and neutralizes PCR duplicate noise.
* **(Optional) Cell Filtering (Columns):** Remove cells with extreme peak counts (e.g., $<1,000$ or $>10,000$) to filter out dead cells and multiplets.
* **(Optional) Peak Filtering (Rows):** Remove rare peaks appearing in fewer than 10 cells to retain statistical power and reduce dimensionality.

---

## Immediate Goal: Pseudo-Bulk Aggregation (Per Cell-Type)
We want to create count matrices representing pseudo-bulk profiles for downstream Differential Accessibility (DA) analysis. Crucially, to allow for variance calculation and statistical testing (e.g., using PyDESeq2), we must maintain the biological replicates.

### Matrix Structure (Generated separately for each Cell-Type)
* **Rows ($M$):** $peak_1, peak_2, \dots, peak_M$ (where $M$ is the size of the unified consensus peak list for this specific cell-type).
* **Columns (12 Biological Replicates):** `rep1_Normal`, `rep2_Normal`, ..., `rep6_Normal`, `rep1_MASH`, ..., `rep6_MASH`.
* **Values:** Each entry $(i, j)$ contains the aggregated sum of the binarized counts of $peak_i$ from all valid single cells belonging to replicate $j$.