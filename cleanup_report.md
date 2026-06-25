# Cleanup report

## 1. Final package structure

`preprocessing/` now contains `abc/`, `eqtl/`, `matrix_splitting/`, `peak_filtering/`, `peaks/`, `qc/`, and `seurat/`. `peaQTL/` now contains `differential_accessibility/`, `eqtl_filtering/`, `peak_to_gene/`, and `significant_peaks/`. Each subpackage has `__init__.py` and groups modules by retained pipeline stage.

## 2. Files moved or renamed

| Old path | New path | Reason | Updated references |
|---|---|---|---|
| preprocessing/snATACseq_analysis.py | preprocessing/qc/snatacseq_analysis.py | sample/input QC package | module command/report |
| preprocessing/process_atac_seq.py | preprocessing/qc/process_atac_seq.py | snATAC QC/lineage package | module command/report |
| preprocessing/create_concensus_peaks.py | preprocessing/peaks/create_concensus_peaks.py | peak construction package | module command/report |
| preprocessing/filter_cell_type_peaks.py | preprocessing/peaks/filter_cell_type_peaks.py | peak filtering package | module command/report |
| preprocessing/process_eQTLs.py | preprocessing/eqtl/process_eqtls.py | eQTL inspection package | module command/report |
| preprocessing/ABCpreproccessing.py | preprocessing/abc/preprocess_abc.py | ABC preprocessing package and clearer spelling | module command/report |
| preprocessing/handle_metadata.py | preprocessing/seurat/handle_metadata.py | Seurat/metadata package | module command/report |
| preprocessing/inspecting_seurat/* | preprocessing/seurat/* | Seurat inspection package | constants.R_SCRIPT_PATH; documentation strings |
| peaQTL/find_differential_peaks.py | peaQTL/differential_accessibility/find_differential_peaks.py | differential-accessibility package | module commands |
| peaQTL/analyze_deseq2_results.py | peaQTL/differential_accessibility/analyze_deseq2_results.py | DESeq2 diagnostics package | module commands |
| peaQTL/plot_deseq2_results.py | peaQTL/differential_accessibility/plot_deseq2_results.py | DESeq2 plotting package | module commands |
| peaQTL/plot_differential_peak_sample_support.py | peaQTL/differential_accessibility/plot_differential_peak_sample_support.py | sample-support diagnostics package | module commands |
| peaQTL/find_significant_peaks_per_ct.py | peaQTL/significant_peaks/find_significant_peaks_per_ct.py | significant-peak selection package | module commands |
| peaQTL/find_soft_significant_peaks_per_ct.py | peaQTL/significant_peaks/find_soft_significant_peaks_per_ct.py | soft significant-peak selection package | module commands |
| peaQTL/drop_eqtl_from_significant_peaks.py | peaQTL/eqtl_filtering/drop_eqtl_from_significant_peaks.py | eQTL overlap removal package | docstring commands |
| peaQTL/create_peaks_without_eqt.py | peaQTL/eqtl_filtering/create_peaks_without_eqt.py | eQTL filtering package | module commands |
| peaQTL/peak2gene.py | peaQTL/peak_to_gene/peak2gene.py | peak-to-gene package | docstring commands |
| peaQTL/create_consensus_peaks_bioframe.py | peaQTL/peak_to_gene/create_consensus_peaks_bioframe.py | peak-to-gene support package | module commands |
| peaQTL/create_table.py | peaQTL/peak_to_gene/create_table.py | peak-to-gene table package | module commands |
| peaQTL/validate_tables.py | peaQTL/peak_to_gene/validate_tables.py | peak-to-gene validation package | sys.path root update |


## 3. Files deleted

| Deleted path | Category | Producer/source | Consumers checked | Why deletion is safe | Regeneration command |
|---|---|---|---|---|---|
| peaQTL/enrichment_analysis.py | Enrichment source | manual enrichment scripts | repository-wide rg and doc command checks | exclusive discarded enrichment workflow | not retained |
| peaQTL/masld_gene_enrichment.py | Enrichment source | manual enrichment scripts | repository-wide rg and doc command checks | exclusive discarded enrichment workflow | not retained |
| results/peaQTL/masld_enrichment/masld_enrichment_default.csv | Enrichment result | removed enrichment modules | repository-wide rg; no retained readers | output of discarded enrichment workflow | deleted workflow only |
| results/peaQTL/masld_enrichment/masld_enrichment_old.csv | Enrichment result | removed enrichment modules | repository-wide rg; no retained readers | output of discarded enrichment workflow | deleted workflow only |
| results/peaQTL/masld_enrichment/masld_enrichment_soft.csv | Enrichment result | removed enrichment modules | repository-wide rg; no retained readers | output of discarded enrichment workflow | deleted workflow only |
| results/peaQTL/pathway_enrichment/Endothelial_significant_pathways.csv | Enrichment result | removed enrichment modules | repository-wide rg; no retained readers | output of discarded enrichment workflow | deleted workflow only |
| results/peaQTL/pathway_enrichment/Hepatocytes_significant_pathways.csv | Enrichment result | removed enrichment modules | repository-wide rg; no retained readers | output of discarded enrichment workflow | deleted workflow only |
| results/peaQTL/pathway_enrichment/Kupffer_significant_pathways.csv | Enrichment result | removed enrichment modules | repository-wide rg; no retained readers | output of discarded enrichment workflow | deleted workflow only |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/Sample_support_plot_review.md | Redundant result | peaQTL.differential_accessibility.plot_differential_peak_sample_support | canonical files remain in sibling per_cell_type directories | review/export duplicate would become stale after color changes | python -m peaQTL.differential_accessibility.plot_differential_peak_sample_support |

## 4. Data-lineage audit

| Path or pattern | Role | Producer/source | Consumers | Kept/deleted | Reason |
|---|---|---|---|---|---|
| data/snatac_gse281367/** | Raw/source snATAC matrices, fragments, peaks, metadata, Seurat object | External GSE281367 download | preprocessing.qc, preprocessing.seurat, matrix splitting, peak construction | Kept | Raw/reference data; no approval requested |
| data/snrna_reference_gse189600/** | Raw/reference snRNA matrices and metadata | External GSE189600 download | not conclusively consumed by retained scripts | Kept uncertain | Potential reference data; no deletion without approval |
| data/reference_grch38/** | Reference GTF/blacklist | External references | preprocessing.peaks.create_concensus_peaks | Kept | Required reference inputs |
| data/ABC/** | ABC source/reference and lifted dictionary | External ABC, preprocessing.abc.preprocess_abc | peaQTL.peak_to_gene.peak2gene | Kept | Required peak-to-gene inputs |
| data/eqtl_gtex_v8/** | GTEx eQTL reference | External GTEx | peaQTL.eqtl_filtering | Kept | Required eQTL filtering input |
| data/derived/gse281367_rds_inspection/** | Derived Seurat metadata/annotations | preprocessing.seurat | matrix splitting and annotation checks | Kept | Retained downstream inputs |
| data/derived/gse281367_cell_type_matrices/** | Derived split matrices | preprocessing.matrix_splitting.split_matrices_by_cell_type | binarization | Kept | Retained downstream input |
| data/derived/gse281367_cell_type_matrices_binarized/** | Derived binary matrices | preprocessing.matrix_splitting.binarize_cell_type_matrices | peak filtering | Kept | Retained downstream input |
| data/derived/gse281367_cell_type_matrices_filtered/** | Derived filtered matrices | preprocessing.peaks.filter_cell_type_peaks | pseudobulk / differential accessibility | Kept | Retained downstream input |
| data/derived/pseudobulk/** | Final derived count tables | retained pseudobulk workflow | peaQTL.differential_accessibility.find_differential_peaks and sample support plotting | Kept | Required downstream input |
| data/derived/significant_peaks/** | Derived significant BED outputs | peaQTL.significant_peaks and eQTL filtering | peak-to-gene mapping | Kept | Required downstream input |
| data/derived/consensus* and reproducible_consensus* | Derived peak references | preprocessing.peaks | filtering / reproducibility | Kept | Required or uncertain retained intermediates |

No non-enrichment file under `data/` was deleted. Potentially unused non-enrichment data were retained pending explicit approval.

## 5. Enrichment removal

Removed scripts: `peaQTL/enrichment_analysis.py`, `peaQTL/masld_gene_enrichment.py`. Removed outputs: `results/peaQTL/pathway_enrichment/*` and `results/peaQTL/masld_enrichment/*`. Intended dependency cleanup: `gseapy` is still present in `requirements.txt` for this text-only PR because Git reports that UTF-16 file as binary; remove it in the binary-file follow-up. No enrichment-specific data files under `data/` were found and deleted. The remaining word "enrichment" in `analyze_deseq2_results.py` is a descriptive statistical column name for p-value excess over a uniform null, not an active enrichment workflow, and was retained to preserve output schema.

## 6. Unused code removal

Removed only conclusively unused enrichment modules and redundant review/export result metadata. `gseapy` remains in `requirements.txt` for this text-only PR because Git reports that file as binary, but repository-wide inspection found it was imported only by the deleted Enrichr script. Other uncertain source modules were retained.

## 7. Plot-color changes

Healthy: #076913
MASLD:   #BC490D

Changed `constants.py` to define shared `HEALTHY_PLOT_COLOR` and `MASLD_PLOT_COLOR`. Updated `peaQTL.differential_accessibility.plot_deseq2_results` so positive MASLD-direction significant peaks use MASLD and negative Healthy-direction peaks use Healthy. Updated `peaQTL.differential_accessibility.plot_differential_peak_sample_support` so sample distributions, heatmap condition annotations, and legends use the shared colors. Canonical plots require regeneration, but current result CSV schemas in this checkout prevented successful regeneration.

## 8. Math-rendering changes

Representative changes: `log2 fold change` -> `$\log_2(\mathrm{FC})$`; `-log10(raw p-value)` -> `$-\log_{10}(p)$`; `|log2FC| ≥ threshold` -> `$|\log_2(\mathrm{FC})| \geq threshold$`; `log2(normalized count + 1)` -> `$\log_2(\mathrm{normalized\ count} + 1)$`.

## 9. Results regenerated

No retained plot outputs were regenerated successfully. Commands attempted are recorded in validation; they failed because tracked result CSV files in this checkout lack the columns expected by retained plotting/analysis code.

## 10. Validation

Baseline: `git status --short`, `git branch --show-current`, and `git rev-parse HEAD` were recorded before edits; branch was `work` at `67e3e11cf20586dfc534a5b5ee3e3345a7791d60`. Baseline `python -m compileall preprocessing peaQTL` passed and no `tests/` directory existed. Post-change `python -m compileall preprocessing peaQTL` passed. Plot regeneration failed on missing expected result columns. Stale-reference searches were run for enrichment terms, old module paths, and old color literals. Data-safety confirmation: no non-enrichment `data/` paths were deleted or moved.

## Binary-file follow-up

The branch was compared against `67e3e11cf20586dfc534a5b5ee3e3345a7791d60` (`HEAD~2` at the time of this follow-up) because no remote/base branch is configured in this checkout. Git reported `requirements.txt` as a binary modification, and the generated `.png` review-copy deletions are binary-file-path changes that are not suitable for a text-only Codex PR. These changes were reverted to their base-branch state in this branch. They still need to be handled after the text-only PR merges.

| Path | Intended action | Reason | Exact command needed after merging the text-only PR |
|---|---|---|---|
| requirements.txt | replace/edit | Remove enrichment-only `gseapy`, but this UTF-16 requirements file is reported as binary by Git/Codex PR tooling. | `python - <<'PY'\nfrom pathlib import Path\np=Path('requirements.txt')\ns=p.read_text(encoding='utf-16')\ns='\\n'.join(line for line in s.splitlines() if not line.startswith('gseapy'))+'\\n'\np.write_text(s, encoding='utf-16')\nPY` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Cholangiocyte_effect_vs_sample_support.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Cholangiocyte_effect_vs_sample_support.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Cholangiocyte_top_peak_sample_distributions.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Cholangiocyte_top_peak_sample_distributions.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Cholangiocyte_top_peak_sample_heatmap.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Cholangiocyte_top_peak_sample_heatmap.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Endothelial_effect_vs_sample_support.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Endothelial_effect_vs_sample_support.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Endothelial_top_peak_sample_distributions.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Endothelial_top_peak_sample_distributions.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Endothelial_top_peak_sample_heatmap.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Endothelial_top_peak_sample_heatmap.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Hepatocytes_effect_vs_sample_support.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Hepatocytes_effect_vs_sample_support.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Hepatocytes_top_peak_sample_distributions.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Hepatocytes_top_peak_sample_distributions.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Hepatocytes_top_peak_sample_heatmap.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Hepatocytes_top_peak_sample_heatmap.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Kupffer_effect_vs_sample_support.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Kupffer_effect_vs_sample_support.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Kupffer_top_peak_sample_distributions.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Kupffer_top_peak_sample_distributions.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Kupffer_top_peak_sample_heatmap.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/Kupffer_top_peak_sample_heatmap.png` |
| results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/all_cell_types_candidate_sample_support.png | delete | Redundant generated review/export copy; canonical per-cell-type plot remains and plotting code changes are retained. | `git rm results/peaQTL/differential_peaks/analysis/sample_support/plots/sample_support_plot_review/images/all_cell_types_candidate_sample_support.png` |

## 11. Remaining uncertainties

Potentially unused but retained pending approval: `data/snrna_reference_gse189600/**`, `results/peaQTL/peak2gene_before_hg38_fix/**`, and any derived `data/` intermediate not conclusively proven unused. The review-copy image deletion was reverted for the text-only PR; no `data/` files were deleted.

## A. Required reruns after this cleanup

1. `python -m peaQTL.differential_accessibility.find_differential_peaks` — inputs: `data/derived/pseudobulk/`; outputs: `results/peaQTL/differential_peaks/*_deseq2_results.csv`; reason: restore result tables with expected schema; Codex ran successfully: no.
2. `python -m peaQTL.differential_accessibility.analyze_deseq2_results` — inputs: DESeq2 CSVs; outputs: `results/peaQTL/differential_peaks/analysis/*.csv`; reason: plotting summaries depend on these tables; Codex ran successfully: no, failed on current stale schema.
3. `python -m peaQTL.differential_accessibility.plot_deseq2_results` — inputs: DESeq2 CSVs and analysis summaries; outputs: `results/peaQTL/differential_peaks/analysis/plots/`; reason: refresh Healthy/MASLD colors and math text; Codex ran successfully: no, blocked by stale schema.
4. `python -m peaQTL.differential_accessibility.plot_differential_peak_sample_support` — inputs: DESeq2 CSVs and `data/derived/pseudobulk/`; outputs: `results/peaQTL/differential_peaks/analysis/sample_support/`; reason: refresh Healthy/MASLD colors and math text; Codex ran successfully: no, blocked by stale schema.
5. `python -m peaQTL.significant_peaks.find_significant_peaks_per_ct` — inputs: DESeq2 CSVs and pseudobulk; outputs: `data/derived/significant_peaks/unfiltered_significant_peaks/`; reason: update downstream selected peaks if DESeq2 results were restored; Codex ran successfully: no.
6. `python -m peaQTL.eqtl_filtering.drop_eqtl_from_significant_peaks --overwrite` — inputs: significant peaks and GTEx eQTLs; outputs: filtered significant peaks; reason: update eQTL-filtered peaks; Codex ran successfully: no.
7. `python -m peaQTL.peak_to_gene.peak2gene` — inputs: filtered peaks and ABC dictionary; outputs: `results/peaQTL/peak2gene/`; reason: update peak-to-gene outputs after upstream reruns; Codex ran successfully: no.

## B. Full pipeline order from scratch

1. `python -m preprocessing.seurat.create_inspection_script`; inputs: constants and Seurat RDS path; outputs: `preprocessing/seurat/inspect_seurat_object.R`; prerequisite: R/Seurat available; mandatory for Seurat export.
2. `Rscript preprocessing/seurat/inspect_seurat_object.R`; inputs: `data/snatac_gse281367/GSE281367_seurat_clustered.rds.gz`; outputs: `data/derived/gse281367_rds_inspection/`; prerequisite: R packages; mandatory.
3. `python -m preprocessing.seurat.recreate_standardized_annotations`; inputs: Seurat inspection output; outputs: standardized annotations; mandatory.
4. `python -m preprocessing.seurat.validate_annotation_barcode_matching`; inputs: annotations and raw barcodes; outputs: validation summaries; recommended QC.
5. `python -m preprocessing.matrix_splitting.split_matrices_by_cell_type`; inputs: raw sample matrices and annotations; outputs: `data/derived/gse281367_cell_type_matrices/`; mandatory.
6. `python -m preprocessing.matrix_splitting.binarize_cell_type_matrices`; inputs: split matrices; outputs: `data/derived/gse281367_cell_type_matrices_binarized/`; mandatory.
7. `python -m preprocessing.matrix_splitting.report_matrix_unique_values`; inputs: split/binarized matrices; outputs: matrix unique-value report; QC optional.
8. `python -m preprocessing.matrix_splitting.plot_binarization_distributions`; inputs: split/binarized matrices; outputs: binarization plots; QC optional.
9. `python -m preprocessing.peaks.create_concensus_peaks`; inputs: raw peaks and blacklist; outputs: consensus/reproducible peaks; mandatory.
10. `python -m preprocessing.peak_filtering.analyze_peak_filtering_hyperparameters`; inputs: binarized matrices; outputs: `results/preprocessing/peak_filtering/`; diagnostic optional.
11. `python -m preprocessing.peaks.filter_cell_type_peaks`; inputs: binarized matrices; outputs: filtered matrices and summaries; mandatory.
12. `python -m preprocessing.abc.preprocess_abc`; inputs: ABC predictions, chain, liftOver; outputs: `data/ABC/Stanford_ABC_Liver_Dictionary.csv`; mandatory for peak-to-gene if dictionary absent.
13. `python -m peaQTL.differential_accessibility.find_differential_peaks`; inputs: pseudobulk tables; outputs: DESeq2 CSVs; mandatory.
14. `python -m peaQTL.differential_accessibility.analyze_deseq2_results`; inputs: DESeq2 CSVs; outputs: analysis summaries; mandatory for plots.
15. `python -m peaQTL.differential_accessibility.plot_deseq2_results`; inputs: DESeq2 and analysis summaries; outputs: DESeq2 plots; mandatory for refreshed plots.
16. `python -m peaQTL.differential_accessibility.plot_differential_peak_sample_support`; inputs: DESeq2 and pseudobulk; outputs: sample-support diagnostics; mandatory for refreshed plots.
17. `python -m peaQTL.significant_peaks.find_significant_peaks_per_ct`; inputs: DESeq2 and pseudobulk; outputs: unfiltered significant peaks; mandatory.
18. `python -m peaQTL.significant_peaks.find_soft_significant_peaks_per_ct`; inputs: DESeq2; outputs: soft significant peaks; optional soft workflow.
19. `python -m peaQTL.eqtl_filtering.drop_eqtl_from_significant_peaks --overwrite`; inputs: significant peaks and GTEx eQTLs; outputs: eQTL-filtered peaks; mandatory before strict peak-to-gene.
20. `python -m peaQTL.peak_to_gene.peak2gene`; inputs: eQTL-filtered peaks and ABC dictionary; outputs: `results/peaQTL/peak2gene/`; mandatory final mapping.
21. `python -m peaQTL.peak_to_gene.peak2gene --soft`; inputs: soft filtered peaks and ABC dictionary; outputs: `results/peaQTL/soft_peak2gene/`; optional soft mapping.
