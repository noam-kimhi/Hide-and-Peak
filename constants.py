from typing import Final, Sequence, Dict
from pathlib import Path

PLOT_DPI: Final[int] = 300

BASE_DIR: Final[Path] = Path.cwd()
DATA_DIR: Final[Path] = BASE_DIR / "data"

EQTL_DIR: Final[Path] = DATA_DIR / "eqtl_gtex_v8"
ATAC_SEQ_DIR: Final[Path] = DATA_DIR / "snatac_gse281367"
ATAC_SEQ_METADATA: Final[Path] = ATAC_SEQ_DIR / "GSE281367_metadata.csv"
MASH_REP_1_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619363_MASH_rep1"
MASH_REP_2_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619364_MASH_rep2"
MASH_REP_3_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619365_MASH_rep3"
MASH_REP_4_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619366_MASH_rep4"
MASH_REP_5_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619367_MASH_rep5"
MASH_REP_6_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619368_MASH_rep6"
NORMAL_REP_1_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619369_Normal_rep1"
NORMAL_REP_2_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619370_Normal_rep2"
NORMAL_REP_3_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619371_Normal_rep3"
NORMAL_REP_4_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619372_Normal_rep4"
NORMAL_REP_5_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619373_Normal_rep5"
NORMAL_REP_6_DIR: Final[Path] = ATAC_SEQ_DIR / "GSM8619374_Normal_rep6"
ATAC_SEQ_DIRS: Final[Sequence[Path]] = [
    MASH_REP_1_DIR,
    MASH_REP_2_DIR,
    MASH_REP_3_DIR,
    MASH_REP_4_DIR,
    MASH_REP_5_DIR,
    MASH_REP_6_DIR,
    NORMAL_REP_1_DIR,
    NORMAL_REP_2_DIR,
    NORMAL_REP_3_DIR,
    NORMAL_REP_4_DIR,
    NORMAL_REP_5_DIR,
    NORMAL_REP_6_DIR
]
DERIVED_DATA_DIR: Final[Path] = DATA_DIR / "derived"
DERIVED_DATA_DIR.mkdir(exist_ok=True, parents=True)

RESULTS_DIR: Final[Path] = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True, parents=True)
PREPROCESSING_RES_DIR: Final[Path] = RESULTS_DIR / "preprocessing"
PREPROCESSING_RES_DIR.mkdir(exist_ok=True, parents=True)
ATAC_SEQ_PREPROCESSING_RES_DIR: Final[Path] = PREPROCESSING_RES_DIR / "atac_seq"
ATAC_SEQ_PREPROCESSING_RES_DIR.mkdir(exist_ok=True, parents=True)

SAMPLE_QC_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "sample_qc_summary.csv"
PAIRWISE_PEAK_OVERLAP_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "pairwise_peak_overlap.csv"
SAMPLE_PEAK_SHARING_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "pairwise_peak_overlap.csv"
CONDITION_QC_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "condition_qc_summary.csv"

BARCODES_SUFFIX: Final[str] = "_barcodes.tsv.gz"
PEAKS_SUFFIX: Final[str] = "_peaks.bed.gz"
MATRIX_SUFFIX: Final[str] = "_matrix.mtx.gz"
FRAGMENTS_SUFFIX: Final[str] = "_fragments.tsv.gz"

EXPECTED_SAMPLE_COUNT: Final[int] = 12
EXPECTED_SAMPLES_PER_CONDITION: Final[Dict[str, int]] = {
    "healthy": 6,
    "MASLD": 6,
}

HEALTHY_KEYWORDS: Final[Sequence[str]] = ("healthy", "normal", "control")
MASLD_KEYWORDS: Final[Sequence[str]] = ("mash", "masld", "nash", "disease")

CONSENSUS_PEAKS_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "consensus_peaks.bed"
CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "consensus_peaks_annotations.csv"
CONSENSUS_PEAKS_SUMMARY_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "consensus_peaks_summary.csv"

# True follows the default behavior of bedtools merge:
# overlapping and directly adjacent ("book-ended") intervals are merged.
MERGE_BOOKENDED_PEAKS = True

REFERENCE_GRCH38_DIR: Final[Path] = DATA_DIR / "reference_grch38"
BLACKLIST_BED_PATH: Final[Path] = REFERENCE_GRCH38_DIR / "ENCFF356LFX.bed.gz"
REPRODUCIBLE_CONSENSUS_PEAKS_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "reproducible_consensus_peaks.bed"
REPRODUCIBLE_CONSENSUS_PEAKS_ANNOTATION_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "reproducible_consensus_peaks_annotations.csv"
PEAK_TO_CONSENSUS_MAP_OUTPUT_PATH: Final[Path] = DERIVED_DATA_DIR / "peak_to_consensus_map.csv"

MIN_REPRODUCIBLE_SAMPLE_SUPPORT = 2

RDS_PATH = ATAC_SEQ_DIR / "GSE281367_seurat_clustered.rds.gz"
PREPROCESS_DIR = BASE_DIR / "preprocessing"
SEURAT_INSPECTION_DIR = PREPROCESS_DIR / "inspecting_seurat"
INSPECTION_DIR = DERIVED_DATA_DIR / "gse281367_rds_inspection"
INSPECTION_DIR.mkdir(parents=True, exist_ok=True)
R_SCRIPT_PATH = SEURAT_INSPECTION_DIR / "inspect_seurat_object.R"
SEURAT_METADATA_PATH: Final[Path] = INSPECTION_DIR / "GSE281367_cell_metadata.csv.gz"
ANNOTATIONS_PATH: Final[Path] = INSPECTION_DIR / "GSE281367_cell_annotations.csv.gz"
ANNOTATION_DIR: Final[Path] = INSPECTION_DIR / "gse281367_annotations"
ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)
CELL_ANNOTATIONS_PATH: Final[Path] = ANNOTATION_DIR / "GSE281367_cell_annotations.csv.gz"

CELL_TYPE_STANDARDIZATION = {
    "Hepatocytes": "Hepatocytes",
    "EC": "Endothelial",
    "NK": "T_NK_B",
    "Kupffer": "Kupffer",
    "Cholangiocyte": "Cholangiocyte",
    "Stellate": "Stellate",
    "Unknow": "Unknown",
}

SAMPLE_MAPPING_PATH: Final[Path] = ANNOTATION_DIR / "sample_directory_mapping.csv"

MATCHING_ANNOTATED_COLUMNS_DIR: Final[Path] = DERIVED_DATA_DIR / "gse281367_matched_columns"

CELL_TYPE_MATRICES_DIR: Final[Path] = DERIVED_DATA_DIR / "gse281367_cell_type_matrices"
CELL_TYPE_MATRICES_DIR.mkdir(parents=True, exist_ok=True)
CELL_TYPE_MATRIX_SPLIT_SUMMARY_PATH: Final[Path] = CELL_TYPE_MATRICES_DIR / "split_summary.csv"

BINARIZED_CELL_TYPE_MATRICES_DIR: Final[Path] = DERIVED_DATA_DIR / "gse281367_cell_type_matrices_binarized"
BINARIZED_CELL_TYPE_MATRICES_DIR.mkdir(parents=True, exist_ok=True)
MATRIX_BINARIZATION_SUMMARY_PATH: Final[Path] = BINARIZED_CELL_TYPE_MATRICES_DIR / "binarization_summary.csv"

PSEUDOBULK_OUTPUT_DIR: Final[Path] = DERIVED_DATA_DIR / "pseudobulk"
PSEUDOBULK_SUMMARY_PATH: Final[Path] = PSEUDOBULK_OUTPUT_DIR / "pseudobulk_summary.csv"

PLOTS_MAIN_COLOR: Final[str] = '#f08423'
PLOTS_SECOND_COLOR: Final[str] = 'red'

MATRIX_BINARIZATION_PLOTS_DIR: Final[Path] = PREPROCESSING_RES_DIR / "matrix_splitting" / "binarization"
MATRIX_BINARIZATION_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

PEAK_FILTERING_RESULTS_DIR: Final[Path] = PREPROCESSING_RES_DIR / "peak_filtering"
PEAK_FILTERING_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PEAK_FILTER_MIN_GROUP_SIZES: Final[Sequence[int]] = (
    1,
    10,
    20,
    30,
    50,
    100,
)

PEAK_FILTER_MIN_CELL_SUPPORTS: Final[Sequence[int]] = (
    1,
    2,
    3,
    5,
    10,
)

PEAK_FILTER_MIN_CELL_FRACTIONS: Final[Sequence[float]] = (
    0.0,
    0.005,
    0.01,
    0.02,
    0.05,
)

FILTERED_CELL_TYPE_MATRICES_DIR: Final[Path] = DERIVED_DATA_DIR / "gse281367_cell_type_matrices_filtered"
FILTERED_CELL_TYPE_MATRICES_DIR.mkdir(parents=True, exist_ok=True)
PEAK_FILTERING_SUMMARY_PATH: Final[Path] = FILTERED_CELL_TYPE_MATRICES_DIR / "peak_filtering_summary.csv"
MIN_PEAK_FILTER_GROUP_SIZE: Final[int] = 30
MIN_PEAK_FILTER_CELL_SUPPORT: Final[int] = 3
MIN_PEAK_FILTER_CELL_FRACTION: Final[float] = 0.005

SIG_PEAKS_DIR: Final[Path] = DERIVED_DATA_DIR / "significant_peaks"
SIG_PEAKS_DIR.mkdir(parents=True, exist_ok=True)
UNFILTERED_SIG_PEAKS_DIR: Final[Path] = SIG_PEAKS_DIR / "unfiltered_significant_peaks"
UNFILTERED_SIG_PEAKS_DIR.mkdir(parents=True, exist_ok=True)
FILTERED_SIG_PEAKS_DIR: Final[Path] = SIG_PEAKS_DIR / "filtered_significant_peaks"
FILTERED_SIG_PEAKS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Exact GTEx liver eQTL overlap with filtered cell-type peaks
# ---------------------------------------------------------------------------

LIVER_EQTL_SIGNIFICANT_PAIRS_PATH: Final[Path] = (
    EQTL_DIR / "Liver.v8.signif_variant_gene_pairs.txt.gz"
)

UNFILTERED_SIGNIFICANT_PEAKS_SUFFIX: Final[str] = (
    "_significant_peaks.bed.gz"
)

FILTERED_SIGNIFICANT_PEAKS_SUFFIX: Final[str] = (
    "_significant_peaks_without_eqtl.bed.gz"
)

SIGNIFICANT_PEAKS_EQTL_FILTERING_SUMMARY_PATH: Final[Path] = (
    FILTERED_SIG_PEAKS_DIR
    / "significant_peaks_eqtl_filtering_summary.csv"
)

# ---------------------------------------------------------------------------
# DESeq2 differential-accessibility result analysis
# ---------------------------------------------------------------------------

PEAQTL_RESULTS_DIR: Final[Path] = RESULTS_DIR / "peaQTL"

DESEQ2_RESULTS_DIR = PEAQTL_RESULTS_DIR / "differential_peaks"
DESEQ2_RESULTS_SUFFIX = "_deseq2_results.csv"

DESEQ2_ANALYSIS_DIR = DESEQ2_RESULTS_DIR / "analysis"
DESEQ2_SUMMARY_CSV = (
    DESEQ2_ANALYSIS_DIR / "deseq2_all_cell_types_summary.csv"
)
DESEQ2_THRESHOLD_GRID_CSV = (
    DESEQ2_ANALYSIS_DIR / "deseq2_significance_threshold_grid.csv"
)

DESEQ2_PLOTS_DIR = DESEQ2_ANALYSIS_DIR / "plots"
DESEQ2_PER_CELL_PLOTS_DIR = DESEQ2_PLOTS_DIR / "per_cell_type"

# Thresholds examined in the summary tables.
DESEQ2_PVALUE_THRESHOLDS = (
    0.05,
    0.01,
    0.001,
    0.0001,
)

DESEQ2_PADJ_THRESHOLDS = (
    0.10,
    0.05,
    0.01,
)

DESEQ2_ABS_LOG2FC_THRESHOLDS = (
    0.0,
    0.25,
    0.50,
    0.58,
    1.00,
)

DESEQ2_BASE_MEAN_THRESHOLDS = (
    1.0,
    5.0,
    10.0,
)

# Default criterion used only for highlighting points and summary plots.
# Final threshold selection can be changed after inspecting the results.
DESEQ2_DEFAULT_PADJ_THRESHOLD = 0.05
DESEQ2_DEFAULT_ABS_LOG2FC_THRESHOLD = 0.58

# Plot settings.
DESEQ2_PLOT_FORMAT = "png"
DESEQ2_PVALUE_HISTOGRAM_BINS = 40
DESEQ2_TOP_PEAK_LABELS = 8
DESEQ2_MAX_NEG_LOG10_FOR_PLOTS = 50.0
DESEQ2_MA_MAX_ABS_LOG2FC = 10.0

DESEQ2_THRESHOLD_PLOT_ABS_LOG2FC_THRESHOLDS: Final[Sequence[float]] = (
    0.0,
    0.58,
    1.0,
)

DESEQ2_TOP_FOREST_PEAKS: Final[int] = 20

DESEQ2_POSITIVE_LFC_LABEL: Final[str] = "More accessible in MASLD"
DESEQ2_NEGATIVE_LFC_LABEL: Final[str] = "More accessible in healthy"

# ---------------------------------------------------------------------------
# Sample-level support for differential peaks
# ---------------------------------------------------------------------------

DESEQ2_SAMPLE_SUPPORT_DIR: Final[Path] = (
    DESEQ2_ANALYSIS_DIR / "sample_support"
)

DESEQ2_SAMPLE_SUPPORT_PLOTS_DIR: Final[Path] = (
    DESEQ2_SAMPLE_SUPPORT_DIR / "plots"
)

DESEQ2_SAMPLE_SUPPORT_PER_CELL_PLOTS_DIR: Final[Path] = (
    DESEQ2_SAMPLE_SUPPORT_PLOTS_DIR / "per_cell_type"
)

DESEQ2_SAMPLE_SUPPORT_SUMMARY_CSV: Final[Path] = (
    DESEQ2_SAMPLE_SUPPORT_DIR
    / "differential_peak_sample_support.csv"
)

DESEQ2_SAMPLE_SUPPORT_SIZE_FACTORS_CSV: Final[Path] = (
    DESEQ2_SAMPLE_SUPPORT_DIR
    / "pseudobulk_visualization_size_factors.csv"
)

# These samples must not participate in any statistics or plots.
DESEQ2_EXCLUDED_SAMPLE_TOKENS: Final[Sequence[str]] = (
    "Normal_rep4",
    "Normal_rep6",
)

DESEQ2_EXPECTED_HEALTHY_SAMPLE_COUNT: Final[int] = 4
DESEQ2_EXPECTED_MASLD_SAMPLE_COUNT: Final[int] = 6

# Number of candidates included in the detailed plots.
DESEQ2_SAMPLE_SUPPORT_TOP_PEAKS: Final[int] = 12
DESEQ2_SAMPLE_SUPPORT_HEATMAP_PEAKS: Final[int] = 30
DESEQ2_SAMPLE_SUPPORT_MAX_PEAK_LABELS: Final[int] = 8

# Used only when calculating a descriptive log2 ratio of group means.
DESEQ2_SAMPLE_SUPPORT_PSEUDOCOUNT: Final[float] = 0.5

# Descriptive consistency thresholds. These do not replace padj.
DESEQ2_SAMPLE_SUPPORT_STRONG_THRESHOLD: Final[float] = 0.80
DESEQ2_SAMPLE_SUPPORT_MODERATE_THRESHOLD: Final[float] = 0.65

DESEQ2_SAMPLE_SUPPORT_LOO_STRONG_THRESHOLD: Final[float] = 1.00
DESEQ2_SAMPLE_SUPPORT_LOO_MODERATE_THRESHOLD: Final[float] = 0.80

# Flag a group when one sample contributes more than this fraction of
# the group's total normalized signal for a peak.
DESEQ2_SAMPLE_SUPPORT_MAX_SAMPLE_SHARE_WARNING: Final[float] = 0.50