from typing import Final, Sequence, Dict
from pathlib import Path

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

MIN_REPRODUCIBLE_SAMPLE_SUPPORT = 2

RDS_PATH = ATAC_SEQ_DIR / "GSE281367_seurat_clustered.rds.gz"
PREPROCESS_DIR = BASE_DIR / "preprocessing"
SEURAT_INSPECTION_DIR = PREPROCESS_DIR / "inspecting_seurat"
INSPECTION_DIR = DERIVED_DATA_DIR / "gse281367_rds_inspection"
INSPECTION_DIR.mkdir(parents=True, exist_ok=True)
R_SCRIPT_PATH = SEURAT_INSPECTION_DIR / "inspect_seurat_object.R"
ANNOTATIONS_PATH: Final[Path] = INSPECTION_DIR / "GSE281367_cell_metadata.csv.gz"
ANNOTATION_DIR: Final[Path] = INSPECTION_DIR / "gse281367_annotations"
ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)

CELL_TYPE_STANDARDIZATION = {
    "Hepatocytes": "Hepatocytes",
    "EC": "Endothelial",
    "NK": "T_NK_B",
    "Kupffer": "Kupffer",
    "Cholangiocyte": "Cholangiocyte",
    "Stellate": "Stellate",
    "Unknow": "Unknown",
}
