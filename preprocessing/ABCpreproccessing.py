import logging
import pandas as pd
from constants import ABC_PRED_FILE, ABC_DATA_DIR

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = ["chr", "start", "end", "TargetGene", "CellType", "ABC.Score"]
OUTPUT_PATH = ABC_DATA_DIR / "Stanford_ABC_Liver_Dictionary.csv"
ABC_SCORE_THRESHOLD = 0.02
LIVER_KEYWORDS = ("liver", "hepatocyte")


def load_abc_predictions(path) -> pd.DataFrame:
    logger.info("Loading file: %s", path)
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=REQUIRED_COLUMNS,
        low_memory=False,
    )
    logger.info("Loaded %d rows.", len(df))
    return df


def filter_liver_cell_types(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Filtering for liver/hepatocyte cell types (case-insensitive)...")
    pattern = "|".join(LIVER_KEYWORDS)
    mask = df["CellType"].str.contains(pattern, case=False, na=False, regex=True)
    filtered = df[mask].copy()
    logger.info(
        "Rows after tissue filter: %d (removed %d).",
        len(filtered),
        len(df) - len(filtered),
    )
    return filtered


def filter_by_abc_score(df: pd.DataFrame) -> pd.DataFrame:
    logger.info(
        "Filtering rows with ABC.Score >= %.4f...", ABC_SCORE_THRESHOLD
    )
    filtered = df[df["ABC.Score"] >= ABC_SCORE_THRESHOLD].copy()
    logger.info(
        "Rows after score filter: %d (removed %d).",
        len(filtered),
        len(df) - len(filtered),
    )
    return filtered


def rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Renaming 'chr' -> 'chrom' for bioframe compatibility...")
    return df.rename(columns={"chr": "chrom"})


def save_output(df: pd.DataFrame, path) -> None:
    logger.info("Saving %d rows to: %s", len(df), path)
    df.to_csv(path, index=False)
    logger.info("Done.")


def main() -> None:
    df = load_abc_predictions(ABC_PRED_FILE)
    df = filter_liver_cell_types(df)
    df = filter_by_abc_score(df)
    df = rename_columns(df)
    save_output(df, OUTPUT_PATH)


if __name__ == "__main__":
    main()
