import pandas as pd
from constants import INSPECTION_DIR

CELL_METADATA_PATH = (
    INSPECTION_DIR / "GSE281367_cell_metadata.csv.gz"
)

COLUMN_SUMMARY_PATH = (
    INSPECTION_DIR / "metadata_column_summary.csv"
)

cell_metadata = pd.read_csv(
    CELL_METADATA_PATH,
    low_memory=False,
)

column_summary = pd.read_csv(COLUMN_SUMMARY_PATH)

print(f"Metadata shape: {cell_metadata.shape}")
print(cell_metadata.head())
print(column_summary)