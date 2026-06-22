import re
from constants import INSPECTION_DIR
import pandas as pd

annotation_pattern = re.compile(
    r"cell|type|annot|cluster|ident|sample|orig|"
    r"condition|disease|group|sub",
    flags=re.IGNORECASE,
)

CELL_METADATA_PATH = (
    INSPECTION_DIR / "GSE281367_cell_metadata.csv.gz"
)

cell_metadata = pd.read_csv(
    CELL_METADATA_PATH,
    low_memory=False,
)

candidate_columns = [
    column
    for column in cell_metadata.columns
    if annotation_pattern.search(column)
]

print("Candidate metadata columns:")
for column in candidate_columns:
    print(f"  {column}")

for column in candidate_columns:
    print(f"\n{'=' * 80}")
    print(column)
    print(cell_metadata[column].value_counts(dropna=False).head(30))

print("All columns:")
print(cell_metadata.columns.tolist())

cell_metadata["raw_barcode_guess"] = (
    cell_metadata["seurat_cell_id"]
    .astype(str)
    .str.extract(r"([ACGT]+-\d+)$", expand=False)
)

print(
    cell_metadata[
        ["seurat_cell_id", "raw_barcode_guess"]
    ]
    .head(30)
    .to_string(index=False)
)

print(
    "\n\nIdentifiers for which a barcode could not be extracted:",
    cell_metadata["raw_barcode_guess"].isna().sum(),
)

