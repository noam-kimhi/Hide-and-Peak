"""Print basic shape, columns, and example rows from the GTEx liver eQTL table.

This exploratory script reads the compressed significant variant-gene-pair file
configured by ``EQTL_DIR`` and only prints summary information.
"""

from __future__ import annotations

import gzip
from constants import EQTL_DIR

file_path = EQTL_DIR / "Liver.v8.signif_variant_gene_pairs.txt.gz"

# Print the entire table shape
with gzip.open(file_path, 'rt') as f:
    header = f.readline().strip().split('\t')
    row_count = sum(1 for _ in f)
print(f"Shape: {row_count} rows x {len(header)} columns")

with gzip.open(file_path, 'rt') as f:
    lines = [next(f) for _ in range(6)]

header = lines[0].strip().split('\t')
print("Columns:", header)
print(f"\nTotal columns: {len(header)}")
print("\nFirst 5 rows:")
for line in lines[1:]:
    fields = line.strip().split('\t')
    row = dict(zip(header, fields))
    for k, v in row.items():
        print(f"  {k}: {v}")
    print()

