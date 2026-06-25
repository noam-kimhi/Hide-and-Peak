"""Remove single-valued columns from the snATAC metadata CSV.

The script reads the configured GSE281367 metadata file and writes a cleaned
CSV beside it. Running this script modifies a file under ``data/``.
"""

from __future__ import annotations

import pandas as pd
from constants import ATAC_SEQ_METADATA, ATAC_SEQ_DIR

def main() -> None:
    """
    Load metadata, drop columns with at most one unique value, and write CSV.
    """
    # import /data/snatac_gse281367/GSE281367_metadata.csv 
    # if a column doesn't have more then one unique value, drop it
    df = pd.read_csv(ATAC_SEQ_METADATA)
    for col in df.columns:
        if df[col].nunique() <= 1:
            df.drop(col, axis=1, inplace=True)
    # save the new dataframe to /data/snatac_gse281367/GSE281367_metadata_cleaned.csv
    df.to_csv(ATAC_SEQ_DIR / 'GSE281367_metadata_cleaned.csv', index=False)
    


if __name__ == '__main__':
    main()