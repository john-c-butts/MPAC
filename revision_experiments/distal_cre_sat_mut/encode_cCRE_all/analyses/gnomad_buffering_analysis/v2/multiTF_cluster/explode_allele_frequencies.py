#!/usr/bin/env python3
"""
Explode Open Targets alleleFrequencies into separate columns per population.
Output: variantId + one column per population allele frequency.
"""

import pyarrow.parquet as pq
import pandas as pd
from pathlib import Path

# Input directory with parquet files
VARIANT_DIR = Path("/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/raw_data/openTargets/variant/")

# Output file
OUTPUT_FILE = Path(__file__).parent / "openTargets_allele_frequencies.tsv.gz"


def explode_allele_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert alleleFrequencies list of dicts to separate columns.

    Input: DataFrame with variantId and alleleFrequencies columns
    Output: DataFrame with variantId and one column per population AF
    """
    records = []

    for _, row in df.iterrows():
        record = {'variantId': row['variantId']}
        af_list = row['alleleFrequencies']

        if af_list is not None:
            for entry in af_list:
                if entry is not None:
                    pop_name = entry.get('populationName')
                    af_value = entry.get('alleleFrequency')
                    if pop_name is not None:
                        # Create column name like 'AF_afr_adj'
                        record[f'AF_{pop_name}'] = af_value

        records.append(record)

    return pd.DataFrame(records)


def main():
    # Get all parquet files
    parquet_files = sorted(VARIANT_DIR.glob("*.parquet"))
    print(f"Found {len(parquet_files)} parquet files")

    # Process each file and collect results
    all_dfs = []

    for i, pf_path in enumerate(parquet_files):
        print(f"Processing {i+1}/{len(parquet_files)}: {pf_path.name}")

        # Read only the columns we need
        table = pq.read_table(pf_path, columns=['variantId', 'alleleFrequencies'])
        df = table.to_pandas()

        # Explode allele frequencies
        exploded = explode_allele_frequencies(df)
        all_dfs.append(exploded)

        # Progress indicator
        if (i + 1) % 5 == 0:
            print(f"  Processed {i+1} files...")

    # Combine all dataframes
    print("Combining all dataframes...")
    result = pd.concat(all_dfs, ignore_index=True)

    # Reorder columns: variantId first, then sorted AF columns
    af_cols = sorted([c for c in result.columns if c.startswith('AF_')])
    result = result[['variantId'] + af_cols]

    print(f"\nFinal shape: {result.shape}")
    print(f"Columns: {result.columns.tolist()}")

    # Save to compressed TSV
    print(f"\nSaving to {OUTPUT_FILE}...")
    result.to_csv(OUTPUT_FILE, sep='\t', index=False, compression='gzip')
    print("Done!")


if __name__ == "__main__":
    main()
