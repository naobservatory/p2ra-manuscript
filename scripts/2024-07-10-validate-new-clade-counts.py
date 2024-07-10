#!/usr/bin/env python3

import pandas as pd

# Import the data from the specified file
df = pd.read_csv(
    "bioprojects/Rothman-PRJNA729801/hv_hits_putative_filtered.tsv", sep="\t"
)


for sample in [
    "SRR14530844",
    "SRR14530845",
    "SRR18341118",
    "SRR18341008",
    "SRR18341003",
    "SRR18341011",
    "SRR18341122",
    "SRR18341120",
    "SRR18341013",
    "SRR14530833",
]:
    sample_df = df[df["sample"] == sample]
    print(sample)
    print(
        f"Number of rows in hv_hits_putative_filtered.tsv: {sample_df.shape[0]}"
    )
    unique_seq_ids = sample_df["seq_id"].nunique()
    print(f"Number of unique seq_id values: {unique_seq_ids}")
