import csv
from dataclasses import dataclass
from collections import defaultdict, namedtuple
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import gmean

MODEL_OUTPUT_DIR = "../model_output"

TABLE_OUTPUT_DIR = "../tables"


BIOPROJECT_DIR = "bioprojects"

if not os.path.exists(f"../bioprojects"):
    os.mkdir("../bioprojects")


TARGET_STUDY_METADATA = {
    "Brinch 2020": ["PRJEB13832", "PRJEB34633"],
    "Spurbeck 2023": ["PRJNA924011"],
    "CC 2021": ["PRJNA661613"],
    "Rothman 2021": ["PRJNA729801"],
}


sample_files = [
    "taxonomic_composition",
    "hv_clade_counts",
    "qc_basic_stats",
    "sample-metadata",
]

taxids = {
    11974: "Caliciviridae",
    151341: "Polyomaviridae",
    11320: "Influenza A",
    10508: "Adenoviridae",
}


def parse_hv_clade_counts():
    hv_family_data = []
    for study, bioprojects in TARGET_STUDY_METADATA.items():
        study_author = study.split()[0]
        for bioproject in bioprojects:
            study_bioproject = f"{study_author}-{bioproject}"

            metadata_samples = {}
            with open(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/sample-metadata.csv",
                mode="r",
                encoding="utf-8-sig",
            ) as file:
                reader = csv.DictReader(file)
                for row in reader:
                    sample = row.pop("sample")
                    metadata_samples[sample] = row

            fine_counts = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/kraken_reports.tsv",
                sep="\t",
            )

            fine_counts_dfs = {
                sample: df for sample, df in fine_counts.groupby("sample")
            }

            hv_clade_counts = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/hv_clade_counts.tsv",
                sep="\t",
            )

            qc_basic_stats = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/qc_basic_stats.tsv",
                sep="\t",
            )

            sample_read_pairs = dict(
                zip(qc_basic_stats["sample"], qc_basic_stats["n_read_pairs"])
            )

            samples = metadata_samples.keys()
            modified_study = study

            for sample in samples:

                # Dropping enriched samples
                if (
                    metadata_samples[sample].get("enrichment") == "enriched"
                    or metadata_samples[sample].get("enrichment") == "1"
                ):
                    continue

                sample_hv_family_counts = hv_clade_counts[
                    (hv_clade_counts["sample"] == sample)
                    & (hv_clade_counts["rank"] == "family")
                ]
                hv_family_counts_dict = dict(
                    zip(
                        sample_hv_family_counts["name"],
                        sample_hv_family_counts["n_reads_clade"],
                    )
                )

                hv_family_data.append(
                    {
                        "study": modified_study,
                        "sample": sample,
                        **hv_family_counts_dict,
                    }
                )

    hv_family_df = pd.DataFrame(hv_family_data)
    return hv_family_df


def shape_hv_family_df(hv_family_df: pd.DataFrame) -> pd.DataFrame:
    # Turn reads into relative abundances (where the denominator is the sum of all human-infecting virus family reads)
    hv_family_df = hv_family_df.groupby(["study"]).sum().drop(columns="sample")
    total_reads = hv_family_df.sum(axis=1)
    hv_family_df = hv_family_df.div(total_reads, axis=0)

    # Identify the top 9 families by relative abundance aggregated across all four studies
    N_BIGGEST_FAMILIES = 9
    total_relative_abundances = hv_family_df.sum()
    top_families = total_relative_abundances.nlargest(
        N_BIGGEST_FAMILIES
    ).index.tolist()
    # Sum the abundances of all other families for each study and add them to hv_family_df
    other_families_sum = hv_family_df.drop(columns=top_families).sum(axis=1)
    hv_family_df = hv_family_df[top_families]
    hv_family_df["Other Viral Families"] = other_families_sum
    print(hv_family_df)
    return hv_family_df


def start():
    hv_family_df = parse_hv_clade_counts()
    hv_family_df = shape_hv_family_df(hv_family_df)
    import tabulate

    print(tabulate.tabulate(hv_family_df, headers="keys", tablefmt="tsv"))


if __name__ == "__main__":
    start()
