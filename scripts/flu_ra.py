#!/usr/bin/env python3

import os
import subprocess
import csv
import pandas as pd
from scipy.stats import gmean

BIOPROJECT_DIR = "bioprojects"

if not os.path.exists(f"../bioprojects"):
    os.mkdir("../bioprojects")


TARGET_STUDY_METADATA = {
    "Rothman 2021": ["PRJNA729801"],
}


sample_files = [
    "taxonomic_composition",
    "hv_clade_counts",
    "kraken_reports",
    "qc_basic_stats",
    "sample-metadata",
]

# Downloading files from S3 in case they are not already present
for study, bioprojects in TARGET_STUDY_METADATA.items():
    for bioproject in bioprojects:
        for fname in sample_files:
            study_author = study.split()[0]
            if fname == "sample-metadata":
                if not os.path.exists(
                    f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.csv"
                ):
                    subprocess.run(
                        [
                            "aws",
                            "s3",
                            "cp",
                            f"s3://nao-mgs-wb/{study_author}-{bioproject}/output/{fname}.csv",
                            f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.csv",
                        ]
                    )
            else:
                if not os.path.exists(
                    f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.tsv"
                ):
                    subprocess.run(
                        [
                            "aws",
                            "s3",
                            "cp",
                            f"s3://nao-mgs-wb/{study_author}-{bioproject}/output/{fname}.tsv.gz",
                            f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.tsv.gz",
                        ]
                    )
                    subprocess.run(
                        [
                            "gzip",
                            "-d",
                            f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.tsv.gz",
                        ]
                    )

# Defining target taxa to cover all viruses, and DNA and RNA viruses
INFLUENZA_A_TAXID = 11320


def flu_a_level() -> pd.DataFrame:
    for study, bioprojects in TARGET_STUDY_METADATA.items():
        flu_a_ras = []
        n_samples_pos = 0
        total_samples = 0
        total_flu_reads = 0
        study_author = study.split()[0]
        for bioproject in bioprojects:
            study_bioproject = f"{study_author}-{bioproject}"
            metadata_samples = {}
            # Loading and reshaping files
            with open(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/sample-metadata.csv",
                mode="r",
                encoding="utf-8-sig",
            ) as file:
                reader = csv.DictReader(file)
                for row in reader:
                    sample = row.pop("sample")
                    metadata_samples[sample] = row
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

            for sample in samples:
                if (
                    metadata_samples[sample].get("enrichment") == "enriched"
                    or metadata_samples[sample].get("enrichment") == "1"
                ):

                    total_samples += 1
                    # Calculating influenza A relative abundance
                    total_reads = sample_read_pairs[sample]
                    flu_a_reads = hv_clade_counts[
                        (hv_clade_counts["taxid"] == INFLUENZA_A_TAXID)
                        & (hv_clade_counts["sample"] == sample)
                    ]["n_reads_clade"]
                    if flu_a_reads.empty:
                        continue
                    n_samples_pos += 1
                    flu_a_reads = flu_a_reads.iloc[0]
                    total_flu_reads += flu_a_reads
                    flu_a_ra = flu_a_reads / total_reads
                    flu_a_ras.append(flu_a_ra)
                else:
                    continue

        gmean_flu_a = gmean(flu_a_ras)
        share_pos = n_samples_pos / total_samples
        print(
            f"{study}: {gmean_flu_a}, {n_samples_pos} out of {total_samples} were positive ({share_pos:.2%})"
            f"({total_flu_reads} reads)"
        )


def start():
    flu_a_level()


if __name__ == "__main__":
    start()
