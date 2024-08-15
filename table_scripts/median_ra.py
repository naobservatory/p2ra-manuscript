#!/usr/bin/env python3

import os
import subprocess
import csv
import pandas as pd


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
target_taxa = {
    10239: ("viruses", "Viruses"),
    2731341: ("Duplodnaviria", "DNA Viruses"),
    2732004: ("Varidnaviria", "DNA Viruses"),
    2731342: ("Monodnaviria", "DNA Viruses"),
    2559587: ("Riboviria", "RNA Viruses"),
}


def assemble_df() -> pd.DataFrame:
    viral_composition_data = []
    for study, bioprojects in TARGET_STUDY_METADATA.items():
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
                total_reads = sample_read_pairs[sample]

                # Calculating human-infecting virus abundance
                total_hv_reads = hv_clade_counts[
                    (hv_clade_counts["taxid"] == 10239)
                    & (hv_clade_counts["sample"] == sample)
                ]["n_reads_clade"]

                if not total_hv_reads.empty:
                    total_hv_reads = total_hv_reads.values[0]
                else:
                    total_hv_reads = 0

                hv_relative_abundance = total_hv_reads / total_reads

                # Calculating abundances of RNA, DNA and total viruses
                sample_fine_counts = fine_counts_dfs[sample]

                tax_reads_dict = sample_fine_counts.set_index("taxid")[
                    "n_reads_clade"
                ].to_dict()

                taxa_abundances = {
                    "DNA Viruses": 0,
                    "RNA Viruses": 0,
                    "Viruses": 0,
                }

                for taxid in target_taxa.keys():
                    virus_group = target_taxa[taxid][1]

                    tax_reads = tax_reads_dict.get(taxid, 0)
                    tax_relative_abundance = tax_reads / total_reads

                    taxa_abundances[virus_group] += tax_relative_abundance

                viral_composition_data.append(
                    {
                        "study": modified_study,
                        "sample": sample,
                        "Human-Infecting Viruses": hv_relative_abundance,
                        **taxa_abundances,
                    }
                )

    viral_composition_df = pd.DataFrame(viral_composition_data)

    return viral_composition_df


def shape_vir_comp_df(viral_composition_df: pd.DataFrame) -> pd.DataFrame:

    viral_composition_df = viral_composition_df.melt(
        id_vars=["study", "sample"],
        value_vars=[
            "Viruses",
            "Human-Infecting Viruses",
            "RNA Viruses",
            "DNA Viruses",
        ],
        var_name="group",
        value_name="Relative Abundance",
    )

    return viral_composition_df


def start():
    viral_composition_df = assemble_df()
    viral_composition_df = shape_vir_comp_df(viral_composition_df)

    medians = {}

    for study in viral_composition_df["study"].unique():
        medians[study] = {}
        for group in viral_composition_df["group"].unique():
            medians[study][group] = viral_composition_df[
                (viral_composition_df["study"] == study)
                & (viral_composition_df["group"] == group)
            ]["Relative Abundance"].median()
        print(study + ":")
        for group in medians[study].keys():
            print(group, f"{medians[study][group]:.2e}")
        print("\n")


if __name__ == "__main__":
    start()
