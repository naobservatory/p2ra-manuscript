#!/usr/bin/env python3

import gzip
import json
import os
import subprocess
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

current_directory = os.getcwd()
if not current_directory.endswith("figures"):
    raise Exception(
        "Current directory is not 'figures'. Please change to the 'figures' directory to run this script."
    )

BIOPROJECT_DIR = "bioprojects"

if not os.path.exists(f"../bioprojects"):
    os.mkdir("../bioprojects")

TARGET_STUDY_METADATA = {
    "Bengtsson-Palme 2016": ["PRJEB14051"],
    "Munk 2022": [
        "PRJEB13831",
        "PRJEB27054",
        "PRJEB27621",
        "PRJEB40798",
        "PRJEB40815",
        "PRJEB40816",
        "PRJEB51229",
    ],
    "Brinch 2020": ["PRJEB13832", "PRJEB34633"],
    "Ng 2019": ["PRJNA438174"],
    "Maritz 2019": ["PRJEB28033"],
    "Brumfield 2022": ["PRJNA812772"],
    "Yang 2020": ["PRJNA645711"],
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

# ... [Keep the file downloading code from composite_fig_1.py] ...

target_taxa = {
    2731341: ("duplodnaviria", "DNA Viruses"),
    2732004: ("varidnaviria", "DNA Viruses"),
    2731342: ("monodnaviria", "DNA Viruses"),
    2559587: ("riboviria", "RNA Viruses"),
}


def assemble_plotting_data():
    plotting_data = []
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

            kraken_reports = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/kraken_reports.tsv",
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
            sequencing_type = "DNA"  # Default to DNA

            if study == "Bengtsson-Palme 2016":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"].startswith(
                        "Inlet"
                    )
                ]

            if study == "Ng 2019":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"] == "Influent"
                ]

            for sample in samples:
                if study == "Brumfield 2022":
                    sequencing_type = metadata_samples[sample].get(
                        "na_type", "DNA"
                    )

                if (
                    metadata_samples[sample].get("enrichment") == "enriched"
                    or metadata_samples[sample].get("enrichment") == "1"
                ):
                    continue

                total_reads = sample_read_pairs[sample]

                sample_kraken = kraken_reports[
                    kraken_reports["sample"] == sample
                ]

                for taxid, (_, virus_type) in target_taxa.items():
                    tax_reads = sample_kraken[sample_kraken["taxid"] == taxid][
                        "n_reads_clade"
                    ].sum()
                    relative_abundance = tax_reads / total_reads

                    plotting_data.append(
                        {
                            "study": study,
                            "sample": sample,
                            "sequencing_type": sequencing_type,
                            "virus_na_type": virus_type,
                            "relative_abundance": relative_abundance,
                        }
                    )

    return pd.DataFrame(plotting_data)


def plot_and_analyze(df):
    df["seq_na_virus_combo"] = (
        df["sequencing_type"] + " Sequencing / " + df["virus_na_type"]
    )

    seq_na_virus_combo_ordered = [
        "DNA Sequencing / DNA Viruses",
        "RNA Sequencing / DNA Viruses",
        "DNA Sequencing / RNA Viruses",
        "RNA Sequencing / RNA Viruses",
    ]

    df["seq_na_virus_combo"] = pd.Categorical(
        df["seq_na_virus_combo"],
        categories=seq_na_virus_combo_ordered,
        ordered=True,
    )

    df = df.sort_values("seq_na_virus_combo")

    combinations = list(
        zip(seq_na_virus_combo_ordered[:-1], seq_na_virus_combo_ordered[1:])
    )

    p_values = []
    for combo1, combo2 in combinations:
        _, p_value = mannwhitneyu(
            df[df["seq_na_virus_combo"] == combo1]["relative_abundance"],
            df[df["seq_na_virus_combo"] == combo2]["relative_abundance"],
        )
        p_values.append(p_value)
        print(f"p_value when comparing '{combo1}' and '{combo2}' = {p_value}")

    df["log_relative_abundance"] = np.log10(df["relative_abundance"])

    plt.figure(figsize=(8, 4))
    sns.boxplot(x="log_relative_abundance", y="seq_na_virus_combo", data=df)

    plt.xlabel("Logged Relative Abundance")
    plt.ylabel("")
    ax = plt.gca()

    for i, p_value in enumerate(p_values):
        ax.text(
            0.2,
            i + 0.5,
            f"p < {0.0001}" if p_value < 0.0001 else f"p = {p_value:.3f}",
        )

    plt.tight_layout()
    plt.savefig("../fig/supplement_figure_1.pdf")
    plt.savefig("../fig/supplement_figure_1.png", dpi=300)


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    figdir.mkdir(exist_ok=True)

    df = assemble_plotting_data()
    plot_and_analyze(df)


if __name__ == "__main__":
    start()
