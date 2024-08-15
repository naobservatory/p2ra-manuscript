#!/usr/bin/env python3

import gzip
import json
import os
import subprocess
import csv
from pathlib import Path
from scipy.stats import gmean
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as ticker  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
from matplotlib.gridspec import GridSpec  # type: ignore
from collections import defaultdict

dashboard = os.path.expanduser("~/code/mgs-pipeline/dashboard/")

DEBUG = None

with open(os.path.join(dashboard, "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

current_directory = os.getcwd()
if not current_directory.endswith("figures"):
    raise Exception(
        "Current directory is not 'figures'. Please change to the 'figures' directory to run this script."
    )


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

prettier_labels = {
    "Brinch 2020": "Brinch 2020",
    "Spurbeck 2023": "Spurbeck 2023",
    "Rothman 2021": "Rothman 2021",
    "CC 2021": "Crits-Christoph 2021",
}


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


target_taxa = {
    10239: ("viruses", "Viruses"),
    2731341: ("Duplodnaviria", "DNA Viruses"),
    2732004: ("Varidnaviria", "DNA Viruses"),
    2731342: ("Monodnaviria", "DNA Viruses"),
    2559587: ("Riboviria", "RNA Viruses"),
    9999999999: ("human viruses", "Viruses"),
}


def assemble_plotting_dfs() -> tuple[pd.DataFrame, pd.DataFrame]:
    viral_composition_data = []
    hv_family_data = []
    for study, bioprojects in TARGET_STUDY_METADATA.items():
        print(study, bioprojects)
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

                if (
                    metadata_samples[sample].get("enrichment") == "enriched"
                    or metadata_samples[sample].get("enrichment") == "1"
                ):
                    continue
                total_reads = sample_read_pairs[sample]

                total_hv_reads = hv_clade_counts[
                    (hv_clade_counts["taxid"] == 10239)
                    & (hv_clade_counts["sample"] == sample)
                ]["n_reads_clade"]

                if not total_hv_reads.empty:
                    total_hv_reads = total_hv_reads.values[0]
                else:
                    total_hv_reads = 0

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

                sample_fine_counts = fine_counts_dfs[sample]

                tax_reads_dict = sample_fine_counts.set_index("taxid")[
                    "n_reads_clade"
                ].to_dict()

                hv_relative_abundance = total_hv_reads / total_reads

                taxa_abundances = {
                    "DNA Viruses": 0,
                    "RNA Viruses": 0,
                    "Viruses": 0,
                }

                for taxid in target_taxa.keys():
                    nucleic_acid_type = target_taxa[taxid][1]
                    tax_reads = tax_reads_dict.get(taxid, 0)
                    tax_relative_abundance = tax_reads / total_reads

                    taxa_abundances[
                        nucleic_acid_type
                    ] += tax_relative_abundance

                viral_composition_data.append(
                    {
                        "study": modified_study,
                        "sample": sample,
                        "Human-Infecting Viruses": hv_relative_abundance,
                        **taxa_abundances,
                    }
                )

    viral_composition_df = pd.DataFrame(viral_composition_data)

    hv_family_df = pd.DataFrame(hv_family_data)

    viral_composition_df.to_csv("viral_composition_df.csv", index=False)
    hv_family_df.to_csv("hv_family_df.csv", index=False)

    return viral_composition_df, hv_family_df


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


def shape_vir_comp_df(viral_composition_df: pd.DataFrame) -> pd.DataFrame:

    viral_composition_df = viral_composition_df.melt(
        id_vars=["study", "sample"],
        value_vars=[
            "Viruses",
            "Human-Infecting Viruses",
            "RNA Viruses",
            "DNA Viruses",
        ],
        var_name="Group",
        value_name="Relative Abundance",
    )

    viral_composition_df["Relative Abundance"] = np.log10(
        viral_composition_df["Relative Abundance"]
    )

    return viral_composition_df


def order_df(
    df: pd.DataFrame, study_nucleic_acid_mapping: dict[str, str]
) -> pd.DataFrame:
    df = df.reset_index()
    df["na_type"] = df["study"].map(study_nucleic_acid_mapping)

    na_type_order = ["DNA", "RNA"]

    df = df.sort_values(
        by="na_type",
        key=lambda col: col.map({k: i for i, k in enumerate(na_type_order)}),
    )
    return df


def boxplot(
    ax: plt.Axes,
    viral_composition_df: pd.DataFrame,
) -> plt.Axes:

    order = [
        "Brinch 2020",
        "Spurbeck 2023",
        "Rothman 2021",
        "CC 2021",
    ]

    sns.boxplot(
        data=viral_composition_df,
        y="study",
        x="Relative Abundance",
        hue="Group",
        order=order,
        width=0.7,
        showfliers=False,
        ax=ax,
    )
    if DEBUG:
        # Calculate and print median relative abundance of human-infecting viruses for each study
        median_abundances = (
            viral_composition_df[
                viral_composition_df["Group"] == "Human-Infecting Viruses"
            ]
            .groupby("study")["Relative Abundance"]
            .median()
        )

        exp_median_abundances = median_abundances.apply(lambda x: 10**x)

        print(
            "Median relative abundance of human-infecting viruses for each study:"
        )
        for study, abundance in exp_median_abundances.items():
            print(f"{study}: {abundance:.2e}")

    ax_title = ax.set_title("a", fontweight="bold")
    ax_title.set_position((-0.15, 0))

    ax.set_xlabel("Relative abundance among all reads")

    ax.set_ylabel("")
    ax.tick_params(left=False, labelright=True, labelleft=False)
    labels = []
    for label in ax.get_yticklabels():
        pretty_label = prettier_labels[label.get_text()]
        label.set_text(pretty_label)
        labels.append(label)
    ax.set_yticklabels(labels)

    ax.yaxis.set_label_position("right")
    formatter = ticker.FuncFormatter(
        lambda y, _: "${{10^{{{:d}}}}}$".format(int(y))
    )

    ax.xaxis.set_major_formatter(formatter)

    ax.set_xlim(right=0, left=-8)

    sns.despine(top=True, right=True, left=True, bottom=False)

    studies = viral_composition_df["study"].unique()

    ax.legend(
        loc=(0.00, -0.54),
        columnspacing=2.2,
        ncol=4,
        title="",
        fontsize=10,
        frameon=False,
    )

    for i in range(-7, 0):
        ax.axvline(i, color="grey", linewidth=0.3, linestyle=":")

    for i in range(1, len(studies)):
        if i == 1:
            ax.axhline(i - 0.5, color="black", linewidth=1, linestyle="-")

        else:
            ax.axhline(i - 0.5, color="grey", linewidth=0.3, linestyle=":")

    ax.text(-8.01, 0.3, "DNA\nSequencing", ha="right")
    ax.text(-8.01, 1.3, "RNA\nSequencing", ha="right")

    return ax


def get_study_nucleic_acid_mapping() -> dict[str, str]:
    study_nucleic_acid_mapping = {
        "Brinch 2020": "DNA",
        "Spurbeck 2023": "RNA",
        "Rothman 2021": "RNA",
        "CC 2021": "RNA",
    }

    if "Brumfield 2022" in study_nucleic_acid_mapping:
        study_nucleic_acid_mapping["Brumfield 2022\n(DNA Subset)"] = "DNA"
        study_nucleic_acid_mapping["Brumfield 2022\n(RNA Subset)"] = "RNA"
        del study_nucleic_acid_mapping["Brumfield 2022"]
    return study_nucleic_acid_mapping


def return_study_order(viral_composition_df: pd.DataFrame) -> list[str]:
    study_nucleic_acid_mapping = get_study_nucleic_acid_mapping()
    viral_composition_df["na_type"] = viral_composition_df["study"].map(
        study_nucleic_acid_mapping
    )
    study_order = (
        viral_composition_df[viral_composition_df["na_type"] == "DNA"]["study"]
        .unique()
        .tolist()
        + viral_composition_df[viral_composition_df["na_type"] == "RNA"][
            "study"
        ]
        .unique()
        .tolist()
    )

    return study_order


def barplot(
    ax: plt.Axes, hv_family_df: pd.DataFrame, study_order: list
) -> plt.Axes:
    ten_color_palette = [
        "#8dd3c7",
        "#f1c232",
        "#bebada",
        "#fb8072",
        "#80b1d3",
        "#fdb462",
        "#b3de69",
        "#fccde5",
        "#bc80bd",
        "#d9d9d9",
    ]

    order = [
        "Brinch 2020",
        "Spurbeck 2023",
        "Rothman 2021",
        "CC 2021",
    ]

    hv_family_df.set_index("study", inplace=True)
    hv_family_df.loc[order].plot(
        kind="barh",
        stacked=True,
        color=ten_color_palette,
        ax=ax,
    )

    ax.invert_yaxis()
    ax_title = ax.set_title("b", fontweight="bold")
    ax_title.set_position((-0.15, 0))
    ax.set_xlabel("Relative abundance among human-infecting viruses")
    ax.set_ylabel("")
    ax.tick_params(left=False, labelright=True, labelleft=False)

    labels = []
    for label in ax.get_yticklabels():
        pretty_label = prettier_labels[label.get_text()]
        label.set_text(pretty_label)
        labels.append(label)
    ax.set_yticklabels(labels)

    ax.axhline(0.5, color="black", linewidth=1, linestyle="-")
    ax.text(-0.01, 0.2, "DNA\nSequencing", ha="right")
    ax.text(-0.01, 1.2, "RNA\nSequencing", ha="right")
    ax.set_xlim(right=1, left=0)

    ax.legend(
        loc=(0.015, -0.8),
        ncol=4,
        fontsize=9.1,
        frameon=False,
    )

    sns.despine(top=True, right=True, left=True, bottom=False)

    return ax


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight", dpi=900)


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    figdir.mkdir(exist_ok=True)
    viral_composition_df, hv_family_df = assemble_plotting_dfs()

    viral_composition_df = shape_vir_comp_df(viral_composition_df)

    hv_family_df = shape_hv_family_df(hv_family_df)

    study_nucleic_acid_mapping = get_study_nucleic_acid_mapping()

    viral_composition_df = order_df(
        viral_composition_df, study_nucleic_acid_mapping
    )
    hv_family_df = order_df(hv_family_df, study_nucleic_acid_mapping)

    fig = plt.figure(
        figsize=(9, 5),
    )

    gs = GridSpec(2, 2, height_ratios=[14, 14], figure=fig, hspace=0.8)

    boxplot(
        fig.add_subplot(gs[0, :]),
        viral_composition_df,
    )

    study_order = return_study_order(viral_composition_df)

    barplot(fig.add_subplot(gs[1, :]), hv_family_df, study_order)
    ##
    plt.tight_layout()
    # print("is relative abundance based on filtered reads as denominator?")
    # print("Clean up get_study_nucleic_acid_mapping")
    # print(
    #     "fix thing where the composition dfs are being saved and reused as CSVs"
    # )
    # print("stop usage of metadata papers file")

    save_plot(fig, figdir, "composite_fig_1")


if __name__ == "__main__":
    start()
