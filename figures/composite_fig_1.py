#!/usr/bin/env python3

import gzip
import json
import os
import subprocess
import csv
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as ticker  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
from matplotlib.gridspec import GridSpec  # type: ignore
from collections import defaultdict

dashboard = os.path.expanduser("~/code/mgs-pipeline/dashboard/")


with open(os.path.join(dashboard, "human_virus_sample_counts.json")) as inf:
    human_virus_sample_counts = json.load(inf)

with open(os.path.join(dashboard, "metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.join(dashboard, "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)

with open(os.path.join(dashboard, "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

with open(os.path.join(dashboard, "taxonomic_names.json")) as inf:
    taxonomic_names = json.load(inf)

if not os.path.exists("../taxonomy"):
    os.mkdir("../taxonomy")


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
    # "Rothman 2021": ["PRJNA729801"], # not yet run through the pipeline
}
#
#    "Bengtsson-Palme 2016": {
#        "bioprojects": ["PRJEB14051"],
#    },
#    "Munk 2022": {
#        "bioprojects": ["PRJEB13831", "PRJEB27054", "PRJEB27621",
#        "PRJEB40798", "PRJEB40815", "PRJEB40816", "PRJEB51229"],
#    },
#    "Brinch 2020": {
#        "bioprojects": ["PRJEB13832", "PRJEB34633"],
#    },
#    "Ng 2019": {
#        "bioprojects": ["PRJNA438174"],
#    },
#    "Maritz 2019": {
#        "bioprojects": ["PRJEB28033"],
#    },
#    "Brumfield 2022": {
#        "bioprojects": ["PRJNA812772"],
#    },
#    "Yang 2020": {
#        "bioprojects": ["PRJNA645711"],
#    },
#    "Spurbeck 2023": {
#        "bioprojects": ["PRJNA924011"],
#    },
#    "CC 2021": {
#        "bioprojects": ["PRJNA661613"],
#    },
#    #"Rothman 2021": { # not yet run through the pipeline
#    #    "bioprojects": ["PRJNA729801"],
#

sample_files = [
    "taxonomic_composition",
    "hv_clade_counts",
    "kraken_reports",
    "qc_basic_stats",
    "sample-metadata",
]


for study, bioprojects in TARGET_STUDY_METADATA.items():
    for bioproject in bioprojects:
        for fname in sample_files:
            study_author = study.split()[0]
            if fname == "sample-metadata":
                if not os.path.exists(
                    f"../taxonomy/{study_author}-{bioproject}/{fname}.csv"
                ):
                    subprocess.run(
                        [
                            "aws",
                            "s3",
                            "cp",
                            f"s3://nao-mgs-wb/{study_author}-{bioproject}/output/{fname}.csv",
                            f"../taxonomy/{study_author}-{bioproject}/{fname}.csv",
                        ]
                    )
            else:
                if not os.path.exists(
                    f"../taxonomy/{study_author}-{bioproject}/{fname}.tsv"
                ):
                    subprocess.run(
                        [
                            "aws",
                            "s3",
                            "cp",
                            f"s3://nao-mgs-wb/{study_author}-{bioproject}/output/{fname}.tsv.gz",
                            f"../taxonomy/{study_author}-{bioproject}/{fname}.tsv.gz",
                        ]
                    )
                    subprocess.run(
                        [
                            "gzip",
                            "-d",
                            f"../taxonomy/{study_author}-{bioproject}/{fname}.tsv.gz",
                        ]
                    )


current_directory = os.getcwd()
if not current_directory.endswith("figures"):
    raise Exception(
        "Current directory is not 'figures'. Please change to the 'figures' directory to run this script."
    )


target_taxa = {
    10239: ("viruses", "Viruses"),
    2731341: ("Duplodnaviria", "DNA Viruses"),
    2732004: ("Varidnaviria", "DNA Viruses"),
    2731342: ("Monodnaviria", "DNA Viruses"),
    2559587: ("Riboviria", "RNA Viruses"),
    # 2842242: ("ribozyviria", "RNA"),  # not in hv_clade_counts
    # 687329: ("anelloviridae", "DNA"), # not in hv_clade_counts
    # 2840022: ("adnaviria", "DNA"), # not in hv_clade_counts
    9999999999: ("human viruses", "Viruses"),
}


def assemble_plotting_dfs() -> tuple[pd.DataFrame, pd.DataFrame]:
    viral_composition_data = []
    hv_family_data = []
    for study, bioprojects in TARGET_STUDY_METADATA.items():
        study_author = study.split()[0]
        for bioproject in bioprojects:
            study_bioproject = f"{study_author}-{bioproject}"

            metadata_samples = {}
            with open(
                f"../taxonomy/{study_bioproject}/sample-metadata.csv",
                mode="r",
                encoding="utf-8-sig",
            ) as file:
                reader = csv.DictReader(file)
                for row in reader:
                    sample = row.pop("sample")
                    metadata_samples[sample] = row

            fine_taxonomy = pd.read_csv(
                f"../taxonomy/{study_bioproject}/kraken_reports.tsv", sep="\t"
            )
            fine_taxonomy_dfs = {
                sample: df for sample, df in fine_taxonomy.groupby("sample")
            }

            hv_clade_counts = pd.read_csv(
                f"../taxonomy/{study_bioproject}/hv_clade_counts.tsv", sep="\t"
            )

            qc_basic_stats = pd.read_csv(
                f"../taxonomy/{study_bioproject}/qc_basic_stats.tsv", sep="\t"
            )

            sample_read_pairs = dict(
                zip(qc_basic_stats["sample"], qc_basic_stats["n_read_pairs"])
            )

            samples = metadata_samples.keys()

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

                if metadata_samples[sample].get("enrichment") == "panel":
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
                    {"study": study, "sample": sample, **hv_family_counts_dict}
                )

                sample_fine_taxonomy = fine_taxonomy_dfs[sample]

                tax_reads_dict = sample_fine_taxonomy.set_index("taxid")[
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
                        "study": study,
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

    N_BIGGEST_FAMILIES = 9

    hv_family_df = hv_family_df.groupby(["study"]).sum().drop(columns="sample")
    hv_family_df = hv_family_df.div(hv_family_df.sum(axis=1), axis=0)

    top_families = (
        hv_family_df.apply(lambda x: np.mean(x), axis=0)
        .nlargest(N_BIGGEST_FAMILIES)
        .index.tolist()
    )
    other_families_sum = hv_family_df.loc[
        :, ~hv_family_df.columns.isin(top_families)
    ].sum(axis=1)
    hv_family_df = hv_family_df[top_families]

    hv_family_df["Other Viral Families"] = other_families_sum

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
        var_name="Identifier",
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

    df["study"] = df["study"].str.replace("-", "-\n")

    return df


def boxplot(
    ax: plt.Axes,
    viral_composition_df: pd.DataFrame,
) -> plt.Axes:

    # order = [
    #    "Bengtsson-\nPalme 2016",
    #    "Munk 2022",
    #    "Brinch 2020",
    #    "Ng 2019",
    #    "Maritz 2019",
    #    "Brumfield 2022\n(DNA Subset)",
    #    "Brumfield 2022\n(RNA Subset)",
    #    "Rothman 2021",
    #    "Yang 2020",
    #    "Spurbeck 2023",
    #    "Crits-\nChristoph 2021",
    # ]

    order = [
        "Bengtsson-\nPalme 2016",
        "Munk 2022",
        "Brinch 2020",
        "Ng 2019",
        "Maritz 2019",
        "Brumfield 2022",
        # "Brumfield 2022\n(RNA Subset)",
        # "Rothman 2021",
        "Yang 2020",
        "Spurbeck 2023",
        "CC 2021",
    ]

    sns.boxplot(
        data=viral_composition_df,
        y="study",
        x="Relative Abundance",
        hue="Identifier",
        order=order,
        width=0.7,
        showfliers=False,
        ax=ax,
    )

    ax_title = ax.set_title("a", fontweight="bold")
    ax.set_xlabel("Relative abundance among all reads")

    ax_title.set_position((-0.15, 0))

    ax.set_ylabel("")
    ax.tick_params(left=False, labelright=True, labelleft=False)
    for label in ax.get_yticklabels():
        label.set_ha("left")

    ax.yaxis.set_label_position("right")
    formatter = ticker.FuncFormatter(
        lambda y, _: "${{10^{{{:d}}}}}$".format(int(y))
    )

    ax.xaxis.set_major_formatter(formatter)

    ax.set_xlim(right=0, left=-8)

    sns.despine(top=True, right=True, left=True, bottom=False)

    studies = viral_composition_df["study"].unique()

    ax.legend(
        loc=(0.00, -0.17),
        columnspacing=2.2,
        ncol=4,
        title="",
        fontsize=10,
        frameon=False,
    )
    # change x labels to log scale (8 -> 10^8)

    for i in range(-7, 0):
        ax.axvline(i, color="grey", linewidth=0.3, linestyle=":")

    for i in range(1, len(studies)):
        if i == 6:
            ax.axhline(i - 0.5, color="black", linewidth=1, linestyle="-")

        else:
            ax.axhline(i - 0.5, color="grey", linewidth=0.3, linestyle=":")

    ax.text(-8.1, 0.3, "DNA \nSequencing", ha="right")
    ax.text(-8.1, 6.3, "RNA \nSequencing", ha="right")

    return ax


def get_study_nucleic_acid_mapping() -> dict[str, str]:
    study_nucleic_acid_mapping = {
        study: metadata["na_type"]
        for study, metadata in metadata_papers.items()
    }

    # if "Brumfield 2022" in study_nucleic_acid_mapping:
    #    study_nucleic_acid_mapping["Brumfield 2022\n(DNA Subset)"] = "DNA"
    #    study_nucleic_acid_mapping["Brumfield 2022\n(RNA Subset)"] = "RNA"
    #    del study_nucleic_acid_mapping["Brumfield 2022"]
    return study_nucleic_acid_mapping


def return_study_order(viral_composition_df: pd.DataFrame) -> list[str]:
    study_nucleic_acid_mapping = get_study_nucleic_acid_mapping()
    viral_composition_df["na_type"] = viral_composition_df["study"].map(
        study_nucleic_acid_mapping
    )
    order = (
        viral_composition_df[viral_composition_df["na_type"] == "DNA"][
            "study"
        ].unique()
        + viral_composition_df[viral_composition_df["na_type"] == "RNA"][
            "study"
        ].unique()
    )


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

    hv_family_df.set_index("study", inplace=True)

    hv_family_df.loc[study_order].plot(
        kind="barh",
        stacked=True,
        color=ten_color_palette,
        ax=ax,
    )

    ax.invert_yaxis()

    ax_title = ax.set_title("b", fontweight="bold")

    ax_title.set_position((-0.15, 0))

    ax.set_xlabel("Relative abundance among human-infecting viruses")

    ax.tick_params(left=False)

    ax.set_ylabel("")
    ax.tick_params(left=False, labelright=True, labelleft=False)
    for label in ax.get_yticklabels():
        label.set_ha("left")

    ax.axhline(5.5, color="black", linewidth=1, linestyle="-")

    ax.text(-0.01, 0.5, "DNA \nSequencing", ha="right")
    ax.text(-0.01, 6.5, "RNA \nSequencing", ha="right")

    ax.set_xlim(right=1, left=0)

    ax.legend(
        loc=(0.015, -0.32),
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
    figdir = Path(parent_dir / "figures")
    figdir.mkdir(exist_ok=True)

    # Load the DataFrames from CSV files if they exist #FIXME
    if os.path.exists("viral_composition_df.csv") and os.path.exists(
        "hv_family_df.csv"
    ):
        viral_composition_df = pd.read_csv("viral_composition_df.csv")
        hv_family_df = pd.read_csv("hv_family_df.csv")
    else:
        viral_composition_df, hv_family_df = assemble_plotting_dfs()

    viral_composition_df = shape_vir_comp_df(viral_composition_df)

    hv_family_df = shape_hv_family_df(hv_family_df)

    study_nucleic_acid_mapping = get_study_nucleic_acid_mapping()

    viral_composition_df = order_df(
        viral_composition_df, study_nucleic_acid_mapping
    )
    hv_family_df = order_df(hv_family_df, study_nucleic_acid_mapping)

    fig = plt.figure(
        figsize=(9, 11),
    )
    ##
    gs = GridSpec(2, 2, height_ratios=[9, 7], figure=fig)
    ##
    boxplot_ax = boxplot(
        fig.add_subplot(gs[0, :]),
        viral_composition_df,
    )

    study_order = [text.get_text() for text in boxplot_ax.get_yticklabels()]

    barplot(fig.add_subplot(gs[1, :]), hv_family_df, study_order)
    ##
    plt.tight_layout()
    print("Watch out, this doesn't include Rothman atm.")
    print("Watch out, Brumfield has not been split.")
    print("You still need to fix the unknown virus issue in hv_clade_counts.")
    print("is relative abundance based on filtered reads as denominator?")
    print("Clean up get_study_nucleic_acid_mapping")
    print(
        "fix thing where the composition dfs are being saved and reused as CSVs"
    )

    save_plot(fig, figdir, "composite_fig_1")


if __name__ == "__main__":
    start()
