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
    "Bengtsson-Palme": {
        "bioprojects": ["PRJEB14051"],
        "year": 2016,
    },
    "Munk": {
        "bioprojects": ["PRJEB13831", "PRJEB27054", "PRJEB27621", "PRJEB40798", "PRJEB40815", "PRJEB40816", "PRJEB51229"],
        "year": 2022,
    },
    "Brinch": {
        "bioprojects": ["PRJEB13832", "PRJEB34633"],
        "year": 2020,
    },
    "Ng": {
        "bioprojects": ["PRJNA438174"],
        "year": 2019,
    },
    "Maritz": {
        "bioprojects": ["PRJEB28033"],
        "year": 2019,
    },
    "Brumfield": {
        "bioprojects": ["PRJNA812772"],
        "year": 2022,
    },
    "Yang": {
        "bioprojects": ["PRJNA645711"],
        "year": 2020,
    },
    "Spurbeck": {
        "bioprojects": ["PRJNA924011"],
        "year": 2023,
    },
    "CC": {
        "bioprojects": ["PRJNA661613"],
        "year": 2021,
    },
    #"Rothman": { # not yet run through the pipeline
    #    "bioprojects": ["PRJNA729801"],
    #    "year": 2021,
    #},
}


taxonomic_files = ["taxonomic_composition", "hv_clade_counts", "kraken_reports"]
metadata_file = "sample-metadata"


for study, metadata_study in TARGET_STUDY_METADATA.items():
    for bioproject in metadata_study["bioprojects"]:
        for fname in taxonomic_files:
            if not os.path.exists(f"../taxonomy/{study}-{bioproject}/{fname}.tsv"):
                subprocess.run(['aws', 's3', 'cp', f's3://nao-mgs-wb/{study}-{bioproject}/output/{fname}.tsv.gz', f"../taxonomy/{study}-{bioproject}/{fname}.tsv.gz"])
                subprocess.run(['gzip', '-d', f"../taxonomy/{study}-{bioproject}/{fname}.tsv.gz"])
        if not os.path.exists(f"../taxonomy/{study}-{bioproject}/{metadata_file}.csv"):
            subprocess.run(['aws', 's3', 'cp', f's3://nao-mgs-wb/{study}-{bioproject}/output/{metadata_file}.csv', f"../taxonomy/{study}-{bioproject}/{metadata_file}.csv"])



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
    #2842242: ("ribozyviria", "RNA"),  # not in hv_clade_counts
    #687329: ("anelloviridae", "DNA"), # not in hv_clade_counts
    #2840022: ("adnaviria", "DNA"), # not in hv_clade_counts
    9999999999: ("human viruses", "Viruses"),
}




def assemble_plotting_dfs() -> tuple[pd.DataFrame, pd.DataFrame]:
    viral_composition_data = []
    hv_genus_data = []
    for study, metadata_study in TARGET_STUDY_METADATA.items():
        for bioproject in metadata_study["bioprojects"]:
            study_bioproject = f"{study}-{bioproject}"

            metadata_samples = {}
            with open(f"../taxonomy/{study_bioproject}/sample-metadata.csv", mode='r', encoding='utf-8-sig') as file:
                reader = csv.DictReader(file)
                for row in reader:
                    sample = row.pop('sample')
                    metadata_samples[sample] = row


            fine_taxonomy = pd.read_csv(f"../taxonomy/{study_bioproject}/kraken_reports.tsv", sep="\t")
            fine_taxonomy_dict = {sample: df for sample, df in fine_taxonomy.groupby("sample")}


            high_level_taxonomy = pd.read_csv(f"../taxonomy/{study_bioproject}/taxonomic_composition.tsv", sep="\t")

            hv_clade_counts = pd.read_csv(f"../taxonomy/{study_bioproject}/hv_clade_counts.tsv", sep="\t")



            samples = metadata_samples.keys()

            if study == "Bengtsson-Palme":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"].startswith(
                        "Inlet"
                    )
                ]

            if study == "Ng":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"] == "Influent"
                ]


            hv_virus_species = hv_clade_counts[hv_clade_counts["rank"] == "species"]["rank"].unique()

            for sample in samples:

                if metadata_samples[sample].get("enrichment") == "panel":
                    continue
                sample_hv_clade_counts = hv_clade_counts[hv_clade_counts["sample"] == sample]
                hv_clade_counts_dict = dict(zip(sample_hv_clade_counts['taxid'], sample_hv_clade_counts['n_reads_clade']))
                hv_genus_data.append({
                    "study": study,
                    "sample": sample,
                    **hv_clade_counts_dict
                })


                total_filtered_reads = high_level_taxonomy[(high_level_taxonomy["sample"] == sample) & (high_level_taxonomy["classification"] == "Filtered")]["n_reads"].iloc[0]

                sample_fine_taxonomy = fine_taxonomy_dict[sample]

                tax_relative_abundance_dict = sample_fine_taxonomy.set_index("taxid")["pc_reads_total"].to_dict()

                hv_reads = hv_clade_counts_dict.get(10239, 0) # if no human infecting viruses, return 0
                hv_relative_abundance = hv_reads / total_filtered_reads


                taxa_abundances = {
                    "DNA Viruses": 0,
                    "RNA Viruses": 0,
                    "Viruses": 0,
                }

                for taxid in target_taxa.keys():
                    nucleic_acid_type = target_taxa[taxid][1]
                    tax_relative_abundance_perc = tax_relative_abundance_dict.get(taxid, 0)
                    tax_relative_abundance = tax_relative_abundance_perc / 100 # % to fraction

                    taxa_abundances[nucleic_acid_type] += tax_relative_abundance

                viral_composition_data.append(
                    {
                        "study": study,
                        "sample": sample,
                        "Human-Infecting Viruses": hv_relative_abundance,
                        **taxa_abundances,
                    }
                )

    viral_composition_df = pd.DataFrame(viral_composition_data)
    #viral_composition_df = shape_vir_comp_df(viral_composition_df)

    hv_genus_df = pd.DataFrame(hv_genus_data)
    #hv_genus_df = shape_hv_genus_df(hv_genus_df)

    #study_nucleic_acid_mapping = get_study_nucleic_acid_mapping()

    #viral_composition_df = order_df(viral_composition_df, study_nucleic_acid_mapping)
    #hv_genus_df = order_df(hv_genus_df, study_nucleic_acid_mapping)


    return viral_composition_df, hv_genus_df




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


def shape_vir_comp_df(viral_composition_df: pd.DataFrame) -> pd.DataFrame:
    viral_composition_df = viral_composition_df[
        (viral_composition_df.drop(["study", "sample"], axis=1) != 0).all(1)
    ].melt(
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

    viral_composition_df["Relative Abundance"] = viral_composition_df["Relative Abundance"].apply(
        np.log10
    )

    return viral_composition_df


def shape_hv_genus_df(hv_genus_df: pd.DataFrame) -> pd.DataFrame:
    taxid_parents = load_taxonomic_data()

    unique_numeric_cols = set(
        col for col in hv_genus_df.columns if isinstance(col, int)
    )

    species_family_mapping = {
        col: get_family(col, taxid_parents) for col in unique_numeric_cols
    }

    hv_genus_df.rename(columns=species_family_mapping, inplace=True)

    melted_df = hv_genus_df.melt(
        id_vars=["study", "sample"],
        value_vars=[
            col
            for col in hv_genus_df.columns
            if col not in ["study", "sample"] and not pd.isna(col)
        ],
        var_name="family_taxid",
        value_name="read_count",
    )

    grouped_df = melted_df.groupby(["study", "family_taxid"]).read_count.sum()

    grouped_df_w_o_zeroes = grouped_df[grouped_df != 0]

    df_normalized = grouped_df_w_o_zeroes.unstack(level=-1).fillna(0)

    df_normalized = df_normalized.div(df_normalized.sum(axis=1), axis=0)

    N_BIGGEST_FAMILIES = 9

    family_mean_across_studies = df_normalized.apply(
        lambda x: np.mean(x), axis=0
    )

    top_families = family_mean_across_studies.nlargest(
        N_BIGGEST_FAMILIES
    ).index

    df_normalized["Other Viral Families"] = df_normalized.loc[
        :, ~df_normalized.columns.isin(top_families)
    ].sum(axis=1)

    df_normalized = df_normalized.loc[
        :, list(top_families) + ["Other Viral Families"]
    ]
    dict_family_name = {}

    for taxid in top_families.tolist():
        dict_family_name[taxid] = get_taxid_name(taxid, taxonomic_names)

    barplot_df = df_normalized.rename(columns=dict_family_name)

    return barplot_df


def get_study_nucleic_acid_mapping() -> dict[str, str]:
    study_nucleic_acid_mapping = {
        study: metadata["na_type"]
        for study, metadata in metadata_papers.items()
    }

    if "Brumfield 2022" in study_nucleic_acid_mapping:
        study_nucleic_acid_mapping["Brumfield 2022\n(DNA Subset)"] = "DNA"
        study_nucleic_acid_mapping["Brumfield 2022\n(RNA Subset)"] = "RNA"
        del study_nucleic_acid_mapping["Brumfield 2022"]
    return study_nucleic_acid_mapping


def load_taxonomic_data() -> dict[int, tuple[str, int]]:
    parents = {}
    with open(os.path.join(dashboard, "nodes.dmp")) as inf:
        for line in inf:
            child_taxid, parent_taxid, child_rank, *_ = line.replace(
                "\t|\n", ""
            ).split("\t|\t")
            parent_taxid = int(parent_taxid)
            child_taxid = int(child_taxid)
            child_rank = child_rank.strip()
            parents[child_taxid] = (child_rank, parent_taxid)
    return parents


def get_family(taxid: int, parents: dict[int, tuple[str, int]]) -> int:
    iteration_count = 0

    original_taxid = taxid

    current_taxid, current_rank, parent_taxid = (
        original_taxid,
    ) + parents.get(taxid, (None, None))

    while current_rank != "family" and parent_taxid is not None:
        current_taxid = parent_taxid
        current_rank, parent_taxid = parents.get(current_taxid, (None, None))
        iteration_count += 1

        if iteration_count > 100:
            break
    else:
        family_taxid = current_taxid
        return family_taxid


def get_taxid_name(
    target_taxid: int, taxonomic_names: dict[str, list[str]]
) -> str:
    tax_name = taxonomic_names[f"{target_taxid}"][0]
    return tax_name


def return_study_order(boxplot_df: pd.DataFrame) -> list[str]:
    study_nucleic_acid_mapping = get_study_nucleic_acid_mapping()
    boxplot_df["na_type"] = boxplot_df["study"].map(study_nucleic_acid_mapping)
    order = (
        boxplot_df[boxplot_df["na_type"] == "DNA"]["study"].unique()
        + boxplot_df[boxplot_df["na_type"] == "RNA"]["study"].unique()
    )


def boxplot(
    ax: plt.Axes,
    boxplot_df: pd.DataFrame,
) -> plt.Axes:

    print(boxplot_df)
    order = [
        "Bengtsson-\nPalme 2016",
        "Munk 2022",
        "Brinch 2020",
        "Ng 2019",
        "Maritz 2019",
        "Brumfield 2022\n(DNA Subset)",
        "Brumfield 2022\n(RNA Subset)",
        "Rothman 2021",
        "Yang 2020",
        "Spurbeck 2023",
        "Crits-\nChristoph 2021",
    ]

    sns.boxplot(
        data=boxplot_df,
        y="study",
        x="Relative Abundance",
        hue="Identifier",
        hue_order=order,
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

    studies = boxplot_df["study"].unique()

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


def barplot(
    ax: plt.Axes, barplot_df: pd.DataFrame, study_order: list
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

    barplot_df.set_index("study", inplace=True)

    barplot_df.loc[study_order].plot(
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

    boxplot_df, barplot_df = assemble_plotting_dfs()

    #fig = plt.figure(
    #    figsize=(9, 11),
    #)
##
    #gs = GridSpec(2, 2, height_ratios=[9, 7], figure=fig)
##
    #boxplot_ax = boxplot(
    #    fig.add_subplot(gs[0, :]),
    #    boxplot_df,
    #)
##
    #study_order = [text.get_text() for text in boxplot_ax.get_yticklabels()]
##
    #barplot(fig.add_subplot(gs[1, :]), barplot_df, study_order)
##
    #plt.tight_layout()
    print("Watch out, this doesn't include Rothman atm.")
    print("Watch out, Bengtsson-Palme has no information on RNA vs DNA.")
    print("You still need to fix the unknown virus issue in hv_clade_counts.")
    print("is relative abundance based on filtered reads as denominator?")

    #save_plot(fig, figdir, "composite_fig_1")


if __name__ == "__main__":
    start()
