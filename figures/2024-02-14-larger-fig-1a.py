#!/usr/bin/env python3

import gzip
import json
import os
import subprocess
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as ticker  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
from matplotlib.gridspec import GridSpec  # type: ignore

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


studies = list(metadata_papers.keys())


target_taxa = {
    9606: ("homo sapiens", "Homo Sapiens"), 
    2: ("bacteria", "Bacteria"),    
    2157: ("archaea", "Archaea"),
    10239: ("viruses", "Viruses"),
    0: ("Unassigned", "Unassigned"),
    }


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


def shape_boxplot_df(boxplot_df: pd.DataFrame) -> pd.DataFrame:
    boxplot_df = boxplot_df[
        (boxplot_df.drop(["study", "sample"], axis=1) != 0).all(1)
    ].melt(
        id_vars=["study", "sample"],
        value_vars=[
            "Homo Sapiens",
            "Bacteria",
            "Archaea",
            "Viruses",
            "Unassigned", 
            
        ],
        var_name="Identifier",
        value_name="Relative Abundance",
    )

    boxplot_df["Relative Abundance"] = boxplot_df["Relative Abundance"].apply(
        np.log10
    )

    return boxplot_df


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


def assemble_plotting_dfs() -> tuple[pd.DataFrame, pd.DataFrame]:
    box_plot_data = []
    for study in studies:
        # Dropping studies that aren't WTP based
        if study not in [
            "Bengtsson-Palme 2016",
            "Munk 2022",
            "Brinch 2020",
            "Ng 2019",
            "Maritz 2019",
            "Brumfield 2022",
            "Rothman 2021",
            "Yang 2020",
            "Spurbeck 2023",
            "Crits-Christoph 2021",
        ]:
            continue

        for bioproject in metadata_papers[study]["projects"]:
            samples = metadata_bioprojects[bioproject]

            if study == "Bengtsson-Palme 2016":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["fine_location"].startswith(
                        "Inlet"
                    )
                ]

            if study == "Ng 2019":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["fine_location"] == "Influent"
                ]

            for sample in samples:
                modified_study = study
                if metadata_samples[sample].get("enrichment") == "panel":
                    continue
                if study == "Brumfield 2022":
                    if metadata_samples[sample]["na_type"] == "DNA":
                        modified_study = "Brumfield 2022\n(DNA Subset)"

                    else:
                        modified_study = "Brumfield 2022\n(RNA Subset)"
                cladecounts = "%s.tsv.gz" % sample
                if not os.path.exists(f"../cladecounts/{cladecounts}"):
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "s3://nao-mgs/%s/cladecounts/%s"
                            % (bioproject, cladecounts),
                            "cladecounts/",
                        ]
                    )
                with gzip.open(f"../cladecounts/{cladecounts}") as inf:
                    taxa_abundances = defaultdict(float)
                    for line in inf:
                        (
                            line_taxid,
                            _,
                            _,
                            clade_assignments,
                            _,
                        ) = line.strip().split()
                        taxid = int(line_taxid)
                        clade_assignments = int(clade_assignments)
                        if taxid in target_taxa:
                            relative_abundance = (
                                clade_assignments / metadata_samples[sample]["reads"]
                            )
                            taxa_abundances[
                                target_taxa[taxid][1]
                            ] += relative_abundance

                    box_plot_data.append(
                        {
                            "study": modified_study,
                            "sample": sample,
                            **taxa_abundances,
                        }
                    )
    boxplot_df = pd.DataFrame(box_plot_data)
    boxplot_df = shape_boxplot_df(boxplot_df)
    study_nucleic_acid_mapping = get_study_nucleic_acid_mapping() 

    boxplot_df = order_df(boxplot_df, study_nucleic_acid_mapping)

    return boxplot_df


def return_study_order(boxplot_df: pd.DataFrame) -> list[str]:
    study_nucleic_acid_mapping = get_study_nucleic_acid_mapping()
    df["na_type"] = df["study"].map(study_nucleic_acid_mapping)
    order = (
        df[df["na_type"] == "DNA"]["study"].unique()
        + df[df["na_type"] == "RNA"]["study"].unique()
    )


def boxplot(
    boxplot_df: pd.DataFrame,
) -> plt.Axes:
    fig, ax = plt.subplots(figsize=(9, 11)) 
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
        # hue_order=order,
        width=0.7,
        showfliers=False,
        ax=ax,
        order=order,
    )

    ax.set_xlabel("Relative abundance among all reads")


    ax.set_ylabel("")
    ax.tick_params(left=False, labelright=True, labelleft=False)
    for label in ax.get_yticklabels():
        label.set_ha("left")

    ax.yaxis.set_label_position("right")
    formatter = ticker.FuncFormatter(
        lambda y, _: "${{10^{{{:d}}}}}$".format(int(y))
    )

    ax.xaxis.set_major_formatter(formatter)

    ax.set_xlim(right=0, left=-6)

    sns.despine(top=True, right=True, left=True, bottom=False)

    studies = boxplot_df["study"].unique()

    ax.legend(
        # location = middle below x axis. expressed in fraction of axes size
        loc=(0.25, -0.1),
        columnspacing=2.2,
        title="",
        ncol=5,
        frameon=False,
    )
    # change x labels to log scale (8 -> 10^8)

    for i in range(-5, 0):
        ax.axvline(i, color="grey", linewidth=0.3, linestyle=":")

    for i in range(1, len(studies)):
        if i == 6:
            ax.axhline(i - 0.5, color="black", linewidth=1, linestyle="-")

        else:
            ax.axhline(i - 0.5, color="grey", linewidth=0.3, linestyle=":")

    ax.text(-6.1, 0.3, "DNA \nSequencing", ha="right")
    ax.text(-6.1, 6.3, "RNA \nSequencing", ha="right")

    return fig
def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight", dpi=900)


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "figures")
    figdir.mkdir(exist_ok=True)

    boxplot_df= assemble_plotting_dfs()
    fig = boxplot(boxplot_df)



    plt.tight_layout()
    plt.show()
    save_plot(fig, figdir, "2024-02-14-larger-fig-1a")


if __name__ == "__main__":
    start()

