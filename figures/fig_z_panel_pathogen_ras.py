#!/usr/bin/env pyth n3

import gzip
import json
import os
import subprocess
from collections import defaultdict

import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as ticker  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
from pathlib import Path

# set print to print all columns
pd.set_option("display.max_columns", None)
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
    2697049: ("SARS-CoV-2"),
    122928: ("Norovirus GI"),
    122929: ("Norovirus GII"),
    11676: ("HIV"),
    493803: ("MCV"),
    }

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


def assemble_plotting_dfs() -> tuple[pd.DataFrame]:
    box_plot_data = []
    for study in studies:
        # Dropping studies that aren't WTP based
        if study not in [
            "Rothman 2021",
            "Spurbeck 2023",
            "Crits-Christoph 2021",
        ]:
            continue

        for bioproject in metadata_papers[study]["projects"]:
            samples = metadata_bioprojects[bioproject]
            for sample in samples:
                modified_study = study
                #if metadata_samples[sample].get("enrichment") == "panel":
                #    continue
                if study in ["Rothman 2021", "Crits-Christoph 2021"]:
                    if metadata_samples[sample]["enrichment"] == "panel":
                        modified_study = f"{study}\n(Panel enriched)"
                    else:
                        modified_study = f"{study}\n(Unenriched)"
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
                                target_taxa[taxid]
                            ] += relative_abundance

                    box_plot_data.append(
                        {
                            "Study": modified_study,
                            "sample": sample,
                            **taxa_abundances,
                        }
                    )
    #print(box_plot_data)
    df = pd.DataFrame(box_plot_data)


    return df

def shape_df(df: pd.DataFrame) -> pd.DataFrame:
    
    df = df.melt(
        id_vars=["Study", "sample"],
        value_vars=[
            "SARS-CoV-2",
            "Norovirus GI",
            "Norovirus GII",
            "HIV",
            "MCV",
            
        ],
        var_name="Pathogen",
        value_name="Relative Abundance",
    )
    
    df["Study"] = df["Study"].str.replace("-", "-\n")
    df["Relative Abundance"] = df["Relative Abundance"].apply(
        np.log10
    )
    df = df.fillna(0)
    print(df.head(20))
    df = df[df["Relative Abundance"] != -np.inf]
    df = df[df["Relative Abundance"] != np.inf]
    df = df[df["Relative Abundance"] != 0]

    return df

def compute_zero_shares(df: pd.DataFrame) -> pd.DataFrame: #FIX FIX
    df = df.melt(
    id_vars=["Study", "sample"],
    value_vars=[
        "SARS-CoV-2",
        "Norovirus GI",
        "Norovirus GII",
        "HIV",
        "MCV",
        
    ],
    var_name="Pathogen",
    value_name="Relative Abundance",
    )
    df.fillna(0)
    # replace inf and -inf with 0
    df = df.replace([np.inf, -np.inf], 0)

    zero_shares = df.assign(is_zero=df["Relative Abundance"] == 0).groupby(["Pathogen", "Study"])["is_zero"].mean()
    print(zero_shares) # FIX FIX FIX





def barplot(
        df: pd.DataFrame,
    )-> plt.Axes:
    sns.barplot(
        data=df,
        y="Study",
        x="Zero Share",
        hue="Pathogen",
        ax=ax,
    )

def boxplot(
    df: pd.DataFrame,
) -> plt.Axes:
    fig, ax = plt.subplots(figsize=(9, 4))
    #sns.boxplot(
    #    data=df,
    #    y="Study",
    #    x="Relative Abundance",
    #    hue="Pathogen",
    #    width=0.7, 
    #    showfliers=False,
    #    ax=ax,
    #)
    
    sns.stripplot(
        data=df,
        y="Study",
        x="Relative Abundance",
        hue="Pathogen",
        dodge=True,
        ax=ax,
        alpha=0.5,
    )
    ax.set_xlabel("Relative abundance among all reads", fontsize=12)

    ax.set_ylabel("")
    ax.tick_params(left=False, labelright=True, labelleft=False)
    for label in ax.get_yticklabels():
        label.set_ha("left")

    ax.yaxis.set_label_position("right")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    formatter = ticker.FuncFormatter(
        lambda y, _: "${{10^{{{:d}}}}}$".format(int(y))
    )

    ax.xaxis.set_major_formatter(formatter)

    ax.set_xlim(right=0, left=-8)

    sns.despine(top=True, right=True, left=True, bottom=False)

    studies = df["Study"].unique()

    ax.legend(
        loc=(0.25, -0.2),
        columnspacing=2.2,
        title="",
        ncol=5,
        frameon=False,
    )
    xmin, xmax = ax.get_xlim()
    for i in range(0, -8, -1):
        ax.axvline(i, color="grey", linewidth=0.3, linestyle=":")

    for i in range(1, len(studies)):
        if i in [2,4]:
            ax.axhline(i - 0.5, color="darkgrey", linewidth=1, linestyle="--")
        else:
            ax.axhline(i - 0.5, color="grey", linewidth=0.3, linestyle=":")

    

    return fig
def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight", dpi=900)


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "figures")
    figdir.mkdir(exist_ok=True)

    df= assemble_plotting_dfs()
    compute_zero_shares(df)
    df = shape_df(df) 
    fig = boxplot(df)



    plt.tight_layout()
    plt.show()
    save_plot(fig, figdir, "fig-1c")


if __name__ == "__main__":
    start()

