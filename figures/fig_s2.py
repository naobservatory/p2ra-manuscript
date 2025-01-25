import csv
import os
from collections import defaultdict
from dataclasses import dataclass
from math import log
from pathlib import Path

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore

PERCENTILES = [5, 25, 50, 75, 95]

MODEL_OUTPUT_DIR = "../model_output"

import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42


def fits_df() -> pd.DataFrame:
    data: defaultdict[str, list] = defaultdict(list)

    for p in PERCENTILES:
        data[f"{p}%"] = []

    with open(os.path.join(MODEL_OUTPUT_DIR, "fits_summary.tsv")) as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            if row["location"] != "Overall":
                continue
            if (row["study"]) == "brinch":
                continue
            data["predictor_type"].append(row["predictor_type"])
            data["virus"].append(row["tidy_name"])
            data["study"].append(row["study"])
            data["location"].append(row["location"])
            data["mean"].append(float(row["mean"]))
            data["std"].append(float(row["std"]))
            data["min"].append(float(row["min"]))
            data["max"].append(float(row["max"]))
            for p in PERCENTILES:
                data[f"{p}%"].append(log(float(row[f"{p}%"]), 10))

    df = pd.DataFrame.from_dict(data)
    return df


def reads_df() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(MODEL_OUTPUT_DIR, "input.tsv"), sep="\t")
    return df


def study_name(study: str) -> str:
    return {
        "brinch": "Brinch (DNA)",
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
        "spurbeck": "Spurbeck",
    }[study]


def count_viral_reads(
    df: pd.DataFrame, by_location: bool = False
) -> pd.DataFrame:
    groups = [
        "pathogen",
        "study",
        "tidy_name",
    ]
    reads_by_study_and_pathogen = (
        df.groupby(groups)[["viral_reads"]].sum().reset_index()
    )
    # create a tsv
    reads_by_study_and_pathogen.to_csv(
        "reads_by_study_and_pathogen.tsv", sep="\t"
    )
    return reads_by_study_and_pathogen


def compute_diffs(df: pd.DataFrame) -> pd.DataFrame:
    viruses = df["virus"].unique()
    results_data = defaultdict(list)
    for virus in viruses:
        virus_df = df[df["virus"] == virus]
        if virus_df["study"].nunique() < 2:
            continue
        if (virus_df["viral_reads"] == 0).sum() >= 2:
            continue
        virus_df = virus_df[virus_df["viral_reads"] != 0]

        predictor_type = virus_df["predictor_type"].unique()

        min_median_index = virus_df["50%"].idxmin()
        max_median_index = virus_df["50%"].idxmax()
        diff_median = (
            virus_df.loc[max_median_index, "50%"]
            - virus_df.loc[min_median_index, "50%"]
        )
        low_diff = (
            virus_df.loc[max_median_index, "5%"]
            - virus_df.loc[min_median_index, "95%"]
        )
        high_diff = (
            virus_df.loc[max_median_index, "95%"]
            - virus_df.loc[min_median_index, "5%"]
        )

        results_data["virus"].append(virus)
        results_data["diff_median"].append(diff_median)
        results_data["low_diff"].append(low_diff)
        results_data["high_diff"].append(high_diff)
        results_data["predictor_type"].append(predictor_type)
        results_data["selected_studies"].append(
            [
                virus_df.loc[min_median_index, "study"],
                virus_df.loc[max_median_index, "study"],
            ]
        )
    df = pd.DataFrame.from_dict(results_data)
    return df


def plot_df(df: pd.DataFrame) -> None:
    df = df.sort_values(by="diff_median").reset_index(drop=True)
    df = df.sort_values(by="predictor_type", ascending=False).reset_index(
        drop=True
    )
    fig, ax = plt.subplots(figsize=(6, 4))
    scatter = ax.scatter(
        x=df["diff_median"],
        y=range(len(df)),
        alpha=0.6,
        edgecolors="w",
    )
    for i in range(len(df)):
        ax.plot(
            [df["low_diff"][i], df["high_diff"][i]],
            [i, i],
            color="k",
        )

    x_min, x_max = ax.get_xlim()

    for i in range(round(x_min), round(x_max)):
        ax.axvline(i, color="k", alpha=0.1, lw=0.5)

    for i in range(len(df)):
        min_study, max_study = df["selected_studies"][i]

        study_combo = f"{study_name(max_study)} <-> {study_name(min_study)}"
        ax.text(
            x_min - 0.18,
            i - 0.25,
            study_combo,
            ha="right",
            va="center",
            fontsize=8,
        )

        if i == (len(df)) - 1:
            ax.text(
                x_max + 0.05,
                i,
                "Incidence\nViruses",
                ha="left",
                va="center",
            )

        if i == len(df) - 1:
            break
        if df["predictor_type"][i] != df["predictor_type"][i + 1]:
            ax.axhline(i + 0.5, color="black", alpha=0.3, linestyle="--")
            ax.text(
                x_max + 0.05,
                i + 0.5,
                "Prevalence\nViruses",
                ha="left",
                va="top",
            )

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["virus"])

    ax.set_title("RA(1%) Difference (median, 90% CI)")
    ax.set_xlabel(
        "OOM difference between highest and lowest study estimate (based on median location)"
    )
    return fig


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(
            figdir / f"{name}.{ext}",
            bbox_inches="tight",
            dpi=600,
        )


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    figdir.mkdir(exist_ok=True)
    reads_data = reads_df()
    fits_data = fits_df()

    viral_counts = count_viral_reads(reads_data)

    fits_data_w_reads = pd.merge(
        fits_data,
        viral_counts,
        how="left",
        left_on=["virus", "study"],
        right_on=["tidy_name", "study"],
    )
    diffs_df = compute_diffs(fits_data_w_reads)
    fig = plot_df(diffs_df)
    save_plot(fig, figdir, "fig_s2")


if __name__ == "__main__":
    start()
