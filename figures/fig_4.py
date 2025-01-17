#!/usr/bin/env python3

import csv
import os
from dataclasses import dataclass
from pathlib import Path
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from matplotlib.lines import Line2D  # type: ignore
from scipy.stats import gmean

PERCENTILES = [5, 25, 50, 75, 95]

MODEL_OUTPUT_DIR = "model_output"

import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42


@dataclass
class SummaryStats:
    mean: float
    std: float
    min: float
    percentiles: dict[int, float]
    max: float


def read_data() -> dict[tuple[str, str, str, str], SummaryStats]:
    data = {}
    with open(
        os.path.join("..", MODEL_OUTPUT_DIR, "fits_summary.tsv")
    ) as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            virus = row["tidy_name"]
            study = row["study"]
            location = row["location"]
            enriched = False
            data[virus, study, location, enriched] = SummaryStats(
                mean=float(row["mean"]),
                std=float(row["std"]),
                min=float(row["min"]),
                percentiles={p: float(row[f"{p}%"]) for p in PERCENTILES},
                max=float(row["max"]),
            )
    with open(
        os.path.join("..", MODEL_OUTPUT_DIR, "panel_fits_summary.tsv")
    ) as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            virus = row["tidy_name"]
            study = row["study"]
            location = row["location"]
            enriched = True
            data[virus, study, location, enriched] = SummaryStats(
                mean=float(row["mean"]),
                std=float(row["std"]),
                min=float(row["min"]),
                percentiles={p: float(row[f"{p}%"]) for p in PERCENTILES},
                max=float(row["max"]),
            )
    return data


def plot_lines(
    ax: plt.Axes,
    median: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    label: str,
    color: str,
    linestyle: str,
    cumulative_incidence: int,
) -> None:
    ax.loglog(
        cumulative_incidence,
        median,
        color=color,
        label=label,
        linestyle=linestyle,
    )

    ax.fill_between(
        cumulative_incidence,
        lower,
        upper,
        color=color,
        alpha=0.2,
    )


def get_reads_required(
    data=dict,
    cumulative_incidence=int,
    detection_threshold=np.ndarray,
    virus=str,
    location=str,
    predictor_type=str,
    study=str,
    enriched=bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    stats = data[virus, study, location, enriched]

    median_reads = detection_threshold / (
        100 * stats.percentiles[50] * cumulative_incidence
    )
    lower_reads = detection_threshold / (
        100 * stats.percentiles[5] * cumulative_incidence
    )
    upper_reads = detection_threshold / (
        100 * stats.percentiles[95] * cumulative_incidence
    )

    return median_reads, lower_reads, upper_reads


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight", dpi=600)


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    os.makedirs(figdir, exist_ok=True)

    data = read_data()

    viruses = ["Norovirus (GII)", "SARS-COV-2"]
    study_labels = {
        "spurbeck": "Spurbeck",
        "crits_christoph": "Crits-Christoph",
        "rothman": "Rothman",
    }
    DETECTION_THRESHOLD = 100

    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(6, 3))

    study_colors = {
        "spurbeck": "#1b9e77",
        "crits_christoph": "#d95f02",
        "rothman": "#7570b3",
    }
    line_styles = ["-", "--"]

    for virus, ax in zip(viruses, axes):
        geomean_dict = {
            "median": [],
            "lower": [],
            "upper": [],
        }
        studies = study_labels.keys()
        i = 0
        for enriched in [False, True]:
            for study in studies:
                if study == "spurbeck" and enriched:
                    continue
                if virus == "Norovirus (GII)" and enriched:
                    continue
                study_median, study_lower, study_upper = get_reads_required(
                    data,
                    cumulative_incidence=np.logspace(-4, -1, 100),
                    detection_threshold=DETECTION_THRESHOLD,
                    virus=virus,
                    location="Overall",
                    predictor_type="incidence",
                    study=study,
                    enriched=enriched,
                )

                geomean_dict["median"].append(study_median)
                geomean_dict["lower"].append(study_lower)
                geomean_dict["upper"].append(study_upper)

                cumulative_incidence = np.logspace(-4, -1, 100)

                ax.set_title(
                    f"{virus}",
                    loc="center",
                    fontdict={"fontsize": 10},
                )

                study_label = f"{study_labels[study]}"
                if enriched:
                    study_label += "\n(panel-enriched)"
                else:
                    study_label += "\n(unenriched)"

                if enriched:
                    linestyle = "--"
                else:
                    linestyle = "-"
                plot_lines(
                    ax=ax,
                    median=study_median,
                    lower=study_lower,
                    upper=study_upper,
                    label=study_label,
                    linestyle=linestyle,
                    color=study_colors[study],
                    cumulative_incidence=cumulative_incidence,
                )

                ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1])
                ax.set_xticklabels(["0.01%", "0.1%", "1%", "10%"], fontsize=8)
                ax.set_yticks([1e3, 1e6, 1e9, 1e12, 1e15])
                ax.tick_params(axis="y", labelsize=8)
                ax.set_xlim(1e-4, 1e-1)

                ax.grid(
                    which="major",
                    linestyle="-",
                    linewidth=0.5,
                    color="gray",
                    alpha=0.7,
                )
                i += 1

    fig.subplots_adjust(hspace=0.1, wspace=0.1)

    # Set y-axis label only for the leftmost plot
    axes[0].set_ylabel("Reads per week")

    # Remove y-ticks from right plot
    axes[1].tick_params(
        axis="y", which="both", left=False, right=False, labelleft=False
    )

    # Set x-axis label for both plots
    for ax in axes:
        ax.set_xlabel("Cumulative Incidence")
        ax.tick_params(axis="x", which="minor", bottom=False)


    # Add legend at the bottom, only showing handles from the right plot (SARS-COV-2)
    handles, labels = axes[1].get_legend_handles_labels()
    legend = fig.legend(
        handles=handles,
        labels=labels,
        bbox_to_anchor=(0.52, -0.18),
        loc="lower center",
        ncol=3,
        fontsize=8,
    )

    fig.tight_layout()
    fig.show()
    save_plot(fig, figdir, "fig_4")


if __name__ == "__main__":
    start()
