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

# Constants for fig_s8
NOVASEQ_LANE_COST = 2494
NOVASEQ_LANE_DEPTH = 1.4e9
NOVASEQ_CELL_COST = 18902
NOVASEQ_CELL_DEPTH = 23.4e9
MISEQ_COST = 788
MISEQ_DEPTH = 2e6
NEXTSEQ_COST = 1397
NEXTSEQ_DEPTH = 45e6
CUM_INC_1_PERC = 0.01
CUM_INC_001_PERC = 0.0001
DETECTION_THRESHOLD = 100

study_labels = {
    "spurbeck": "Spurbeck",
    "crits_christoph": "Crits-Christoph",
    "rothman": "Rothman",
}


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


def get_depth_and_costs():
    novaseq_lane_costs = [NOVASEQ_LANE_COST]
    novaseq_lane_depth = [0]
    novaseq_cell_costs = [NOVASEQ_CELL_COST]
    novaseq_cell_depth = [0]
    miseq_costs = [MISEQ_COST]
    miseq_depth = [0]
    nextseq_costs = [NEXTSEQ_COST]
    nextseq_depth = [0]
    for i in range(1, 1000000):
        novaseq_lane_costs.append(i * NOVASEQ_LANE_COST)
        novaseq_lane_depth.append(i * NOVASEQ_LANE_DEPTH)
        novaseq_cell_costs.append(i * NOVASEQ_CELL_COST)
        novaseq_cell_depth.append(i * NOVASEQ_CELL_DEPTH)
        miseq_costs.append(i * MISEQ_COST)
        miseq_depth.append(i * MISEQ_DEPTH)
        nextseq_costs.append(i * NEXTSEQ_COST)
        nextseq_depth.append(i * NEXTSEQ_DEPTH)

    return (
        novaseq_lane_costs,
        novaseq_lane_depth,
        novaseq_cell_costs,
        novaseq_cell_depth,
        miseq_costs,
        miseq_depth,
        nextseq_costs,
        nextseq_depth,
    )


def get_cost(virus, cumulative_incidence):
    seq_depths = {}

    for study in study_labels:
        for enriched in [True, False]:
            if study == "spurbeck" and enriched == True:
                continue

            seq_depths[virus, study, enriched] = get_reads_required(
                data,
                cumulative_incidence,
                DETECTION_THRESHOLD,
                virus,
                "Overall",
                "incidence",
                study,
                enriched,
            )[
                0
            ]  # Only take median value

    seq_costs = {}
    for (virus, study, enriched), depth in seq_depths.items():
        nova_seq_lane_cost = NOVASEQ_LANE_COST * np.ceil(
            depth / NOVASEQ_LANE_DEPTH
        )
        nova_seq_cell_cost = NOVASEQ_CELL_COST * np.ceil(
            depth / NOVASEQ_CELL_DEPTH
        )
        miseq_cost = MISEQ_COST * np.ceil(depth / MISEQ_DEPTH)
        nextseq_cost = NEXTSEQ_COST * np.ceil(depth / NEXTSEQ_DEPTH)
        lowest_cost = min(
            nova_seq_lane_cost,
            nova_seq_cell_cost,
            miseq_cost,
            nextseq_cost,
        )

        cost = lowest_cost
        sequencer = (
            "NovaSeq (lane)"
            if lowest_cost == nova_seq_lane_cost
            else (
                "NovaSeq (cell)"
                if lowest_cost == nova_seq_cell_cost
                else ("MiSeq" if lowest_cost == miseq_cost else "NextSeq")
            )
        )

        seq_costs[virus, study, enriched] = (sequencer, cost, depth)
    return seq_costs


def plot_fig_4(fig, gs):
    axes = [fig.add_subplot(gs[0, i]) for i in range(2)]  # Place in first row
    viruses = ["Norovirus (GII)", "SARS-COV-2"]
    study_colors = {
        "spurbeck": "#1b9e77",
        "crits_christoph": "#d95f02",
        "rothman": "#7570b3",
    }
    fig.text(
        0.05,
        0.91,
        "a ",
        fontsize=12,
        fontweight="bold",
    )

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

                if i == 0:
                    ax.set_title(
                        f"{virus}",
                        loc="left",
                        x=0,
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
                ax.set_ylim(1e1, 1e15)

                ax.grid(
                    which="major",
                    linestyle="-",
                    linewidth=0.5,
                    color="gray",
                    alpha=0.7,
                )

                i += 1

    # Set y-axis label only for the leftmost plot
    axes[0].set_ylabel("Reads per week")

    # # Remove y-ticks from right plot
    # axes[1].tick_params(
    #     axis="y", which="both", left=False, right=False, labelleft=False
    # )

    # Set x-axis label for both plots
    for ax in axes:
        ax.set_xlabel("Cumulative Incidence")
        ax.tick_params(axis="x", which="minor", bottom=False)

    # Add legend under both plots
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles=handles,
        labels=labels,
        bbox_to_anchor=(1.06, 0.68),  # Moved to right side
        loc="lower right",
        ncol=1,
        fontsize=8,
    )


def plot_fig_s8(fig, gs):
    axs = []
    # Create subplots in order: [noro_1%, noro_0.01%, scv2_1%, scv2_0.01%]
    for i in range(2):  # For each virus
        gs_sub = gs[1, i].subgridspec(
            2, 1
        )  # Create sub-gridspec for each column
        for j in range(2):  # For each incidence level
            axs.append(
                fig.add_subplot(gs_sub[j, 0])
            )  # Stack plots vertically in each column

    Y_MIN = 10**2
    Y_MAX = 10**8

    fig.text(
        0.05,
        0.565,
        "b",
        fontsize=12,
        fontweight="bold",
    )


    i = 0
    for virus in ["Norovirus (GII)", "SARS-COV-2"]:
        for cumulative_incidence in [CUM_INC_1_PERC, CUM_INC_001_PERC]:
            ax = axs[i]

            (
                novaseq_lane_costs,
                novaseq_lane_depth,
                novaseq_cell_costs,
                novaseq_cell_depth,
                miseq_costs,
                miseq_depth,
                nextseq_costs,
                nextseq_depth,
            ) = get_depth_and_costs()

            ax.step(
                miseq_depth,
                miseq_costs,
                where="pre",
                color="black",
                linestyle="-.",
                linewidth=1,
            )
            ax.step(
                novaseq_lane_depth,
                novaseq_lane_costs,
                where="pre",
                color="black",
                linestyle="--",
                linewidth=1,
            )
            ax.step(
                nextseq_depth,
                nextseq_costs,
                where="pre",
                color="black",
                linestyle="-",
                linewidth=1,
            )
            ax.step(
                novaseq_cell_depth,
                novaseq_cell_costs,
                where="pre",
                color="black",
                linestyle=":",
                linewidth=1,
            )

            seq_costs = get_cost(virus, cumulative_incidence)

            study_colors = {
                "crits_christoph": "#d95f02",
                "rothman": "#7570b3",
                "spurbeck": "#1b9e77",
            }
            enriched_icons = {
                True: "x",
                False: "o",
            }

            for (virus_key, study, enriched), (
                sequencer,
                cost,
                depth,
            ) in seq_costs.items():
                if study == "spurbeck" and enriched == True:
                    continue
                if virus == "Norovirus (GII)" and enriched == True:
                    continue
                study_label = f"{study_labels[study]}"
                if enriched:
                    study_label += "\n(panel-enriched)"
                else:
                    study_label += "\n(unenriched)"
                ax.scatter(
                    depth,
                    cost,
                    color=study_colors[study],
                    marker=enriched_icons[enriched],
                    label=study_label,
                    zorder=10,
                )

            if i == 0:
                for cost, name in zip(
                    [
                        MISEQ_COST,
                        NEXTSEQ_COST,
                        NOVASEQ_LANE_COST,
                        NOVASEQ_CELL_COST,
                    ],
                    ["MiSeq", "NextSeq", "NovaSeq (lane)", "NovaSeq (cell)"],
                ):
                    ax.text(
                        1e3, cost, name, ha="left", va="bottom", fontsize=7
                    )

            if i == 0:
                ax.set_title(
                    f"{virus}, Cumulative incidence (CI): {cumulative_incidence:.2%}",
                    fontsize=10,
                    loc="left",
                    x=0,
                )
            else:
                ax.set_title(
                    f"{virus}, CI: {cumulative_incidence:.2%}",
                    fontsize=10,
                    loc="left",
                    x=0,
                )

            ax.set_xscale("log")
            ax.set_yscale("log")
            if i == 1 or i == 3:  # Bottom row
                ax.set_xlabel("Required Read Depth per week", fontsize=10)
            if i == 0 or i == 1:  # Left column
                ax.set_ylabel("Sequencing Cost (weekly, $)", fontsize=10)

            for y in np.arange(np.log10(Y_MIN), np.log10(Y_MAX), 0.5):
                ax.axhline(
                    10**y, color="black", alpha=0.2, ls="--", linewidth=0.3
                )

            _, x_max = ax.get_xlim()
            x_tick_positions = []
            for x in np.arange(2, np.ceil(np.log10(x_max)), 2):
                x_tick = 10**x
                ax.axvline(
                    x_tick, color="black", alpha=0.2, ls="--", linewidth=0.3
                )
                x_tick_positions.append(x_tick)

            ax.set_xticks(x_tick_positions)
            ax.tick_params(axis="both", which="major", labelsize=8)
            ax.tick_params(axis="both", which="minor", length=0)
            ax.set_xlim(1e2, 1e14)
            ax.set_ylim(Y_MIN, Y_MAX)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            i += 1

    # Collect all unique handles and labels from all subplots
    all_handles = []
    all_labels = []
    seen_labels = set()
    for ax in axs:
        handles, labels = ax.get_legend_handles_labels()
        for h, l in zip(handles, labels):
            if l not in seen_labels:
                seen_labels.add(l)
                all_handles.append(h)
                all_labels.append(l)

    # Add combined legend centered under the right plots
    fig.legend(
        all_handles,
        all_labels,
        loc="lower right",
        bbox_to_anchor=(1.06, 0.27),  # Moved to right side
        ncol=1,
        fontsize=8,
    )


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    os.makedirs(figdir, exist_ok=True)

    global data
    data = read_data()

    # Create figure with 2x3 layout
    fig = plt.figure(
        figsize=(10, 12)
    )  # Increased width to accommodate right legends
    gs = plt.GridSpec(
        2,
        2,
        width_ratios=[1, 1],
        height_ratios=[2, 3.5],
        hspace=0.3,
        wspace=0.15,
    )

    # Plot both figures
    plot_fig_4(fig, gs)
    plot_fig_s8(fig, gs)

    fig.tight_layout()
    save_plot(fig, figdir, "fig_4")


if __name__ == "__main__":
    start()
