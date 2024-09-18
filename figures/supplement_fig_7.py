import matplotlib.pyplot as plt
import numpy as np
import csv
import os
from dataclasses import dataclass
from pathlib import Path

MODEL_OUTPUT_DIR = "model_output"

WEEKS_PER_YEAR = 52


CUM_INC_1_PERC = 0.01
CUM_INC_001_PERC = 0.0001
DETECTION_THRESHOLD = 100
PERCENTILES = [5, 25, 50, 75, 95]
NOVASEQ_LANE_COST_PER_YEAR = 2494 * WEEKS_PER_YEAR
NOVASEQ_LANE_DEPTH = 1.4e9
NOVASEQ_CELL_COST_PER_YEAR = 18902 * WEEKS_PER_YEAR
NOVASEQ_CELL_DEPTH = 23.4e9
MISEQ_COST_PER_YEAR = 788 * WEEKS_PER_YEAR
MISEQ_DEPTH = 2e6
NEXTSEQ_COST_PER_YEAR = 1397 * WEEKS_PER_YEAR
NEXTSEQ_DEPTH = 45e6


study_labels = {
    "crits_christoph": "Crits-Christoph",
    "rothman": "Rothman",
    "spurbeck": "Spurbeck",
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
            # predictor_type = row["predictor_type"]
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
            # predictor_type = row["predictor_type"]
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


data = read_data()


def get_reads_required(
    data,
    cumulative_incidence,
    detection_threshold,
    virus,
    study,
    enriched,
) -> float:
    location = "Overall"
    stats = data[virus, study, location, enriched]

    median_reads = detection_threshold / (
        100 * stats.percentiles[50] * cumulative_incidence
    )
    return median_reads


def get_cost(virus, cumulative_incidence):

    # predictor_type = "incidence"
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
                study,
                enriched,
            )

    seq_costs = {}
    for (virus, study, enriched), depth in seq_depths.items():

        nova_seq_lane_cost = NOVASEQ_LANE_COST_PER_YEAR * np.ceil(
            depth / NOVASEQ_LANE_DEPTH
        )
        nova_seq_cell_cost = NOVASEQ_CELL_COST_PER_YEAR * np.ceil(
            depth / NOVASEQ_CELL_DEPTH
        )
        miseq_cost = MISEQ_COST_PER_YEAR * np.ceil(depth / MISEQ_DEPTH)
        nextseq_cost = NEXTSEQ_COST_PER_YEAR * np.ceil(depth / NEXTSEQ_DEPTH)
        lowest_cost = min(
            nova_seq_lane_cost,
            nova_seq_cell_cost,
            miseq_cost,
            nextseq_cost,
        )

        yearly_cost = lowest_cost

        sequencer = (
            "NovaSeq (lane)"
            if lowest_cost == nova_seq_lane_cost
            else (
                "NovaSeq (cell)"
                if lowest_cost == nova_seq_cell_cost
                else ("MiSeq" if lowest_cost == miseq_cost else "NextSeq")
            )
        )

        seq_costs[virus, study, enriched] = (sequencer, yearly_cost, depth)
    return seq_costs


def get_depth_and_costs() -> tuple[
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
]:

    novaseq_lane_costs = [NOVASEQ_LANE_COST_PER_YEAR]
    novaseq_lane_depth = [0]
    novaseq_cell_costs = [NOVASEQ_CELL_COST_PER_YEAR]
    novaseq_cell_depth = [0]
    miseq_costs = [MISEQ_COST_PER_YEAR]
    miseq_depth = [0]
    nextseq_costs = [NEXTSEQ_COST_PER_YEAR]
    nextseq_depth = [0]
    for i in range(1, 1000000):

        novaseq_lane_costs.append(i * NOVASEQ_LANE_COST_PER_YEAR)
        novaseq_lane_depth.append(i * NOVASEQ_LANE_DEPTH)
        novaseq_cell_costs.append(i * NOVASEQ_CELL_COST_PER_YEAR)
        novaseq_cell_depth.append(i * NOVASEQ_CELL_DEPTH)
        miseq_costs.append(i * MISEQ_COST_PER_YEAR)
        miseq_depth.append(i * MISEQ_DEPTH)
        nextseq_costs.append(i * NEXTSEQ_COST_PER_YEAR)
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


def plot_steps_and_dots():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    os.makedirs(figdir, exist_ok=True)

    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)  # , dpi=150)
    fig_numeration = ["a", "b", "c", "d"]
    axs = axs.flatten()
    i = 0
    for virus in ["SARS-COV-2", "Norovirus (GII)"]:
        for cumulative_incidence in [CUM_INC_1_PERC, CUM_INC_001_PERC]:
            ax = axs[i]

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
                "crits_christoph": "red",
                "rothman": "blue",
                "spurbeck": "green",
            }
            enriched_icons = {
                True: "o",
                False: "x",
            }

            for (virus, study, enriched), (
                sequencer,
                cost,
                depth,
            ) in seq_costs.items():
                if study == "spurbeck" and enriched == True:
                    continue
                if virus == "Norovirus (GII)" and enriched == True:
                    continue
                if enriched:
                    enriched_label = "enriched"
                else:
                    enriched_label = "unenriched"
                ax.scatter(
                    depth,
                    cost,
                    color=study_colors[study],
                    marker=enriched_icons[enriched],
                    label=f"{study_labels[study]} ({enriched_label})",
                    zorder=10,
                )

            max_depth = max(
                depth
                for (virus, study, enriched), (
                    sequencer,
                    cost,
                    depth,
                ) in seq_costs.items()
            )
            costs = [
                cost
                for (virus, study, enriched), (
                    sequencer,
                    cost,
                    depth,
                ) in seq_costs.items()
            ]

            max_cost = max(costs)
            min_cost = min(costs)
            for cost, name in zip(
                [
                    MISEQ_COST_PER_YEAR,
                    NEXTSEQ_COST_PER_YEAR,
                    NOVASEQ_LANE_COST_PER_YEAR,
                    NOVASEQ_CELL_COST_PER_YEAR,
                ],
                ["MiSeq", "NextSeq", "NovaSeq (lane)", "NovaSeq (cell)"],
            ):
                if i == 2 and name == "NovaSeq (cell)":
                    continue
                ax.text(1e3, cost, name, ha="left", va="bottom", fontsize=10)

            ax.set_xscale("log")
            ax.set_yscale("log")
            if i in [2, 3]:
                ax.set_xlabel("RequiredRead Depth", fontsize=12)
            if i in [0, 2]:
                ax.set_ylabel("Yearly Sequencing Cost ($)", fontsize=12)
            ax.set_title(
                f"{fig_numeration[i]}) {virus} Cumulative incidence: {cumulative_incidence:.2%}",
                fontsize=14,
                loc="left",
            )

            # ax.set_title("Sequencing Cost vs. Read Depth", fontsize=14)

            # ax.grid(True, which="major", axis="both", ls="-", alpha=0.2)

            for y in np.arange(
                np.floor(np.log10(min_cost)), np.ceil(np.log10(max_cost)), 0.5
            ):
                ax.axhline(
                    10**y, color="black", alpha=0.2, ls="--", linewidth=0.3
                )
            _, x_max = ax.get_xlim()
            x_tick_positions = []
            for x in np.arange(2, np.ceil(np.log10(x_max)), 1):
                x_tick = 10**x
                ax.axvline(
                    x_tick,
                    color="black",
                    alpha=0.2,
                    ls="--",
                    linewidth=0.3,
                )
                x_tick_positions.append(x_tick)

            ax.set_xticks(x_tick_positions)
            ax.tick_params(axis="both", which="major", labelsize=10)
            ax.tick_params(axis="both", which="minor", length=0)
            # Set x-axis ticks at every order of magnitude
            ax.set_xlim(1e2, max_depth * 1.1)

            if max_cost < 1e4:
                ax.set_ylim(10**4, 1e4)
            else:
                ax.set_ylim(10**4, max_cost * 1.1)
            ax.set_xlim(1e2, 1e14)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if i == 0:  # Only add legend to the first subplot
                handles, labels = ax.get_legend_handles_labels()
                fig.legend(
                    handles,
                    labels,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 0),
                    ncol=5,
                    fontsize=10,
                )
            i += 1

    plt.tight_layout()
    plt.savefig(figdir / "supplement_fig_7.png", bbox_inches="tight")


plot_steps_and_dots()
