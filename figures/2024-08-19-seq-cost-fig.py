import matplotlib.pyplot as plt
import numpy as np
import csv
import os
from dataclasses import dataclass

MODEL_OUTPUT_DIR = "model_output"

CUM_INC_1_PERC = 0.01
CUM_INC_001_PERC = 0.0001
DETECTION_THRESHOLD = 100
PERCENTILES = [5, 25, 50, 75, 95]
NOVASEQ_COST = 3492
NOVASEQ_DEPTH = 1.4e9
MISEQ_COST = 1104
MISEQ_DEPTH = 2e6
NEXTSEQ_COST = 1956
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
            if study == "spurbeck":  # and enriched == True:
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

        nova_seq_cost = NOVASEQ_COST * np.ceil(depth / NOVASEQ_DEPTH)
        miseq_cost = MISEQ_COST * np.ceil(depth / MISEQ_DEPTH)
        nextseq_cost = NEXTSEQ_COST * np.ceil(depth / NEXTSEQ_DEPTH)
        lowest_cost = min(nova_seq_cost, miseq_cost, nextseq_cost)

        sequencer = (
            "NovaSeq"
            if lowest_cost == nova_seq_cost
            else "MiSeq" if lowest_cost == miseq_cost else "NextSeq"
        )

        seq_costs[virus, study, enriched] = (sequencer, lowest_cost, depth)
    return seq_costs


def get_depth_and_costs() -> tuple[
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
]:

    novaseq_costs = []
    novaseq_depth = []
    miseq_costs = []
    miseq_depth = []
    nextseq_costs = []
    nextseq_depth = []
    for i in range(0, 1000000):

        if i == 0:
            novaseq_costs.append(NOVASEQ_COST)
            novaseq_depth.append(0)
            miseq_costs.append(MISEQ_COST)
            miseq_depth.append(0)
            nextseq_costs.append(NEXTSEQ_COST)
            nextseq_depth.append(0)

        else:
            novaseq_costs.append(i * NOVASEQ_COST)
            novaseq_depth.append(i * NOVASEQ_DEPTH)
            miseq_costs.append(i * MISEQ_COST)
            miseq_depth.append(i * MISEQ_DEPTH)
            nextseq_costs.append(i * NEXTSEQ_COST)
            nextseq_depth.append(i * NEXTSEQ_DEPTH)

    return (
        novaseq_costs,
        novaseq_depth,
        miseq_costs,
        miseq_depth,
        nextseq_costs,
        nextseq_depth,
    )


(
    novaseq_costs,
    novaseq_depth,
    miseq_costs,
    miseq_depth,
    nextseq_costs,
    nextseq_depth,
) = get_depth_and_costs()


def plot_steps_and_dots():
    fig, axs = plt.subplots(2, 2, figsize=(15, 13))  # , dpi=150)
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
                novaseq_depth,
                novaseq_costs,
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
                if enriched:
                    enriched_label = "enriched"
                else:
                    enriched_label = "unenriched"
                ax.scatter(
                    depth,
                    cost,
                    color=study_colors[study],
                    marker=enriched_icons[enriched],
                    label=f"{study_labels[study]} | {enriched_label}",
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
            max_cost = max(
                cost
                for (virus, study, enriched), (
                    sequencer,
                    cost,
                    depth,
                ) in seq_costs.items()
            )
            if max_cost < 3e4:
                for cost, name in zip(
                    [MISEQ_COST, NEXTSEQ_COST, NOVASEQ_COST],
                    ["MiSeq", "NextSeq", "NovaSeq"],
                ):
                    ax.text(
                        1e4, cost, name, ha="left", va="bottom", fontsize=10
                    )
            else:
                for depth, name, cost in zip(
                    [MISEQ_DEPTH, NEXTSEQ_DEPTH, NOVASEQ_DEPTH],
                    ["MiSeq", "NextSeq", "NovaSeq"],
                    [MISEQ_COST, NEXTSEQ_COST, NOVASEQ_COST],
                ):
                    y_pos = max_cost / 2
                    num_seq_runs = round(y_pos / cost)
                    seq_depth = num_seq_runs * depth
                    ax.text(
                        seq_depth * 1.1,
                        y_pos,
                        name,
                        ha="left",
                        va="bottom",
                        fontsize=10,
                    )

            ax.set_xscale("log")
            ax.set_xlabel("Read Depth", fontsize=12)
            ax.set_ylabel("Cost ($)", fontsize=12)
            ax.set_title(
                f"{fig_numeration[i]}) | {virus} Cummulative incidence: {cumulative_incidence:.2%}",
                fontsize=14,
                loc="left",
            )

            # ax.set_title("Sequencing Cost vs. Read Depth", fontsize=14)

            ax.grid(True, which="major", axis="both", ls="-", alpha=0.2)

            ax.tick_params(axis="both", which="major", labelsize=10)
            ax.tick_params(axis="both", which="minor", length=0)
            ax.set_xlim(1e2, max_depth * 1.1)
            if max_cost < 1e4:
                ax.set_ylim(0, 1e4)
            else:
                ax.set_ylim(0, max_cost * 1.1)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            ax.legend(loc="upper left", fontsize=10)
            i += 1

    plt.tight_layout()
    plt.show()


plot_steps_and_dots()
