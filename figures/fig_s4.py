#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import os
from dataclasses import dataclass
from pathlib import Path
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42

MODEL_OUTPUT_DIR = "model_output"
TABLE_OUTPUT_DIR = "tables"

CUM_INC_1_PERC = 0.01
CUM_INC_001_PERC = 0.0001
DETECTION_THRESHOLD = 100
PERCENTILES = [5, 25, 50, 75, 95]
NOVASEQ_LANE_COST = 2494
NOVASEQ_LANE_DEPTH = 1.4e9
NOVASEQ_CELL_COST = 18902
NOVASEQ_CELL_DEPTH = 23.4e9
MISEQ_COST = 788
MISEQ_DEPTH = 2e6
NEXTSEQ_COST = 1397
NEXTSEQ_DEPTH = 45e6


# https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/respiratory-virus-oligo-panel.html#tabs-0f175ae031-item-7349ca530e-order
# $10,368 for 32 reactions.
ENRICHMENT_COST = 324
# https://bauercore.fas.harvard.edu/fy25-ngs-library-prep#:~:text=1/4%20volume-,%24216,-%24302
# $216 per sample with Harvard Account Code
LIBRARY_PREP_COST = 216
# https://www.illumina.com/products/by-type/molecular-biology-reagents/ribo-zero-plus-rrna-depletion.html#tabs-4d64a43abe-item-354b58d9fa-order
# $910 per sample for 16 samples
RIBODEPLETION_COST = 57

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


def get_sequencing_cost(virus, cumulative_incidence, study, enriched):
    seq_depth = get_reads_required(
        data,
        cumulative_incidence,
        DETECTION_THRESHOLD,
        virus,
        study,
        enriched,
    )

    novaseq_lane_cost = NOVASEQ_LANE_COST * np.ceil(
        seq_depth / NOVASEQ_LANE_DEPTH
    )
    miseq_cost = MISEQ_COST * np.ceil(seq_depth / MISEQ_DEPTH)
    nextseq_cost = NEXTSEQ_COST * np.ceil(seq_depth / NEXTSEQ_DEPTH)
    novaseq_cell_cost = NOVASEQ_CELL_COST * np.ceil(
        seq_depth / NOVASEQ_CELL_DEPTH
    )
    lowest_cost = min(
        novaseq_lane_cost, miseq_cost, nextseq_cost, novaseq_cell_cost
    )

    sequencer = (
        "NovaSeq (lane)"
        if lowest_cost == novaseq_lane_cost
        else (
            "MiSeq"
            if lowest_cost == miseq_cost
            else "NextSeq" if lowest_cost == nextseq_cost else "NovaSeq (cell)"
        )
    )

    seq_cost = lowest_cost
    return seq_cost, sequencer, seq_depth


def get_cost_data():
    virus_study_1_perc = []
    virus_study_001_perc = []
    seq_costs_1_perc = []
    seq_costs_001_perc = []
    processing_costs_1_perc = []
    processing_costs_001_perc = []

    for cumulative_incidence in [CUM_INC_1_PERC, CUM_INC_001_PERC]:
        for enriched in [True, False]:
            for study in study_labels:
                for virus in ["SARS-COV-2", "Norovirus (GII)"]:

                    # No panel-enriched estimates for Spurbeck
                    if study == "spurbeck" and enriched:
                        continue
                    if enriched and virus == "Norovirus (GII)":
                        continue

                    # Computing sequencing cost

                    seq_cost, sequencer, seq_depth = get_sequencing_cost(
                        virus, cumulative_incidence, study, enriched
                    )
                    processing_cost = 0
                    if enriched:
                        processing_cost += ENRICHMENT_COST
                    if study == "crits_christoph" and not enriched:
                        processing_cost += RIBODEPLETION_COST
                    processing_cost += LIBRARY_PREP_COST
                    total_cost = seq_cost + processing_cost
                    yearly_cost = total_cost * 52

                    enrichment_label = (
                        "Panel-enriched" if enriched else "Unenriched"
                    )
                    pretty_study = study_labels[study]

                    study_label = (
                        f"{pretty_study} | {virus} ({enrichment_label}"
                    )

                    if cumulative_incidence == CUM_INC_1_PERC:
                        processing_costs_1_perc.append(processing_cost)
                        seq_costs_1_perc.append(seq_cost)
                        virus_study_1_perc.append(study_label)
                    else:
                        processing_costs_001_perc.append(processing_cost)
                        seq_costs_001_perc.append(seq_cost)
                        virus_study_001_perc.append(study_label)
    costs_1_perc = {
        "Processing Cost": processing_costs_1_perc,
        "Sequencing Cost": seq_costs_1_perc,
    }
    costs_001_perc = {
        "Processing Cost": processing_costs_001_perc,
        "Sequencing Cost": seq_costs_001_perc,
    }
    return (
        virus_study_1_perc,
        virus_study_001_perc,
        costs_1_perc,
        costs_001_perc,
    )


def save_plot(fig, figdir: Path, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(figdir / f"{name}.{ext}", bbox_inches="tight", dpi=600)


def create_fig():

    virus_study_1_perc, virus_study_001_perc, costs_1_perc, costs_001_perc = (
        get_cost_data()
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    left_1_perc = np.zeros(len(virus_study_1_perc))
    left_001_perc = np.zeros(len(virus_study_001_perc))

    for cost_type in costs_1_perc:
        ax1.barh(
            virus_study_1_perc,
            costs_1_perc[cost_type],
            left=left_1_perc,
            label=cost_type,
        )
        left_1_perc += costs_1_perc[cost_type]
        ax1.set_title("1% Cumulative Incidence")

    for cost_type in costs_001_perc:
        ax2.barh(
            virus_study_001_perc,
            costs_001_perc[cost_type],
            left=left_001_perc,
            label=cost_type,
        )
        left_001_perc += costs_001_perc[cost_type]
        ax2.set_title("0.01% Cumulative Incidence")
    ax1.set_xlim(1, 1e4)
    ax1.text(
        1e4 + 500,
        len(virus_study_1_perc),
        "Total Cost per sample",
        va="center",
        ha="left",
        fontsize=12,
    )
    for i, (label, seq_cost, processing_cost) in enumerate(
        zip(
            virus_study_1_perc,
            costs_1_perc["Sequencing Cost"],
            costs_1_perc["Processing Cost"],
        )
    ):

        total_cost = processing_cost + seq_cost
        ax1.text(
            1e4 + 500,
            i,
            f"${total_cost:,.0f}",
            va="center",
            ha="left",
            fontsize=10,
        )

    for i, (label, seq_cost, processing_cost) in enumerate(
        zip(
            virus_study_001_perc,
            costs_001_perc["Sequencing Cost"],
            costs_001_perc["Processing Cost"],
        )
    ):

        total_cost = processing_cost + seq_cost
        ax2.text(
            1e4 + 500,
            i,
            f"${total_cost:,.0f}",
            va="center",
            ha="left",
            fontsize=10,
        )
    # Format x-axis labels with dollar signs

    ax2.xaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: f"${x:,.0f}")
    )

    plt.subplots_adjust(right=0.85)
    ax2.set_xlabel("Processing and sequencing cost per week", fontsize=12)

    ax1.legend()

    plt.tight_layout()
    return fig


def start():
    parent_dir = Path("..")
    figdir = Path(parent_dir / "fig")
    os.makedirs(figdir, exist_ok=True)
    fig = create_fig()
    save_plot(fig, figdir, "fig_s4")


if __name__ == "__main__":
    start()
