#! /usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
from dataclasses import dataclass
from scipy.stats import gmean

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


def round_to_three_digits(num):
    if num == 0:
        return 0
    magnitude = 10 ** (int(np.log10(num)) - 2)
    return round(num / magnitude) * magnitude


def tidy_number(reads_required: float) -> str:
    sci_notation = f"{reads_required:.2e}"

    coefficient, exponent = sci_notation.split("e")

    # Round the coefficient to two decimal places
    coefficient = f"{float(coefficient):.2f}"

    # Remove leading '+' and zeros from exponent
    exponent = exponent.lstrip("+0")

    # Convert exponent to superscript
    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    exponent = exponent.translate(superscript_map)

    return f"{coefficient} × 10{exponent}"


def get_cost_table():
    with open(
        os.path.join("..", TABLE_OUTPUT_DIR, "table_s10.tsv"), "w"
    ) as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(
            (
                "Virus",
                "Enrichment",
                "Cumulative incidence",
                "Study",
                "Processing cost (per sample)",
                "Sequencing cost (per sample)",
                "Total cost (per sample)",
                "Total cost (per year)",
                "Sequencing depth",
                "Sequencer",
            )
        )
        for virus in ["SARS-COV-2", "Norovirus (GII)"]:
            for enriched in [True, False]:
                if enriched and virus == "Norovirus (GII)":
                    continue
                for cumulative_incidence in [CUM_INC_1_PERC, CUM_INC_001_PERC]:
                    seq_costs = []
                    processing_costs = []
                    seq_depths = []
                    total_costs = []
                    yearly_costs = []
                    for study in study_labels:
                        # No panel-enriched estimates for Spurbeck
                        if study == "spurbeck" and enriched:
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

                        seq_costs.append(seq_cost)
                        processing_costs.append(processing_cost)
                        total_costs.append(total_cost)
                        yearly_costs.append(yearly_cost)

                        seq_depths.append(seq_depth)

                        # Prettifying everything
                        total_cost = f"${int(total_cost):,}"
                        yearly_cost = f"${int(yearly_cost):,}"
                        processing_cost = f"${int(processing_cost):,}"
                        seq_cost = f"${int(seq_cost):,}"
                        seq_depth = tidy_number(seq_depth)
                        enrichment_label = (
                            "Panel-enriched" if enriched else "Unenriched"
                        )
                        pretty_study = study_labels[study]
                        pretty_cum_inc = f"{cumulative_incidence:.2%}"

                        writer.writerow(
                            (
                                virus,
                                enrichment_label,
                                pretty_cum_inc,
                                pretty_study,
                                processing_cost,
                                seq_cost,
                                total_cost,
                                yearly_cost,
                                seq_depth,
                                sequencer,
                            )
                        )
                    mean_seq_cost = f"${int(gmean(seq_costs)):,}"
                    mean_processing_cost = f"${int(gmean(processing_costs)):,}"
                    mean_total_cost = f"${int(gmean(total_costs)):,}"
                    mean_yearly_cost = f"${int(gmean(yearly_costs)):,}"

                    mean_seq_depth = tidy_number(gmean(seq_depths))
                    writer.writerow(
                        (
                            virus,
                            enrichment_label,
                            pretty_cum_inc,
                            "Overall (geomean)",
                            mean_processing_cost,
                            mean_seq_cost,
                            mean_total_cost,
                            mean_yearly_cost,
                            mean_seq_depth,
                            "N/A",
                        )
                    )


get_cost_table()
