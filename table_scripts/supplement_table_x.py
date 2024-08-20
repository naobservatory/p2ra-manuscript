import matplotlib.pyplot as plt
import numpy as np
import csv
import os
from dataclasses import dataclass

MODEL_OUTPUT_DIR = "model_output"
TABLE_OUTPUT_DIR = "tables"

CUM_INC_1_PERC = 0.01
CUM_INC_001_PERC = 0.0001
DETECTION_THRESHOLD = 100
PERCENTILES = [5, 25, 50, 75, 95]
NOVASEQ_LANE_COST = 3492
NOVASEQ_LANE_DEPTH = 1.4e9
NOVASEQ_CELL_COST = 26462
NOVASEQ_CELL_DEPTH = 23.4e9
MISEQ_COST = 1104
MISEQ_DEPTH = 2e6
NEXTSEQ_COST = 1956
NEXTSEQ_DEPTH = 45e6

ENRICHMENT_COST = 324
LIBRARY_PREP_COST = 302
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

    # predictor_type = "incidence"

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


def get_cost_table():
    with open(
        os.path.join("..", TABLE_OUTPUT_DIR, "cost_table.tsv"), "w"
    ) as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(
            (
                "Study",
                "Virus",
                "Enrichment",
                "Total cost (per sample)",
                "Total cost (per year)",
                "Sequencing depth",
                "Sequencer",
                "Cumulative incidence",
            )
        )
        for cumulative_incidence in [CUM_INC_1_PERC, CUM_INC_001_PERC]:
            for enriched in [True, False]:
                for study in study_labels:
                    for virus in ["SARS-COV-2", "Norovirus (GII)"]:

                        # No panel-enriched estimates for Spurbeck
                        if study == "spurbeck" and enriched:
                            continue

                        # Computing sequencing cost
                        total_cost = 0
                        seq_cost, sequencer, seq_depth = get_sequencing_cost(
                            virus, cumulative_incidence, study, enriched
                        )
                        if enriched:
                            total_cost += ENRICHMENT_COST
                        if study == "crits_christoph" and not enriched:
                            total_cost += RIBODEPLETION_COST
                        total_cost += LIBRARY_PREP_COST
                        total_cost += seq_cost
                        yearly_cost = total_cost * 52

                        # Prettifying everything
                        total_cost = f"${round(int(total_cost), -2):,}"
                        yearly_cost = f"${round(int(yearly_cost), -2):,}"
                        seq_depth = f"{seq_depth:.2e}"
                        enrichment_label = (
                            "panel-enriched" if enriched else "unenriched"
                        )
                        pretty_study = study_labels[study]
                        pretty_cum_inc = f"{cumulative_incidence:.2%}"
                        writer.writerow(
                            (
                                pretty_study,
                                virus,
                                enrichment_label,
                                total_cost,
                                yearly_cost,
                                seq_depth,
                                sequencer,
                                pretty_cum_inc,
                            )
                        )


get_cost_table()
