#!/usr/bin/env python3

import csv
import numpy as np
import os
from dataclasses import dataclass
from scipy.stats import gmean

PERCENTILES = [5, 25, 50, 75, 95]
MODEL_OUTPUT_DIR = "../model_output"
TABLE_OUTPUT_DIR = "../tables"


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
        os.path.join(MODEL_OUTPUT_DIR, "panel_fits_summary.tsv")
    ) as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            virus = row["tidy_name"]
            predictor_type = row["predictor_type"]
            study = row["study"]
            location = row["location"]
            data[virus, predictor_type, study, location] = SummaryStats(
                mean=float(row["mean"]),
                std=float(row["std"]),
                min=float(row["min"]),
                percentiles={p: float(row[f"{p}%"]) for p in PERCENTILES},
                max=float(row["max"]),
            )
    return data


def get_reads_required(
    data=dict,
    cumulative_incidence=int,
    detection_threshold=np.ndarray,
    virus=str,
    predictor_type=str,
    study=str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    stats = data[virus, predictor_type, study, "Overall"]

    median_reads = detection_threshold / (
        100 * stats.percentiles[50] * cumulative_incidence
    )
    # Low P2RA gives high total reads required. Hence, low p2ra percentile gives higher bound for reads required.
    upper_bound_reads = detection_threshold / (
        100 * stats.percentiles[5] * cumulative_incidence
    )
    lower_bound_reads = detection_threshold / (
        100 * stats.percentiles[95] * cumulative_incidence
    )

    return median_reads, lower_bound_reads, upper_bound_reads


def tidy_number(reads_required=int) -> str:
    sci_notation = f"{reads_required:.2e}"

    coefficient, exponent = sci_notation.split("e")

    exponent = exponent.replace("+", "")
    if exponent.startswith("0") and len(exponent) > 1:
        exponent = exponent[1:]

    exponent = (
        exponent.replace("0", "⁰")
        .replace("1", "¹")
        .replace("2", "²")
        .replace("3", "³")
        .replace("4", "⁴")
        .replace("5", "⁵")
        .replace("6", "⁶")
        .replace("7", "⁷")
        .replace("8", "⁸")
        .replace("9", "⁹")
    )

    return f"{coefficient} × 10{exponent}"


def start():
    data = read_data()
    TARGET_INCIDENCE = 0.01
    TARGET_THRESHOLDS = [10, 100, 1000]
    viruses = ["Norovirus (GII)", "SARS-COV-2", "Influenza A"]
    study_labels = {
        "rothman": "Rothman Panel-enriched",
        "crits_christoph": "Crits-Christoph Panel-enriched",
    }
    with open(
        os.path.join(TABLE_OUTPUT_DIR, "table_s9.tsv"),
        mode="w",
        newline="",
    ) as file:
        tsv_writer = csv.writer(file, delimiter="\t")
        tsv_writer.writerow(
            [
                "Virus",
                "Study",
                "Median",
                "5th Percentile",
                "95th Percentile",
                "Detection Threshold",
            ]
        )
        for detection_threshold in TARGET_THRESHOLDS:
            for virus in viruses:
                geomean_dict = {
                    "median": [],
                    "lower_bound": [],
                    "upper_bound": [],
                }
                studies = study_labels.keys()
                for i, study in enumerate(studies):
                    (
                        median_reads,
                        lower_bound_reads,
                        upper_bound_reads,
                    ) = get_reads_required(
                        data,
                        cumulative_incidence=TARGET_INCIDENCE,
                        detection_threshold=detection_threshold,
                        virus=virus,
                        predictor_type="incidence",
                        study=study,
                    )

                    geomean_dict["median"].append(median_reads)
                    geomean_dict["lower_bound"].append(lower_bound_reads)
                    geomean_dict["upper_bound"].append(upper_bound_reads)

                    tidy_median = tidy_number(median_reads)
                    tidy_lower = tidy_number(lower_bound_reads)
                    tidy_upper = tidy_number(upper_bound_reads)
                    tsv_writer.writerow(
                        [
                            virus,
                            study_labels[study],
                            tidy_median,
                            tidy_lower,
                            tidy_upper,
                            detection_threshold,
                        ]
                    )
                    if i == len(studies) - 1:
                        geomean_median = gmean(geomean_dict["median"])
                        geomean_lower = gmean(geomean_dict["lower_bound"])
                        geomean_upper = gmean(geomean_dict["upper_bound"])

                        tidy_median = tidy_number(geomean_median)
                        tidy_lower = tidy_number(geomean_lower)
                        tidy_upper = tidy_number(geomean_upper)

                        tsv_writer.writerow(
                            [
                                virus,
                                "Mean (geometric)",
                                tidy_median,
                                tidy_lower,
                                tidy_upper,
                                detection_threshold,
                            ]
                        )


if __name__ == "__main__":
    start()
