import csv
from dataclasses import dataclass
from collections import defaultdict, namedtuple
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import os

PERCENTILES = [5, 25, 50, 75, 95]

MODEL_OUTPUT_DIR = "../model_output"
TABLE_OUTPUT_DIR = "../tables"

TARGET_STUDIES = ["rothman", "crits_christoph"]


@dataclass(frozen=False)
class SummaryStats:
    mean: float
    std: float
    min: float
    percentiles: dict[int, float]
    max: float
    fold_change: float


def tidy_number(reads_required=int) -> str:
    sci_notation = f"{reads_required:.2e}"

    coefficient, exponent = sci_notation.split("e")

    is_negative = exponent.startswith("-")
    if is_negative:
        exponent = exponent[1:]

    exponent = exponent.lstrip("0")

    if is_negative:
        exponent = "⁻" + exponent

    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    exponent = exponent.translate(superscript_map)

    return f"{coefficient} × 10{exponent}"


studies = [
    "rothman_unenriched",
    "rothman_panel_enriched",
    "crits_christoph_unenriched",
    "crits_christoph_panel_enriched",
]

pretty_study = {
    "rothman_unenriched": "Rothman Unenriched",
    "crits_christoph_unenriched": "Crits-Christoph Unenriched",
    "rothman_panel_enriched": "Rothman Panel-enriched",
    "crits_christoph_panel_enriched": "Crits-Christoph Panel-enriched",
}


def read_data() -> dict[tuple[str, str, str, str, str], SummaryStats]:
    data = {}

    with open(os.path.join(MODEL_OUTPUT_DIR, "fits_summary.tsv")) as datafile:
        unenriched_medians = {}
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            study = row["study"]
            if study == "rothman":
                study_unenriched = "rothman_unenriched"
            elif study == "crits_christoph":
                study_unenriched = "crits_christoph_unenriched"
            else:
                continue
            virus = row["tidy_name"]
            predictor_type = row["predictor_type"]
            location = row["location"]
            unenriched_median = float(row[f"{50}%"])
            unenriched_medians[virus, predictor_type, study, location] = (
                unenriched_median
            )
            data[virus, predictor_type, study_unenriched, location] = (
                SummaryStats(
                    mean=tidy_number(float(row["mean"])),
                    std=tidy_number(float(row["std"])),
                    min=tidy_number(float(row["min"])),
                    percentiles={
                        p: tidy_number(float(row[f"{p}%"]))
                        for p in PERCENTILES
                    },
                    max=tidy_number(float(row["max"])),
                    fold_change="NA",
                )
            )

    with open(
        os.path.join(MODEL_OUTPUT_DIR, "panel_fits_summary.tsv")
    ) as datafile:
        enriched_medians = {}
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            study = row["study"]
            if study == "rothman":
                study_enriched = "rothman_panel_enriched"
                study_unenriched = "rothman_unenriched"
            elif study == "crits_christoph":
                study_enriched = "crits_christoph_panel_enriched"
                study_unenriched = "crits_christoph_unenriched"
            else:
                continue
            virus = row["tidy_name"]
            predictor_type = row["predictor_type"]
            location = row["location"]
            enriched_median = float(row[f"{50}%"])
            unenriched_median = unenriched_medians[
                virus, predictor_type, study, location
            ]
            fold_change = enriched_median / unenriched_median
            fold_change_enriched = (
                round(fold_change, 2)
                if fold_change > 1
                else tidy_number(fold_change)
            )
            fold_change_unenriched = (
                round(1 / fold_change, 2)
                if fold_change < 1
                else tidy_number(1 / fold_change)
            )

            data[virus, predictor_type, study_enriched, location] = (
                SummaryStats(
                    mean=tidy_number(float(row["mean"])),
                    std=tidy_number(float(row["std"])),
                    min=tidy_number(float(row["min"])),
                    percentiles={
                        p: tidy_number(float(row[f"{p}%"]))
                        for p in PERCENTILES
                    },
                    max=tidy_number(float(row["max"])),
                    fold_change=fold_change_enriched,
                )
            )
            data[
                virus, predictor_type, study_unenriched, location
            ].fold_change = fold_change_unenriched

    return data


def create_tsv():
    data = read_data()
    viruses = set()
    for entry in data.keys():
        virus, predictor_type = entry[:2]
        viruses.add((virus, predictor_type))
    sorted_viruses = sorted(viruses, key=lambda x: (x[1], x[0]))

    headers = [
        "Virus",
        "Study",
        "Median",
        "5th Percentile",
        "95th Percentile",
        "Fold change\n(un- vs. enriched, median)",
    ]

    with open(
        os.path.join(TABLE_OUTPUT_DIR, "table_s7.tsv"),
        "w",
        newline="",
    ) as file:
        writer = csv.DictWriter(file, fieldnames=headers, delimiter="\t")
        writer.writeheader()

        for virus, predictor_type in sorted_viruses:

            for study in studies:

                stats = data[virus, predictor_type, study, "Overall"]
                writer.writerow(
                    {
                        "Virus": virus,
                        "Study": pretty_study[study],
                        "Median": stats.percentiles[50],
                        "5th Percentile": stats.percentiles[5],
                        "95th Percentile": stats.percentiles[95],
                        "Fold change\n(un- vs. enriched, median)": stats.fold_change,
                    }
                )


def start():
    create_tsv()


if __name__ == "__main__":
    start()
