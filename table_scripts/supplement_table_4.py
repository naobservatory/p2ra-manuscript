import csv
from dataclasses import dataclass

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import os
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

    return f"{coefficient} x 10{exponent}"


def read_data() -> dict[tuple[str, str, str, str], SummaryStats]:
    data = {}
    with open(os.path.join(MODEL_OUTPUT_DIR, "fits_summary.tsv")) as datafile:
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


def create_tsv():
    data = read_data()
    viruses = set()
    for entry in data.keys():
        virus, predictor_type = entry[:2]
        viruses.add((virus, predictor_type))

    sorted_viruses = sorted(viruses, key=lambda x: (x[1], x[0]))

    study_tidy = {
        "rothman": "Rothman",
        "crits_christoph": "Crits-Christoph",
        "spurbeck": "Spurbeck",
        "brinch": "Brinch",
    }

    headers = ["Virus", "Study", "Median", "Lower", "Upper"]

    with open(
        os.path.join(TABLE_OUTPUT_DIR, "supplement_table_4.tsv"),
        "w",
        newline="",
    ) as file:
        writer = csv.DictWriter(file, fieldnames=headers, delimiter="\t")
        writer.writeheader()

        for virus, predictor_type in sorted_viruses:
            virus_medians = []
            virus_lowers = []
            virus_uppers = []
            studies = ["rothman", "crits_christoph", "spurbeck"] + (
                ["brinch"] if predictor_type == "prevalence" else []
            )
            for study in studies:
                stats = data[virus, predictor_type, study, "Overall"]
                virus_medians.append(stats.percentiles[50])
                virus_lowers.append(stats.percentiles[5])
                virus_uppers.append(stats.percentiles[95])
                writer.writerow(
                    {
                        "Virus": virus,
                        "Study": study_tidy[study],
                        "Median": tidy_number(stats.percentiles[50]),
                        "Lower": tidy_number(stats.percentiles[5]),
                        "Upper": tidy_number(stats.percentiles[95]),
                    }
                )

            study_median_mean = gmean(virus_medians)
            study_upper_mean = gmean(virus_uppers)
            study_lower_mean = gmean(virus_lowers)
            writer.writerow(
                {
                    "Virus": virus,
                    "Study": "Mean (geometric)",
                    "Median": tidy_number(study_median_mean),
                    "Lower": tidy_number(study_lower_mean),
                    "Upper": tidy_number(study_upper_mean),
                }
            )


def start():
    create_tsv()


if __name__ == "__main__":
    start()
