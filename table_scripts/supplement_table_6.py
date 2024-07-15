import csv
import pandas as pd
import os

from collections import defaultdict
from math import log
from scipy.stats import gmean

PERCENTILES = ["5%", "25%", "50%", "75%", "95%"]
MODEL_OUTPUT_DIR = "../model_output"
TABLE_OUTPUT_DIR = "../tables"


def reads_df() -> pd.DataFrame:
    df = pd.read_csv("input.tsv", sep="\t")
    return df


def spurbeck_fits_data() -> pd.DataFrame:
    data = {
        "predictor_type": [],
        "virus": [],
        "study": [],
        "location": [],
        "ultrafiltrated": [],
    }
    for p in PERCENTILES:
        data[f"{p}"] = []

    with open(os.path.join(MODEL_OUTPUT_DIR, "fits_summary.tsv")) as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            if row["location"] == "Overall":
                continue

            if row["study"] != "spurbeck":
                continue

            if row["location"] in ["E", "F", "G", "H"]:
                data["ultrafiltrated"].append(True)
            else:
                data["ultrafiltrated"].append(False)

            data["predictor_type"].append(row["predictor_type"])
            data["virus"].append(row["tidy_name"])
            data["study"].append(row["study"])
            data["location"].append(row["location"])
            for p in PERCENTILES:
                data[f"{p}"].append(abs(log(float(row[f"{p}"]), 10)))

    df = pd.DataFrame.from_dict(data)
    return df


def compute_geo_mean_ratio(df: pd.DataFrame) -> pd.DataFrame:
    target_viruses = [
        "Norovirus (GI)",
        "Norovirus (GII)",
        "SARS-COV-2",
        "MCV",
        "JCV",
        "BKV",
    ]
    gmean_variance = defaultdict(list)
    for virus in df["virus"].unique():
        if virus not in target_viruses:
            continue
        virus_df = df[df["virus"] == virus]
        ultrafiltrated_virus_df = virus_df[virus_df["ultrafiltrated"]]
        non_ultrafiltrated_virus_df = virus_df[~virus_df["ultrafiltrated"]]
        gmean_variance["virus"].append(virus)
        for quantile in PERCENTILES:
            ultrafiltrated_gm = gmean(
                ultrafiltrated_virus_df[quantile].dropna()
            )
            non_ultrafiltrated_gm = gmean(
                non_ultrafiltrated_virus_df[quantile].dropna()
            )

            variance = float(ultrafiltrated_gm) - float(non_ultrafiltrated_gm)

            gmean_variance[f"Difference at {quantile}"].append(
                round(variance, 2)
            )

    return pd.DataFrame(gmean_variance)


def start():
    df_fits = spurbeck_fits_data()

    variance_df = compute_geo_mean_ratio(df_fits)

    variance_df.to_csv(
        os.path.join(TABLE_OUTPUT_DIR, "supplement_table_6.tsv"),
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    start()
