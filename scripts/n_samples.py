import gzip
import json
import os
import subprocess
import csv
import pandas as pd
from scipy.stats import gmean

if os.path.basename(os.getcwd()) != "table_scripts":
    raise RuntimeError("Run this script from table_scripts/")

from collections import defaultdict, namedtuple

BIOPROJECT_DIR = "bioprojects"
TABLE_DIR = "tables"
if not os.path.exists(os.path.join("..", TABLE_DIR)):
    os.makedirs(os.path.join("..", TABLE_DIR), exist_ok=True)

DEBUG = None


TARGET_STUDY_METADATA = {
    "Brinch 2020": ["PRJEB13832", "PRJEB34633"],
    "Spurbeck 2023": ["PRJNA924011"],
    "CC 2021": ["PRJNA661613"],
    "Rothman 2021": ["PRJNA729801"],
}


def count_samples():
    n_samples_enriched = 0
    n_samples_unenriched = 0

    for study, bioprojects in TARGET_STUDY_METADATA.items():
        study_author = study.split()[0]
        for bioproject in bioprojects:
            study_bioproject = f"{study_author}-{bioproject}"

            metadata_samples = {}
            with open(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/sample-metadata.csv",
                mode="r",
                encoding="utf-8-sig",
            ) as file:
                reader = csv.DictReader(file)
                for row in reader:
                    sample = row.pop("sample")
                    metadata_samples[sample] = row
            samples = metadata_samples.keys()

            if study == "Bengtsson-Palme 2016":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"].startswith(
                        "Inlet"
                    )
                ]

            if study == "Ng 2019":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"] == "Influent"
                ]

            for sample in samples:
                if study_author == "Brinch":
                    n_samples_unenriched += 1
                    continue
                elif study_author == "Spurbeck":
                    n_samples_unenriched += 1
                elif (
                    metadata_samples[sample].get("enrichment") == "enriched"
                    or metadata_samples[sample].get("enrichment") == "1"
                ):
                    n_samples_enriched += 1
                    continue

                elif (
                    metadata_samples[sample].get("enrichment") == "unenriched"
                    or metadata_samples[sample].get("enrichment") == "0"
                ):
                    n_samples_unenriched += 1
                    continue

    print(
        "There are %s enriched samples across Crits-Christoph and Rothman and %s unenriched samples across all four studies (Brinch, Spurbeck, Crits-Christoph, and Rothman)"
        % (n_samples_enriched, n_samples_unenriched)
    )


if __name__ == "__main__":
    count_samples()
