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

sample_files = [
    "taxonomic_composition",
    "hv_clade_counts",
    # "kraken_reports",
    "qc_basic_stats",
    "sample-metadata",
]


def get_data():
    for study, bioprojects in TARGET_STUDY_METADATA.items():
        for bioproject in bioprojects:
            for fname in sample_files:
                study_author = study.split()[0]
                if fname == "sample-metadata":
                    if not os.path.exists(
                        f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.csv"
                    ):
                        subprocess.run(
                            [
                                "aws",
                                "s3",
                                "cp",
                                f"s3://nao-mgs-wb/{study_author}-{bioproject}/output/{fname}.csv",
                                f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.csv",
                            ]
                        )
                else:

                    if not os.path.exists(
                        f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.tsv"
                    ):
                        subprocess.run(
                            [
                                "aws",
                                "s3",
                                "cp",
                                f"s3://nao-mgs-wb/{study_author}-{bioproject}/output/{fname}.tsv.gz",
                                f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.tsv.gz",
                            ]
                        )
                        subprocess.run(
                            [
                                "gzip",
                                "-d",
                                f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}/{fname}.tsv.gz",
                            ]
                        )


def format_rel_abun(rel_abun):
    ROUNDING_DIGITS = 2
    scientific_rel_abun = "{:.9e}".format(rel_abun)
    if "e" in scientific_rel_abun:
        base, exponent = scientific_rel_abun.split("e")
        rounded_base = round(float(base), ROUNDING_DIGITS)
        try:
            target_read_per_n = int(1 / rel_abun)
        except:
            target_read_per_n = "N/A"
        superscript_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
        exponent_unicode = str(int(exponent)).translate(superscript_map)
        return "{} × 10{} (1 in {})".format(
            rounded_base, exponent_unicode, target_read_per_n
        )


def assemble_table():
    get_data()
    virus_hv_rel_abun = namedtuple(
        "virus_hv_rel_abun", ["hv_rel_abun", "virus_rel_abun"]
    )

    table = []
    gmeans_hv = []
    gmeans_virus = []

    for study, bioprojects in TARGET_STUDY_METADATA.items():
        study_relative_abundance = defaultdict(lambda: defaultdict(float))
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

            hv_clade_counts = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/hv_clade_counts.tsv",
                sep="\t",
            ).set_index(["taxid", "sample"])

            taxonomic_composition = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/taxonomic_composition.tsv",
                sep="\t",
            ).set_index(["sample", "classification"])

            qc_basic_stats = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/qc_basic_stats.tsv",
                sep="\t",
            ).set_index(["sample", "stage"])

            samples = metadata_samples.keys()
            modified_study = study

            for sample in samples:
                if study == "CC 2021":
                    # print(metadata_samples[sample]["enrichment"])
                    if metadata_samples[sample]["enrichment"] == "enriched":
                        modified_study = "Crits-Christoph 2021 Panel-enriched"
                    elif (
                        metadata_samples[sample]["enrichment"] == "unenriched"
                    ):
                        modified_study = "Crits-Christoph 2021 Unenriched"

                if study == "Rothman 2021":
                    # print(metadata_samples[sample]["enrichment"])
                    if metadata_samples[sample]["enrichment"] == "0":
                        modified_study = "Rothman 2021 Unenriched"
                    elif metadata_samples[sample]["enrichment"] == "1":
                        modified_study = "Rothman 2021 Panel-enriched"

                try:
                    total_hv_reads = hv_clade_counts.at[
                        (10239, sample), "n_reads_clade"
                    ]
                except KeyError:
                    if DEBUG:
                        print(
                            f"{sample} has no hv hits, thus no entry in hv_clade_counts. Setting total_hv_reads to 0"
                        )
                    total_hv_reads = 0

                total_reads = qc_basic_stats.at[
                    (sample, "raw_concat"), "n_read_pairs"
                ]
                hv_rel_abun = total_hv_reads / total_reads
                virus_rel_abun = taxonomic_composition.at[
                    (sample, "Viral"), "p_reads"
                ]
                study_relative_abundance[modified_study][sample] = (
                    virus_hv_rel_abun(hv_rel_abun, virus_rel_abun)
                )

        for modified_study in study_relative_abundance.keys():
            gmean_virus = gmean(
                [
                    sample.virus_rel_abun
                    for sample in study_relative_abundance[
                        modified_study
                    ].values()
                    if sample.virus_rel_abun > 0
                ]
            )
            gmean_hv = gmean(
                [
                    sample.hv_rel_abun
                    for sample in study_relative_abundance[
                        modified_study
                    ].values()
                    if sample.hv_rel_abun > 0
                ]
            )

            gmeans_hv.append(gmean_hv)

            gmeans_virus.append(gmean_virus)
            table.append(
                [
                    modified_study,
                    format_rel_abun(gmean_virus),
                    format_rel_abun(gmean_hv),
                ]
            )

    gmean_hv_across_studies = gmean(gmeans_hv)

    gmean_virus_across_studies = gmean(gmeans_virus)

    table.append(
        [
            "All studies",
            format_rel_abun(gmean_virus_across_studies),
            format_rel_abun(gmean_hv_across_studies),
        ]
    )
    return table


if __name__ == "__main__":
    table = assemble_table()
    from tabulate import tabulate

    headers = [
        "Study",
        "Geometric Mean Virus Abundance",
        "Geometric Mean HV Abundance",
    ]
    with open(
        os.path.join("..", TABLE_DIR, "supplement_table_1.tsv"), "w"
    ) as f:
        f.write(tabulate(table, headers, tablefmt="tsv"))
