import csv
from dataclasses import dataclass
from collections import defaultdict, namedtuple
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import gmean

MODEL_OUTPUT_DIR = "../model_output"

TABLE_OUTPUT_DIR = "../tables"

additional_taxids = {
    39733: "Astroviridae"
    # to be extended
}


BIOPROJECT_DIR = "bioprojects"

if not os.path.exists(f"../bioprojects"):
    os.mkdir("../bioprojects")


TARGET_STUDY_METADATA = {
    "Bengtsson-Palme 2016": ["PRJEB14051"],
    "Munk 2022": [
        "PRJEB13831",
        "PRJEB27054",
        "PRJEB27621",
        "PRJEB40798",
        "PRJEB40815",
        "PRJEB40816",
        "PRJEB51229",
    ],
    "Brinch 2020": ["PRJEB13832", "PRJEB34633"],
    "Ng 2019": ["PRJNA438174"],
    "Maritz 2019": ["PRJEB28033"],
    "Brumfield 2022": ["PRJNA812772"],
    "Yang 2020": ["PRJNA645711"],
    "Spurbeck 2023": ["PRJNA924011"],
    "CC 2021": ["PRJNA661613"],
    "Rothman 2021": ["PRJNA729801"],
}


sample_files = [
    "taxonomic_composition",
    "hv_clade_counts",
    "kraken_reports",
    "qc_basic_stats",
    "sample-metadata",
]


def get_target_viruses():
    prelim_taxids = {}
    with open(os.path.join(MODEL_OUTPUT_DIR, "fits_summary.tsv")) as datafile:
        reader = csv.DictReader(datafile, delimiter="\t")
        for row in reader:
            virus = row["tidy_name"]
            taxid = row["taxids"]
            prelim_taxids[taxid] = virus
    for taxid, virus in additional_taxids.items():
        prelim_taxids[taxid] = virus
    return prelim_taxids


def parse_hv_clade_counts():
    taxids = get_target_viruses()
    table_data = []

    for study, bioprojects in TARGET_STUDY_METADATA.items():
        study_author = study.split()[0]
        study_rel_abuns = defaultdict(lambda: defaultdict(list))
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
            # print(metadata_samples)
            # fine_counts = pd.read_csv(
            #     f"../{BIOPROJECT_DIR}/{study_bioproject}/kraken_reports.tsv",
            #     sep="\t",
            # )

            # fine_counts_dfs = {
            #     sample: df for sample, df in fine_counts.groupby("sample")
            # }

            hv_clade_counts = (
                pd.read_csv(
                    f"../{BIOPROJECT_DIR}/{study_bioproject}/hv_clade_counts.tsv",
                    sep="\t",
                )
                .set_index(["taxid", "sample"])["n_reads_clade"]
                .to_dict()
            )
            # print(hv_clade_counts)

            qc_basic_stats = pd.read_csv(
                f"../{BIOPROJECT_DIR}/{study_bioproject}/qc_basic_stats.tsv",
                sep="\t",
            )

            sample_read_pairs = dict(
                zip(qc_basic_stats["sample"], qc_basic_stats["n_read_pairs"])
            )

            samples = metadata_samples.keys()
            modified_study = study
            if study == "Bengtsson-Palme 2016":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"].startswith(
                        "Inlet"
                    )
                ]
                modified_study = "Bengtsson-Palme 2016"

            if study == "Ng 2019":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["sample_type"] == "Influent"
                ]

            for sample in samples:
                if study == "CC 2021":
                    # print(metadata_samples[sample]["enrichment"])
                    if metadata_samples[sample]["enrichment"] == "enriched":
                        modified_study = "Crits-Christoph 2023 Panel-enriched"
                    elif (
                        metadata_samples[sample]["enrichment"] == "unenriched"
                    ):
                        modified_study = "Crits-Christoph 2023 Unenriched"

                if study == "Rothman 2021":
                    # print(metadata_samples[sample]["enrichment"])
                    if metadata_samples[sample]["enrichment"] == "0":
                        modified_study = "Rothman 2021 Unenriched"
                    elif metadata_samples[sample]["enrichment"] == "1":
                        modified_study = "Rothman 2021 Panel-enriched"

                for taxid in taxids.keys():
                    rel_abun = (
                        hv_clade_counts.get((int(taxid), str(sample)), 0)
                        / sample_read_pairs[sample]
                    )
                    study_rel_abuns[modified_study][taxid].append(rel_abun)

        for study, taxid_rel_abuns in study_rel_abuns.items():
            for taxid, rel_abuns in taxid_rel_abuns.items():
                rel_abuns = np.array(rel_abuns)
                rel_abun_gmean = gmean(rel_abuns[rel_abuns > 0])
                pos_share = (np.sum(rel_abuns > 0) / len(rel_abuns)).round(2)
                num_pos_num_total = f"{np.sum(rel_abuns > 0)}/{len(rel_abuns)}"
                table_data.append(
                    [
                        study,
                        taxids[taxid],
                        taxid,
                        rel_abun_gmean,
                        pos_share,
                        num_pos_num_total,
                    ]
                )
    from tabulate import tabulate

    table = tabulate(table_data, headers="firstrow")
    with open(
        os.path.join(TABLE_OUTPUT_DIR, "supplement_table_x1.tsv"),
        "w",
        newline="",
    ) as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(
            [
                "Study",
                "Taxon",
                "TaxID",
                "Geometric Mean Relative Abundance",
                "Positive Share",
                "Number of Positive/Total",
            ]
        )
        for row in table_data:
            writer.writerow(row)
    return rel_abuns


# def create_tsv():
#     data = read_data()
#     viruses = set()
#     for entry in data.keys():
#         virus, predictor_type = entry[:2]
#         viruses.add((virus, predictor_type))

#     sorted_viruses = sorted(viruses, key=lambda x: (x[1], x[0]))

#     study_tidy = {
#         "rothman": "Rothman",
#         "crits_christoph": "Crits-Christoph",
#         "spurbeck": "Spurbeck",
#         "brinch": "Brinch",
#     }

#     headers = ["Virus", "Study", "Median", "Lower", "Upper"]

#     with open(
#         os.path.join(TABLE_OUTPUT_DIR, "supplement_table_5.tsv"),
#         "w",
#         newline="",
#     ) as file:
#         writer = csv.DictWriter(file, fieldnames=headers, delimiter="\t")
#         writer.writeheader()

#         for virus, predictor_type in sorted_viruses:
#             studies = ["rothman", "crits_christoph", "spurbeck"] + (
#                 ["brinch"] if predictor_type == "prevalence" else []
#             )
#             for study in studies:
#                 stats = data[virus, predictor_type, study, "Overall"]
#                 writer.writerow(
#                     {
#                         "Virus": virus,
#                         "Study": study_tidy[study],
#                         "Median": stats.percentiles[50],
#                         "Lower": stats.percentiles[5],
#                         "Upper": stats.percentiles[95],
#                     }
#                 )


def start():
    parse_hv_clade_counts()


if __name__ == "__main__":
    start()
