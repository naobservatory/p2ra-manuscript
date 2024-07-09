import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

BIOPROJECT_DIR = "bioprojects"

if not os.path.exists(f"../bioprojects"):
    os.mkdir("../bioprojects")


TARGET_STUDY_METADATA = {
    # "Bengtsson-Palme 2016": ["PRJEB14051"],
    # "Munk 2022": [
    #     "PRJEB13831",
    #     "PRJEB27054",
    #     "PRJEB27621",
    #     "PRJEB40798",
    #     "PRJEB40815",
    #     "PRJEB40816",
    #     "PRJEB51229",
    # ],
    # "Brinch 2020": ["PRJEB13832", "PRJEB34633"],
    # "Ng 2019": ["PRJNA438174"],
    # "Maritz 2019": ["PRJEB28033"],
    # "Brumfield 2022": ["PRJNA812772"],
    # "Yang 2020": ["PRJNA645711"],
    # "Spurbeck 2023": ["PRJNA924011"],
    # "CC 2021": ["PRJNA661613"],
    "Rothman 2021": ["PRJNA729801"],
}

sample_files = [
    "hv_clade_counts_new",
    "hv_clade_counts",
]


def read_clade_counts(file_path):
    df = pd.read_csv(file_path, sep="\t")
    sub_df = df[(df["taxid"] == 10239)]["n_reads_clade", "sample"]
    dict = {}
    for reads, sample in sub_df.itertuples():
        dict[sample] = reads
    return dict


def collect_data():
    data = defaultdict(lambda: defaultdict(dict))
    for study, bioprojects in TARGET_STUDY_METADATA.items():
        for bioproject in bioprojects:
            study_author = study.split()[0]
            bioproject_dir = f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}"

            new_file = os.path.join(bioproject_dir, "hv_clade_counts_new.tsv")
            old_file = os.path.join(bioproject_dir, "hv_clade_counts.tsv")
            with open(new_file, "r") as f:
                next(f)
                for line in f:
                    (
                        taxid,
                        name,
                        rank,
                        parent_taxid,
                        sample,
                        n_reads_direct,
                        n_reads_clade,
                    ) = line.strip().split("\t")
                    # if taxid == "10239":  # Check if taxid is 10239 (Viruses)
                    data[sample][name]["old"] = int(n_reads_clade)
            with open(old_file, "r") as f:
                next(f)
                for line in f:
                    (
                        taxid,
                        name,
                        rank,
                        parent_taxid,
                        sample,
                        n_reads_direct,
                        n_reads_clade,
                    ) = line.strip().split("\t")
                    # if taxid == "10239":  # Check if taxid is 10239 (Viruses)
                    data[sample][name]["new"] = int(n_reads_clade)

    differences = defaultdict(lambda: defaultdict(list))
    for sample in data.keys():
        for name in data[sample].keys():
            new_count = data[sample][name].get("new", 0)
            old_count = data[sample][name].get("old", 0)
            differences[sample][name] = int(new_count - old_count)
    return differences


differences = collect_data()

df = pd.DataFrame(differences).transpose()

df.sort_values(by="Viruses", key=abs, ascending=False, inplace=True)
top_10_df = df.head(10)

for sample in top_10_df.index:
    row = df.loc[sample]
    sorted_row = row.abs().sort_values(ascending=False)

    top_5 = sorted_row.head(15)

    print(f"\nSample: {sample}")
    for name in top_5.index:
        value = row[name]
        print(f"{name}: {value}")
