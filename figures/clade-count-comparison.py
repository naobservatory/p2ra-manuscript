import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

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
    data = defaultdict()
    for study, bioprojects in TARGET_STUDY_METADATA.items():
        for bioproject in bioprojects:
            study_author = study.split()[0]
            bioproject_dir = f"../{BIOPROJECT_DIR}/{study_author}-{bioproject}"

            new_file = os.path.join(bioproject_dir, "hv_clade_counts_new.tsv")
            old_file = os.path.join(bioproject_dir, "hv_clade_counts.tsv")
            with open(new_file, "r") as f:
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
                    if taxid == "10239":  # Check if taxid is 10239 (Viruses)
                        data[sample] = [int(n_reads_clade)]
            with open(old_file, "r") as f:
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
                    if taxid == "10239":  # Check if taxid is 10239 (Viruses)
                        data[sample].append(int(n_reads_clade))
    new_counts = []
    old_counts = []
    for key in data.keys():
        new_counts.append(data[key][0])
        old_counts.append(data[key][1])

    return new_counts, old_counts


def plot_comparison():
    plt.figure(figsize=(10, 10))
    sns.scatterplot(x=old_counts, y=new_counts)
    plt.xlabel("n_reads_clade (old)")
    plt.ylabel("n_reads_clade (new)")
    plt.title(
        "Comparison of n_reads_clade between old and new hv_clade_counts"
    )

    # Add diagonal line
    # max_val = max(data["old"].max(), data["new"].max())
    # plt.plot([0, max_val], [0, max_val], "r--")

    plt.xscale("log")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig("hv_clade_counts_comparison.png", dpi=300)
    plt.close()


new_counts, old_counts = collect_data()

plot_comparison()
