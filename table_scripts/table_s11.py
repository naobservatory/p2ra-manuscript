#! /usr/bin/env python3
import pandas as pd
from collections import defaultdict
import os


def read_csv_file(file_path):
    return pd.read_csv(file_path, sep="\t")


# To create model runs with different pseudocount values:
# Set the pseudocount value in pathogen_properties.py (default is 0.1)
# Alter the name of the output files in fit.py and fit_panel.py, to make it clear
# which pseudocount value was used for the fit (e.g., 0.01, 0.1, 1.0)
# I.e., name the output files:
#    - For enriched: model_output/panel_fits_summary_{pseudocount}.tsv
#    - For unenriched: model_output/fits_summary_{pseudocount}.tsv
# Then, run the scripts for each pseudocount value:
#    - fit.py for unenriched sequencing
#    - fit_panel.py for enriched sequencing
# This should result in files that are imported below.


fits_enriched_001 = read_csv_file("model_output/panel_fits_summary_0.01.tsv")
fits_unenriched_001 = read_csv_file("model_output/fits_summary_0.01.tsv")
fits_enriched_01 = read_csv_file("model_output/panel_fits_summary_0.1.tsv")
fits_unenriched_01 = read_csv_file("model_output/fits_summary_0.1.tsv")
fits_enriched_1 = read_csv_file("model_output/panel_fits_summary_1.tsv")
fits_unenriched_1 = read_csv_file("model_output/fits_summary_1.tsv")

tables_dir = "tables"

items = defaultdict(lambda: defaultdict(float))

pretty_study_names = {
    "crits_christoph": "Crits-Christoph et al. 2021",
    "rothman": "Rothman et al. 2021",
    "spurbeck": "Spurbeck et al. 2023",
}

for df, pseudocount, enrichment in [
    (fits_unenriched_001, 0.01, "unenriched"),
    (fits_enriched_001, 0.01, "enriched"),
    (fits_unenriched_01, 0.1, "unenriched"),
    (fits_enriched_01, 0.1, "enriched"),
    (fits_unenriched_1, 1, "unenriched"),
    (fits_enriched_1, 1, "enriched"),
]:
    for pathogen in [
        "Influenza A",
        "Influenza B",
        "Norovirus (GI)",
        "Norovirus (GII)",
        "SARS-CoV-2",
    ]:
        try:
            pathogen_df = df[
                (df["tidy_name"] == pathogen) & (df["location"] == "Overall")
            ][["study", "tidy_name", "location", "50%"]]
        except Exception as e:
            print(f"No data for {pathogen} {enrichment} {pseudocount}")
            raise

        pathogen_df.rename(columns={"50%": "median"}, inplace=True)
        for line in pathogen_df.itertuples():
            key = (line.study, line.tidy_name, line.location, enrichment)
            items[key][pseudocount] = line.median


def tidy_number(rel_abund=float) -> str:
    sci_notation = f"{rel_abund:.2e}"
    coefficient, exponent = sci_notation.split("e")

    # Convert to superscript, handling negative sign
    superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    exp_num = exponent.replace("+", "")
    if exp_num.startswith("-"):
        exp_num = "⁻" + exp_num[1:].lstrip("0")
    else:
        exp_num = (
            exp_num.lstrip("0") or "0"
        )  # use "0" if all zeros were stripped

    return f"{coefficient} x 10{exp_num.translate(superscript)}"


lines = []
for key, entry in items.items():
    study, tidy_name, location, enrichment = key
    pretty_study_name = pretty_study_names[study]
    pseudocount_001 = entry[0.01]
    pseudocount_01 = entry[0.1]
    pseudocount_1 = entry[1]

    change_01 = 1 - pseudocount_01 / pseudocount_001
    change_1 = 1 - pseudocount_1 / pseudocount_001
    pseudocount_001 = tidy_number(pseudocount_001)
    pseudocount_01 = tidy_number(pseudocount_01)
    pseudocount_1 = tidy_number(pseudocount_1)
    lines.append(
        [
            tidy_name,
            enrichment,
            pretty_study_name,
            pseudocount_001,
            f"{pseudocount_01} ({change_01:.1%})",
            f"{pseudocount_1} ({change_1:.1%})",
        ]
    )

columns = [
    "Pathogen",
    "Enrichment",
    "Study",
    "P2RA at 0.01",
    "P2RA at 0.1 (% Decrease)",
    "P2RA at 1 (% Decrease)",
]
df = pd.DataFrame(lines, columns=columns)
df.sort_values(by=["Pathogen", "Enrichment", "Study"], inplace=True)

df.to_csv(os.path.join(tables_dir, "table_s11.tsv"), sep="\t", index=False)
