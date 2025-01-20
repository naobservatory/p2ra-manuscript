#! /usr/bin/env python3
import pandas as pd
from collections import defaultdict
import os


def read_csv_file(file_path):
    return pd.read_csv(file_path, sep="\t")


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

    exponent = exponent.replace("+", "")
    if exponent.startswith("-"):
        exponent = "⁻" + exponent[2:].lstrip("0")
    else:
        exponent = exponent.lstrip("0")

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

    return f"{coefficient} x 10{exponent}"


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
