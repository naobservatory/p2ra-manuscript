import pandas as pd
from collections import defaultdict


def read_csv_file(file_path):
    return pd.read_csv(file_path, sep="\t")


fits_enriched_001 = read_csv_file("model_output/panel_fits_summary_0.01.tsv")
fits_unenriched_001 = read_csv_file("model_output/fits_summary_0.01.tsv")
fits_enriched_01 = read_csv_file("model_output/panel_fits_summary_0.1.tsv")
fits_unenriched_01 = read_csv_file("model_output/fits_summary_0.1.tsv")
fits_enriched_1 = read_csv_file("model_output/panel_fits_summary_1.tsv")
fits_unenriched_1 = read_csv_file("model_output/fits_summary_1.tsv")


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
    flu_a = df[
        (df["tidy_name"] == "Influenza A") & (df["location"] == "Overall")
    ][["study", "tidy_name", "location", "50%"]]

    flu_b = df[
        (df["tidy_name"] == "Influenza B") & (df["location"] == "Overall")
    ][["study", "tidy_name", "location", "50%"]]

    for flu_df in [flu_a, flu_b]:
        flu_df.rename(columns={"50%": "median"}, inplace=True)

        for line in flu_df.itertuples():
            key = (line.study, line.tidy_name, line.location, enrichment)
            items[key][pseudocount] = line.median


with open("sensitivity_table.tsv", "w") as f:
    f.write(
        "Study\tPathogen\tEnrichment\tP2RA at 0.01\tP2RA at 0.1 (% Decrease)\tP2RA at 1 (% Decrease)\n"
    )
    for key, entry in items.items():
        study, tidy_name, location, enrichment = key
        pretty_study_name = pretty_study_names[study]
        pseudocount_001 = entry[0.01]
        pseudocount_01 = entry[0.1]
        pseudocount_1 = entry[1]

        change_01 = 1 - pseudocount_01 / pseudocount_001
        change_1 = 1 - pseudocount_1 / pseudocount_01

        f.write(
            f"{pretty_study_name}\t{tidy_name}\t{enrichment}\t{pseudocount_001:.2e}\t{pseudocount_01:.2e} ({change_01:.1%})\t{pseudocount_1:.2e} ({change_1:.1%})\n"
        )
