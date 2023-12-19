import os
import json
import gzip
import subprocess
from scipy.stats import gmean

if os.path.basename(os.getcwd()) != "preprint_stats":
    raise RuntimeError("Run this script from preprint_stats/")


dashboard = os.path.expanduser("~/code/mgs-pipeline/dashboard/")

with open(os.path.join(dashboard, "human_virus_sample_counts.json")) as inf:
    human_virus_sample_counts = json.load(inf)

with open(os.path.join(dashboard, "metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.join(dashboard, "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)

with open(os.path.join(dashboard, "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

studies = list(metadata_papers.keys())


def create_table():
    study_abundances = {}

    virus_taxid = 10239

    for study in studies:
        if study not in [
            "Bengtsson-Palme 2016",
            "Munk 2022",
            "Brinch 2020",
            "Ng 2019",
            "Maritz 2019",
            "Brumfield 2022",
            "Rothman 2021",
            "Yang 2020",
            "Spurbeck 2023",
            "Crits-Christoph 2021",
        ]:
            continue

        for bioproject in metadata_papers[study]["projects"]:
            samples = metadata_bioprojects[bioproject]

            if study == "Bengtsson-Palme 2016":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["fine_location"].startswith(
                        "Inlet"
                    )
                ]

            if study == "Ng 2019":
                samples = [
                    sample
                    for sample in samples
                    if metadata_samples[sample]["fine_location"] == "Influent"
                ]

            for sample in samples:
                if metadata_samples[sample].get("enrichment") == "panel":
                    continue

                cladecounts = "%s.tsv.gz" % sample
                if not os.path.exists(f"../cladecounts/{cladecounts}"):
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "s3://nao-mgs/%s/cladecounts/%s"
                            % (bioproject, cladecounts),
                            "../cladecounts/",
                        ]
                    )
                with gzip.open(
                    f"../cladecounts/{cladecounts}", mode="rt"
                ) as inf:
                    for line in inf:
                        (
                            line_taxid,
                            _,
                            _,
                            clade_assignments,
                            _,
                        ) = line.strip().split()
                        taxid = int(line_taxid)
                        clade_hits = int(clade_assignments)
                        if taxid == virus_taxid:
                            if study not in study_abundances:
                                study_abundances[study] = [
                                    [],
                                    [],
                                ]
                            study_abundances[study][1].append(
                                clade_hits / metadata_samples[sample]["reads"]
                            )
                humanreads = "%s.humanviruses.tsv" % sample

                if not os.path.exists(f"../humanviruses/{humanreads}"):
                    print(
                        "Downloading %s from %s" % (humanreads, bioproject),
                        flush=True,
                    )
                    subprocess.check_call(
                        [
                            "aws",
                            "s3",
                            "cp",
                            "s3://nao-mgs/%s/humanviruses/%s"
                            % (bioproject, humanreads),
                            "../humanviruses/",
                        ]
                    )

                with open(f"../humanviruses/{humanreads}") as inf:
                    human_virus_reads = 0
                    for line in inf:
                        (
                            line_taxid,
                            clade_assignments,
                            _,
                        ) = line.strip().split("\t")
                        clade_hits = int(clade_assignments)
                        human_virus_reads += int(clade_hits)
                    human_virus_relative_abundance = (
                        human_virus_reads / metadata_samples[sample]["reads"]
                    )

                    if human_virus_relative_abundance > 0.0:
                        study_abundances[study][0].append(
                            human_virus_relative_abundance
                        )

    gmean_human_virus_shares = []
    gmean_all_virus_shares = []
    table = [
        (
            "Study",
            "Human Viruses",
            "All Viruses",
        )
    ]

    for study, abundances in study_abundances.items():
        gmean_human_virus_shares.append(gmean(abundances[0]))
        gmean_all_virus_shares.append(gmean(abundances[1]))
        table.append(
            (
                study,
                gmean(abundances[0]),
                gmean(abundances[1]),
            )
        )

    table.append(
        (
            "All studies",
            gmean(gmean_human_virus_shares),
            gmean(gmean_all_virus_shares),
        )
    )

    return table


def format_number(num):
    ROUNDING_DIGITS = 2
    if "e" in "{:.9e}".format(num):
        base, exponent = "{:.9e}".format(num).split("e")
        rounded_base = round(float(base), ROUNDING_DIGITS)
        target_read_per_n = int(1 / num)
        return "{} * 10^{} (1 in {})".format(
            rounded_base, int(exponent), target_read_per_n
        )
    else:
        return str(num)


def start():
    table = create_table()
    formatted_table = []
    for row in table:
        formatted_row = []
        for cell in row:
            if isinstance(cell, float):
                formatted_row.append(format_number(cell))
            else:
                formatted_row.append(str(cell))
        formatted_table.append("\t".join(formatted_row))

    with open("supplement_table_1.tsv", "w") as file:
        for line in formatted_table:
            file.write(line + "\n")


if __name__ == "__main__":
    start()
