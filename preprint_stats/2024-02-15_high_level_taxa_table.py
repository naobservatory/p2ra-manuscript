import gzip
import json
import os
import subprocess

from scipy.stats import gmean
from collections import defaultdict
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
    study_abundances = defaultdict(list)

    
    target_taxa = {
        2759: ("eukaryota", "Eukaryota"), 
        2: ("bacteria", "Bacteria"),    
        2157: ("archaea", "Archaea"),
        10239: ("viruses", "Viruses"),
        0: ("Unassigned", "Unassigned"),
        }



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

        if study not in study_abundances: 
            study_abundances[study] = defaultdict(list)

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
                with gzip.open(f"../cladecounts/{cladecounts}") as inf:
                    taxa_abundances = defaultdict(float)
                    for line in inf:
                        (
                            line_taxid,
                            _,
                            _,
                            clade_assignments,
                            _,
                        ) = line.strip().split()
                        taxid = int(line_taxid)
                        clade_assignments = int(clade_assignments)
                        if taxid in target_taxa:
                            relative_abundance = (
                                clade_assignments / metadata_samples[sample]["reads"]
                            )
                            if relative_abundance == 0:
                                continue
                            taxon_name = target_taxa[taxid][1] 
                            study_abundances[study][taxon_name].append(
                                relative_abundance
                            )
                        
    printed_columns = False 
    for study in study_abundances:
        if not printed_columns:
            printed_columns = True
            print("Study", end="\t")
            for taxon in study_abundances[study].keys():
                print(taxon, end="\t")
            print("\n")  
        print(study, end="\t") 
        abundance_in_study = [] 
        for taxon, abundances in study_abundances[study].items():
            abundance = gmean(abundances)
            print(abundance, end="\t")
        print("\n")



def start():
    table = create_table()
if __name__ == "__main__":
    start()
