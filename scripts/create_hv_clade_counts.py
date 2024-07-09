#!/usr/bin/env python3

import os
import subprocess

if not os.path.exists("../bioprojects"):
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
    "Rothman 2021": ["PRJNA729801"],  # not yet run through the pipeline
}


for study, bioprojects in TARGET_STUDY_METADATA.items():
    for bioproject in bioprojects:
        study_author = study.split()[0]

        if not os.path.exists(
            f"../bioprojects/{study_author}-{bioproject}/hv_hits_putative_filtered.tsv"
        ):
            subprocess.run(
                [
                    "aws",
                    "s3",
                    "cp",
                    f"s3://nao-mgs-wb/{study_author}-{bioproject}/output/hv_hits_putative_filtered.tsv.gz",
                    f"../bioprojects/{study_author}-{bioproject}/hv_hits_putative_filtered.tsv.gz",
                ]
            )
            subprocess.run(
                [
                    "gunzip",
                    f"../bioprojects/{study_author}-{bioproject}/hv_hits_putative_filtered.tsv.gz",
                ]
            )
        if not os.path.exists(
            f"../bioprojects/{study_author}-{bioproject}/hv_hits_putative_collapsed.tsv"
        ):
            subprocess.run(
                [
                    "Rscript",
                    "collapseHV.R",
                    f"../bioprojects/{study_author}-{bioproject}/hv_hits_putative_filtered.tsv",
                    f"../bioprojects/{study_author}-{bioproject}/hv_hits_putative_collapsed.tsv",
                ]
            )
        if not os.path.exists(f"../taxonomy/viral-taxids.tsv"):
            subprocess.run(
                [
                    "aws",
                    "s3",
                    "cp",
                    "s3://nao-mgs-wb/ref/results/viral-taxids.tsv.gz",
                    "../taxonomy/viral-taxids.tsv.gz",
                ]
            )
            subprocess.run(
                [
                    "gunzip",
                    "../taxonomy/viral-taxids.tsv.gz",
                ]
            )

        if not os.path.exists(
            f"../bioprojects/{study_author}-{bioproject}/hv_clade_counts_new.tsv"
        ):
            subprocess.run(
                [
                    "Rscript",
                    "count-viral-taxa.R",
                    "--reads",
                    f"../bioprojects/{study_author}-{bioproject}/hv_hits_putative_collapsed.tsv",
                    "--taxa",
                    "../taxonomy/viral-taxids.tsv",
                    "--output",
                    f"../bioprojects/{study_author}-{bioproject}/hv_clade_counts_new.tsv",
                ]
            )
