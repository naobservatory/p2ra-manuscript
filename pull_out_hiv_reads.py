import pandas as pd

TARGET_STUDY_METADATA = {
    "Brinch 2020": ["PRJEB13832", "PRJEB34633"],
    "Spurbeck 2023": ["PRJNA924011"],
    "CC 2021": ["PRJNA661613"],
    "Rothman 2021": ["PRJNA729801"],
}


with open("p2ra_hiv_reads.tsv", "w") as f:
    f.write(
        "study_name\tbio_project_id\tseq_id\ttaxid\tread_sequence\tquery_seq_fwd\tquery_seq_rev\n"
    )
    for study_name, study_metadata in TARGET_STUDY_METADATA.items():

        study_author = study_name.split()[0]
        for bio_project_id in study_metadata:
            putative_reads = pd.read_csv(
                f"bioprojects/{study_author}-{bio_project_id}/hv_hits_putative_filtered.tsv",
                delimiter="\t",
            )
            hiv_reads = putative_reads[putative_reads["taxid"] == 11676]  # [
            #     [
            #         "seq_id",
            #         "taxid",
            #         "query_seq_fwd",
            #         "query_seq_rev",
            #     ]
            # ]

            for _, row in hiv_reads.iterrows():
                f.write(f"{study_name}\t{bio_project_id}\t")
                f.write("\t".join(str(value) for value in row) + "\n")
