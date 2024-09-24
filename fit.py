#!/usr/bin/env python3
from pathlib import Path
import os

import pandas as pd
import stats
from mgs import Enrichment, MGSData, target_bioprojects
from pathogens import predictors_by_taxid

MODEL_OUTPUT_DIR = "model_output"


def summarize_output(coeffs: pd.DataFrame) -> pd.DataFrame:
    return coeffs.groupby(
        [
            "pathogen",
            "tidy_name",
            "taxids",
            "predictor_type",
            "study",
            "location",
        ]
    ).ra_at_1in100.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])


def start(num_samples: int, plot: bool) -> None:
    figdir = os.path.join(MODEL_OUTPUT_DIR, "model_fig")
    if plot:
        os.makedirs(figdir, exist_ok=True)
    mgs_data = MGSData.from_repo()
    input_data = []
    output_data = []
    study_pathogen_rhats = {}
    for (
        pathogen_name,
        tidy_name,
        predictor_type,
        taxids,
        predictors,
    ) in predictors_by_taxid():
        taxids_str = "_".join(str(t) for t in taxids)
        for study, bioprojects in target_bioprojects.items():
            enrichment = None if study == "brinch" else Enrichment.VIRAL
            model = stats.build_model(
                mgs_data,
                bioprojects,
                predictors,
                taxids,
                random_seed=sum(taxids),
                enrichment=enrichment,
            )
            if model is None:
                continue
            model.fit_model(num_samples=num_samples)

            rhat = model.get_rhat()
            study_pathogen_rhats[f"{study}, {tidy_name}"] = rhat

            if plot:
                taxid_str = "-".join(str(tid) for tid in taxids)
                model.plot_figures(
                    path=figdir,
                    prefix=f"{pathogen_name}-{taxid_str}-{predictor_type}-{study}",
                )
            metadata = dict(
                pathogen=pathogen_name,
                tidy_name=tidy_name,
                taxids=taxids_str,
                predictor_type=predictor_type,
                study=study,
            )
            input_data.append(model.input_df.assign(**metadata))
            output_data.append(model.get_coefficients().assign(**metadata))

    input = pd.concat(input_data)
    input.to_csv(
        os.path.join(MODEL_OUTPUT_DIR, "input.tsv"), sep="\t", index=False
    )
    coeffs = pd.concat(output_data)
    coeffs.to_csv(
        os.path.join(MODEL_OUTPUT_DIR, "fits.tsv"), sep="\t", index=False
    )
    summary = summarize_output(coeffs)
    summary.to_csv(
        os.path.join(MODEL_OUTPUT_DIR, "fits_summary.tsv"), sep="\t"
    )
    print("Model fitting for non-enriched samples complete\nR-hat statistics:")
    for pathogen_and_study, rhat in study_pathogen_rhats.items():
        print(f"{pathogen_and_study}: rhat={rhat}")


if __name__ == "__main__":
    # TODO: Command line arguments
    start(num_samples=8000, plot=True)
