# Predicting Virus Relative Abundance in Wastewater
Simon L. Grimm, Jeff T. Kaufman, Daniel Rice, Michael M. McLaren, Charlie Whittaker, William J. Bradshaw & Kevin M. Esvelt

## Abstract
Detecting novel pathogens at an early stage requires robust early warning that is both sensitive and pathogen-agnostic. Wastewater metagenomic sequencing (W-MGS) could enable highly pathogen-agnostic disease monitoring, but its sensitivity and financial feasibility are dependent on the relative abundance of novel pathogen sequences in W-MGS data. Here we collate W-MGS data from a diverse range of studies to characterize the relative abundance of known viruses in wastewater samples. We develop a Bayesian statistical model to integrate a subset of these data with epidemiological estimates of incidence and prevalence, and use it to estimate the expected relative abundance of different viral pathogens for a given prevalence or incidence in the community. Our results reveal pronounced variation between sites and studies, with estimates differing by one to three orders of magnitude for the same pathogen. For example, the expected relative abundance of SARS-CoV-2 at weekly incidence of 1% of the population varied between 10-7 and 10-10. Integrating these estimates with a simple cost model highlights substantial variation in the volume of W-MGS required to detect these pathogens. The mean sequencing cost of identifying 100 reads from a new SARS-CoV-2-like pathogen at 1% cumulative incidence was $2,579,000. A Norovirus-like pathogen which sheds more would require sequencing costs of $18,000. This model, and its parameter estimates, represent an important resource for future investigation into the performance of wastewater MGS, and can be extended to incorporate new wastewater datasets as they become available.

### Repository structure
This work is split across two repositories: This repo contains scripts to collect incidence and prevalence estimates, and to build the model, while determining relative abundances from existing data is in the [INSERT FINAL MGS_PIPELINE REPO}
[mgs-pipeline](https://github.com/naobservatory/mgs-pipeline) repo.

Scripts to create the manuscript figures are found in the `figures/` directory. Scripts to generate epidemiological estimates are in the `pathogens/` directory.

### Working with prevalence data

In python, run `import pathogens` and then iterate over `pathogens.pathogens`.
Each pathogen implements an `estimate_prevalences` method which gives one or
more estimates.

Run `./summarize.py` to get an overview of the data.

### Downloading the data

To download the data for the model, build the first composite figure:

```
$ cd figures/
$ python3 composite_fig_1.py
```

### Statistical model

For an overview of the statistical model see [model.md](model.md).

To fit the model, run `./fit.py`. This will create, under `model_output/`:

* `input.tsv`, a table of the input data to the model
* `fits.tsv`, a table of samples from the posterior distribution of model parameters
* `fits_summary.tsv`, a table listing summary statistics of the posterior distributions of model parameters
* `fig/`, a directory containing a large number plots of posterior distributions and samples from posterior predictive distributions
  (see [model.md](model.md) for details)

Once the model has been fit, run `./plot_summaries.py` to create plots of the posterior distribution of $RA(1\perhundred)$.

* `fig/incidence-violin.{pdf,png}`, posteriors for all incidence viruses
* `fig/prevalence-violin.{pdf,png}`, posteriors for all prevalence viruses
* `fig/by_location_incidence-violin.{pdf,png}`, posteriors separated by location for the most common incidence viruses
* `fig/by_location_prevalence-violin.{pdf,png}`, posteriors separated by location for a subset of prevalence viruses

### Development

#### Testing

Before creating a PR or submitting code, run `./check.sh`.  It will run tests that check your formatting.  This also runs automatically on GitHub when making pull requests, but it's much faster to catch problems locally first.

If `./check.sh` complains about formatting or import sorting, you can fix this automatically with `./check.sh --fix`.

#### Installing pystan

Pystan should be installed along with the other requirements when you run:
```
python -m pip install -r requirements-dev.txt"
```
However, on some non-Linux systems (including M2 Macbooks), one of `pystan`'s dependencies,`httpstan`, may fail to install.
To get around this problem, you can [install httpstan from source](https://httpstan.readthedocs.io/en/latest/installation.html#installation-from-source).
Once it is built and installed, you can then install the requirements file as above.
(Note that you can clone the `httpstan` repo anywhere on your computer.
We recommend doing this outside of the `p2ra` repo directory so that git doesn't try to track it.)
