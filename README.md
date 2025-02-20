# Predicting Virus Relative Abundance in Wastewater
Simon L. Grimm, Jeff T. Kaufman, Daniel Rice, Michael M. McLaren, Charles Whittaker, William J. Bradshaw & Kevin M. Esvelt

## Abstract
### Background
Metagenomic sequencing of wastewater (W-MGS) can in principle detect any known or novel pathogen in a population. We quantify the sensitivity and cost of W-MGS for viral pathogen detection by jointly analysing W-MGS and epidemiological data for a range of human-infecting viruses.
### Methods
Sequencing data from four studies were analysed to estimate the relative abundance (RA) of 11 human-infecting viruses. Corresponding prevalence and incidence estimates were obtained or calculated from academic and public-health reports. These estimates were combined using a hierarchical Bayesian model to predict RA at set prevalence or incidence values, allowing comparison across studies and viruses. These predictions were then used to estimate the sequencing depth and concomitant cost required for pathogen detection using W-MGS with or without use of a hybridization-capture enrichment panel. 
### Findings
After controlling for variation in local infection rates, relative abundance varied by orders of magnitude across studies for a given virus. For instance, a local SARS-CoV-2 weekly incidence of 1% corresponds to predicted SARS-CoV-2 relative abundance ranging from 3·8 × 10-10 to 2·4 × 10-7 across studies, translating to orders-of-magnitude variation in the cost of operating a system able to detect a SARS-CoV-2-like pathogen at a given sensitivity. Use of a respiratory virus enrichment panel in two studies dramatically increased predicted relative abundance of SARS-CoV-2, lowering yearly costs by 24- to 29-fold for a system able to detect a SARS-CoV-2-like pathogen before reaching 0.01% cumulative incidence.
### Interpretation
The large variation in viral relative abundance after controlling for epidemiological factors indicates that other sources of inter-study variation, such as differences in sewershed hydrology and lab protocols, have a substantial impact on the sensitivity and cost of W-MGS. Well-chosen hybridization capture panels can dramatically increase sensitivity and reduce cost for viruses in the panel, but may reduce sensitivity to unknown or unexpected pathogens.
### Funding 
Wellcome Trust; Open Philanthropy; Musk Foundation



### Repository structure
This work is split across two repositories: This repo contains scripts to collect incidence and prevalence estimates and to build the model. The computational pipeline to determine relative abundance from existing data can be found under https://github.com/naobservatory/mgs-workflow/tree/2.1.0.

Scripts to create the manuscript figures are found in the `figures/` directory. Scripts to generate epidemiological estimates are in the `pathogens/` directory.

### Working with prevalence data

In python, run `import pathogens` and then iterate over `pathogens.pathogens`. Each pathogen implements an `estimate_prevalences` method which gives one or more estimates.

Run `./summarize.py` to get an overview of the data.

### Downloading the data

To download the data used in the tables and figures presented in this repo, visit https://doi.org/10.6084/m9.figshare.28395104.v1

### Statistical model

For an overview of the statistical model see [model.md](model.md).

To fit the model, run `./fit.py`. This will create, under `model_output/`:

* `input.tsv`, a table of the input data to the model
* `fits.tsv`, a table of samples from the posterior distribution of model parameters
* `fits_summary.tsv`, a table listing summary statistics of the posterior distributions of model parameters
* `fig/`, a directory containing a large number plots of posterior distributions and samples from posterior predictive distributions
  (see [model.md](model.md) for details)

Once the model has been fit, run `./plot_summaries.py` to create plots of the posterior distribution of $RA(1%)$.

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
