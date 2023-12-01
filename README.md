## Predicting Virus Relative Abundance in Wastewater
Simon L. Grimm, Jeff T. Kaufman, Daniel Rice, Michael M. McLaren, Charlie Whittaker, William J. Bradshaw & Kevin M. Esvelt1,2,3 (insert Charlie)

### Abstract
Detecting novel pathogens at an early stage requires robust early warning that is both sensitive and pathogen-agnostic. Wastewater metagenomic sequencing (W-MGS) could enable highly pathogen-agnostic disease monitoring, but its sensitivity and financial feasibility are dependent on the relative abundance of novel pathogen sequences in W-MGS data. Here we collate W-MGS data from a diverse range of studies to characterize the relative abundance of known viruses in wastewater samples. We develop a Bayesian statistical model to integrate a subset of these data with epidemiological estimates of incidence and prevalence, and use it to estimate the expected relative abundance of different viral pathogens for a given prevalence or incidence in the community. Our results reveal pronounced variation between sites and studies, with estimates differing by one to three orders of magnitude for the same pathogen. For example, the expected relative abundance of SARS-CoV-2 at weekly incidence of 1% of the population varied between 10-7 and 10-10. Integrating these estimates with a simple cost model highlights substantial variation in the volume of W-MGS required to detect these pathogens. The mean sequencing cost of identifying 100 reads from a new SARS-CoV-2-like pathogen at 1% cumulative incidence was $2,579,000. A Norovirus-like pathogen which sheds more would require sequencing costs of $18,000. This model, and its parameter estimates, represent an important resource for future investigation into the performance of wastewater MGS, and can be extended to incorporate new wastewater datasets as they become available.

In this project we're attempting to understand how the prevalence of
human pathogens translates into the relative abundance we see in wastewater
metagenomics.  This work is split across two repositories: in this repo we're
collecting prevalence estimates and building the model, while determining
relative abundances from existing data is in the
[mgs-pipeline](https://github.com/naobservatory/mgs-pipeline) repo.

### Working with prevalence data

In python, run `import pathogens` and then iterate over `pathogens.pathogens`.
Each pathogen implements an `estimate_prevalences` method which gives one or
more estimates.

Run `./summarize.py` to get an overview of the data.

### Statistical model

For an overview of the statistical model see [model.md](model.md).

To fit the model, run `./fit.py`. This will create:

* `input.tsv`, a table of the input data to the model
* `fits.tsv`, a table of samples from the posterior distribution of model parameters
* `fits_summary.tsv`, a table listing summary statistics of the posterior distributions of model parameters
* `fig/`, a directory containing a large number plots of posterior distributions and samples from posterior predictive distributions
  (see [model.md](model.md) for details)

Once the model has been fit, run `./plot_summaries.py` to create plots of the posterior distribution of $RA(1\perthousand)$ for the write-up:

* `fig/incidence-violin.{pdf,png}`, posteriors for all incidence viruses
* `fig/prevalence-violin.{pdf,png}`, posteriors for all prevalence viruses
* `fig/by_location_incidence-violin.{pdf,png}`, posteriors separated by location for the most common incidence viruses

### Development

#### Making changes

Create a branch named yourname-purpose and push your changes to it.  Then
create a pull request.  Use the "request review" feature to ask for a review:
all PRs need to be reviewed by someone else, and for now include Jeff and Simon
on all PRs unless they're OOO.

Once your change has been approved by your reviewers and passes
[presubmit](#presubmit), you can merge it.  Don't merge someone else's PR
without confirming with them: they may have other changes they've realized they
needed to make, or a tricky branch structure that needs to be resolved in a
particular order.

Handle incoming reviews at least twice a day (morning and afternoon) -- slow
reviews add a lot of friction.  As a PR author you can avoid this friction by
creating another branch that diverges from the code you have under review; ask
Jeff to show you how if you're interested.  Configure [notification
routing](https://github.com/settings/notifications/custom_routing) on github so
that work-related notifications go to your work account.

#### Testing

Run `./test.py`

#### Presubmit

Before creating a PR or submitting code, run `./check.sh`.  It will run tests
than check your types and formatting.  This also runs automatically on GitHub
when you make PRs, but it's much faster to catch problems locally first.

If `./check.sh` complains about formatting or import sorting, you can fix this
automatically with `./check.sh --fix`.

#### Installing pystan

Pystan should be installed along with the other requirements when you run:
```
python -m pip install -r requirements-dev.txt"
```
However, on some non-Linux systems (including M2 Macbooks), one of `pystan`'s dependencies,`httpstan`, may fail to install.
To get around this problem, you can [install httpstan from source](https://httpstan.readthedocs.io/en/latest/installation.html#installation-from-source).
Once it is built and installed, you can then install the requirements file as above.
(Note that you can clone the `httpstan` repo anywhere on your computer.
I recommend doing it outside of the `p2ra` repo directory to that git doesn't try to track it.)
