# andikoni_wtp_public

Replication code and data for

> TODO author list (2026). *TODO paper title.* TODO journal.
> TODO DOI

The study asks whether tree planting on private and on community land are
substitutes drawing on one underlying prosocial orientation, or two distinct
behaviours. Willingness to give up trees was elicited from 71 respondents with a
bid ladder for each type of land, and the two outcomes are modelled jointly:
once with an unrestricted bivariate normal residual covariance, and once with a
one-factor latent model whose loadings on the two outcomes are equal and
opposite.

The code is released under the MIT License, so fork and adapt it freely. The
data and figures carry the paper's license, CC BY-NC-ND 4.0: reuse with
attribution, non-commercial, no redistribution of modified versions. Both are in
[LICENSE](LICENSE). Please cite the paper; `CITATION.cff` carries the reference
in machine-readable form.

## Quick start

Install the packages listed under *Software*, then from the repository root:

```
# a fast end-to-end check (minutes)
WTP_TEST=1 Rscript --vanilla 0_RUN_ALL.R

# the paper's settings
Rscript --vanilla 0_RUN_ALL.R
```

That is the whole pipeline: it fits the models, writes the figures and writes
the comparison tables. Each step can also be run on its own, by name, in the
order below.

## What is here

| Path | Contents |
|---|---|
| `0_RUN_ALL.R` | Runs the whole pipeline. Start here. |
| `1_MAKE_*.R` … `4_MAKE_*.R` | The four pipeline steps, numbered in run order (table below). |
| `code/stan_models/` | The three Stan programs. `base_model.stan` is the unrestricted MVN benchmark, `latent.stan` the one-factor model, `latent_scale_model.stan` puts the predictors on the latent scale. |
| `code/prepare_data.R` | Builds the Stan input list from the released CSV. Sourced by steps 1–3 so they cannot disagree about the model inputs. |
| `code/functions/utility.R` | Shared helpers, and the one place the sampler settings and test-mode switch are defined. |
| `code/si_validation_sweep.R` | The supplementary simulation study, driven by step 4. |
| `data/wtp_data.csv` | The analysis dataset: 71 respondents, 7 columns (see *Data*). |
| `data/CODEBOOK.md` | Variable definitions, the elicitation procedure, and what was removed for release. |
| `data/data_dictionary.csv` | One row per column: type, unit, meaning, observed range. |
| `figures/` | The paper's figures as produced by the pipeline on our machine, for comparison with your own run. |
| `LICENSE`, `CITATION.cff` | Terms of reuse and the reference to cite. |

## Pipeline

Run from the repository root. Each step loads its own packages, sources
`code/functions/utility.R`, and reads only files an earlier step wrote.

| Step | Script | Reads | Writes | Time at full settings |
|---|---|---|---|---|
| 1 | `1_MAKE_fit_models.R` | `data/wtp_data.csv`, 3 Stan programs | `data/posterior_{base,latent,latent_pred}.RDS` (not in git) | minutes |
| 2 | `2_MAKE_figures.R` | step 1 output | `figures/{base,predictors,latent}.pdf` | seconds |
| 3 | `3_MAKE_model_comparison.R` | step 1 output | `data/model_comparison.csv`, `data/rank1_summary.csv`, `figures/{rank1_share,rq3}.pdf` | seconds |
| 4 | `4_MAKE_si_validation_sweep.R` | 2 Stan programs (all data simulated) | `si_output/sweep_results.{csv,rds}`, `si_output/fig{1,2,3}_*.png` | days (30 cells x 2 models; see below) |

Steps 1–3 are the main text and finish in about a minute together. Step 4 is
the supplement and is independent of the data — it simulates everything it
needs.

Step 4 is by far the expensive part and runs the cells sequentially. One cell
at lambda = 0.5, N = 250 takes about 45 minutes (roughly 10 minutes for the
latent model and 35 for the MVN benchmark, which is the slower of the two), and
cost grows with N. The full 10 x 3 grid is therefore a one- to two-day run on a
single machine. Do a `WTP_TEST=1` pass first, and consider trimming `N_GRID` in
`code/si_validation_sweep.R` if you only need the shape of the crossover.

Full settings for every fit of the three main models: 4 chains, 250 warmup and
250 sampling iterations (`code/functions/utility.R`). The sweep uses 500/500.

### Test mode

`WTP_TEST=1` runs everything in minutes: every model uses 20 warmup and 20
sampling iterations and the sweep keeps only its corner cells. It proves every
path, file and figure on your machine; the posteriors it produces are not
interpretable. Unset it for the paper's settings.

```
# bash / zsh
WTP_TEST=1 Rscript --vanilla 0_RUN_ALL.R

# PowerShell
$env:WTP_TEST = "1"; Rscript --vanilla 0_RUN_ALL.R
Remove-Item Env:WTP_TEST   # back to full settings
```

### Environment variables

| Variable | Effect |
|---|---|
| `WTP_TEST=1` | Test mode, described above. |
| `WTP_SKIP_SI=1` | Skip step 4. Steps 1–3 reproduce every main-text figure on their own. |

### What to expect while a fit runs

CmdStan prints occasional `Informational Message: The current Metropolis
proposal is about to be rejected` lines during warmup, typically from
`lkj_corr_cholesky_lpdf`. These are expected while the sampler adapts and do not
indicate a problem unless they persist through the sampling phase.

## Software

The pipeline was last run end to end on R 4.6.1 (Windows) with:

| Package | Version |
|---|---|
| CmdStan | 2.37.0 |
| `cmdstanr` | 0.9.0 |
| `posterior` | 1.7.1 |
| `loo` | 2.10.1 |
| `rethinking` | 2.42 |
| `MASS` | 7.3-65 |
| `dplyr` | 1.2.1 |
| `ggplot2` | 4.0.3 |

`rethinking` is not on CRAN; install it from GitHub with
`devtools::install_github("rmcelreath/rethinking")`. CmdStan is installed with
`cmdstanr::install_cmdstan()`. `ggplot2` is optional: step 4 falls back to base
graphics without it.

Step 3 reads a `model` column from `loo::loo_compare()`, which `loo` 2.10
introduced; it falls back to row names on older versions.

## Data

`data/wtp_data.csv` holds 71 respondents and 7 columns.
[`data/CODEBOOK.md`](data/CODEBOOK.md) documents each column, the elicitation
procedure, and how the file was de-identified. `data/data_dictionary.csv` gives
the same column list in machine-readable form.

The release is pseudonymous rather than anonymous: names and the panel and
household identifiers are absent, but the respondents are not k-anonymous on
`{sex, marital status, age, education}`, because coarsening those variables
would prevent the published results from reproducing. The codebook states the
limit in full; respect it if you redistribute or link the file.

## Note on the model comparison

The latent model's `log_lik` is conditional on the per-person latent `T_i`,
which biases a person-level LOO comparison toward the latent model. Step 4 sets
`RECOMPUTE_MARGINAL_LL = TRUE` and rebuilds the marginal per-person log
likelihood in R from the posterior draws. The same caveat applies to the LOO
table written by step 3, which is a descriptive summary rather than a
model-selection criterion.
