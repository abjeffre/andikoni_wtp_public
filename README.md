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
| `data/CODEBOOK.md` | Variable definitions, elicitation, corrections, and what was removed for release. |
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
| 4 | `4_MAKE_si_validation_sweep.R` | 2 Stan programs (all data simulated) | `si_output/sweep_results.{csv,rds}`, `si_output/fig{1,2,3}_*.png` | hours (33 cells x 2 models) |

Steps 1–3 are the main text. Step 4 is the supplement and is independent of the
data — it simulates everything it needs.

Full settings for every fit of the three main models: 4 chains, 250 warmup and
250 sampling iterations (`code/functions/utility.R`). The sweep uses 500/500.

### Running

```
Rscript --vanilla 0_RUN_ALL.R
```

That is the whole pipeline. To run a single step instead, run it by name from
the repository root — each is self-contained.

### Test run first

`WTP_TEST=1` runs everything in minutes: every model uses 20 warmup and 20
sampling iterations and the sweep keeps only its corner cells. This proves
every path, file and figure on your machine; the posteriors it produces are
meaningless. Unset it for the paper's settings.

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

`data/wtp_data.csv` holds 71 respondents and 7 columns. `data/CODEBOOK.md`
documents each one, together with the elicitation procedure and the two rounds
of corrections applied to the wealth variable.

### De-identification and residual risk

The source workbook contains respondent **names**, an individual panel
identifier (`PE`) whose prefix also encodes the household, and a household
identifier (`PESU2018`). None of these are in this repository, in any form, at
any point in its history. Row order was randomised under a fixed seed, because
in the source workbook rows are sorted by `PE` and `PE` order is alphabetical by
name, so row position alone was identifying.

**This release is pseudonymous, not anonymous, and we state the limit plainly.**
Across the quasi-identifiers `{sex, marital status, age, education}` the 71
respondents fall into 61 distinct combinations: 53 respondents are unique on
those four variables, and every respondent is in a class of fewer than five.
Anyone who already knows a participant's age, sex, marital status and schooling
could likely locate their row and read off their wealth and their two bids. We
judged this acceptable because coarsening age and wealth — the two continuous
model predictors — would prevent the published results from reproducing, and
because the study population is not identified at village level in the paper.
If you intend to redistribute or link this file, that is the constraint to
respect.

## Known issues

Three points where the code as published departs from what was probably
intended. All three are carried over **verbatim** from the original analysis
script so that this repository reproduces the published numbers; each is marked
`# FLAGGED` at its location.

1. **Third panel of `figures/predictors.pdf`** (`2_MAKE_figures.R`). The
   education and wealth rows of the "Contrasts (Comm − Priv)" panel index
   `post$bE[2]` and `post$bW[2]` without a comma, so a single posterior draw is
   recycled instead of the second outcome's column (`post$bE[, 2]`). The three
   rows above them use `[, 2]` as intended.
2. **MVN residual reconstruction** (`extract_mu_residuals` in
   `code/functions/utility.R`). The linear predictor is rebuilt from `a`, `bW`,
   `bA` and `bA2` only; `base_model.stan` also contains a sex effect and an
   education effect, which therefore remain in the residuals.
3. **Rank-1 share** (`3_MAKE_model_comparison.R`). The originally reported
   statistic pairs the MVN private residual with the *latent* community
   residual. The step now also prints both within-model pairs alongside it.

A fourth, more consequential point concerns the supplement. The latent model's
`log_lik` is conditional on the per-person latent `T_i`, which biases a
person-level LOO comparison toward the latent model. Step 4 therefore sets
`RECOMPUTE_MARGINAL_LL = TRUE` and rebuilds the marginal per-person log
likelihood in R from the posterior draws. The same caveat applies to the LOO
table written by step 3, which is reported as a descriptive summary rather than
a model-selection criterion.
