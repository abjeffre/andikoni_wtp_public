###################################################################
############ 4. SI VALIDATION SWEEP ###############################
# Reads  code/stan_models/{latent,base_model}.stan  (all data is simulated)
# Writes si_output/sweep_results.{csv,rds},
#        si_output/fig{1,2,3}_*.png
#
# Sweeps the dependence parameter lambda from 0 to 0.9 crossed with sample
# size, feeds the IDENTICAL simulated dataset to the latent and MVN models,
# and compares them by paired PSIS-LOO. At full settings this is 33 cells x 2
# models and takes hours; WTP_TEST=1 reduces it to the corner cells.
#
# This is the CORRECTED sweep: it recomputes the marginal per-person log_lik
# in R (RECOMPUTE_MARGINAL_LL = TRUE) rather than trusting the conditional
# log_lik that latent.stan emits. See README "Known issues".

library(cmdstanr)
library(rethinking)
library(posterior)
library(loo)
library(MASS)
library(dplyr)
source("code/functions/utility.R")

dir.create("si_output", showWarnings = FALSE)
source("code/si_validation_sweep.R")

cat("step 4 complete: si_output/\n")
