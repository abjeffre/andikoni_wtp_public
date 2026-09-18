###################################################################
############ 0. RUN THE WHOLE PIPELINE ############################
# The single entry point. Runs steps 1-4 in order from the repository root:
#
#   Rscript --vanilla 0_RUN_ALL.R
#
# Environment variables:
#   WTP_TEST=1     every model uses 20 warmup / 20 sampling iterations and the
#                  sweep keeps only its corner cells, so the whole pipeline
#                  finishes in minutes. The posteriors it produces are
#                  meaningless; it exists to prove every path and figure.
#   WTP_SKIP_SI=1  skip step 4, the supplementary sweep, which takes hours at
#                  full settings. Steps 1-3 reproduce every figure in the main
#                  text on their own.

steps <- c(
  "1_MAKE_fit_models.R",
  "2_MAKE_figures.R",
  "3_MAKE_model_comparison.R",
  "4_MAKE_si_validation_sweep.R"
)

if (nzchar(Sys.getenv("WTP_SKIP_SI"))) {
  steps <- setdiff(steps, "4_MAKE_si_validation_sweep.R")
}

if (nzchar(Sys.getenv("WTP_TEST"))) {
  message("WTP_TEST is set: running at test settings. Results are meaningless.")
} else {
  message("Running at full settings. Step 4 takes hours; WTP_SKIP_SI=1 skips it.")
}

started <- Sys.time()
for (step in steps) {
  message("\n=== ", step, " ===")
  t0 <- Sys.time()
  # Each step is run in its own environment so that objects cannot leak
  # between them; every step reads what it needs from disk.
  source(step, local = new.env(), echo = FALSE)
  message(sprintf(
    "--- %s finished in %.1f min", step,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}

message(sprintf(
  "\nPipeline complete in %.1f min.",
  as.numeric(difftime(Sys.time(), started, units = "mins"))
))
