############################################################
# SI VALIDATION SWEEP: Latent factor model vs MVN benchmark
# ----------------------------------------------------------
# Design (locked):
#   * 2 continuous outcomes (priv, pub).
#   * ONE matched latent DGP per cell; lambda swept 0 -> 0.9.
#   * The IDENTICAL dataset is fed to BOTH models (paired comparison).
#   * Crossed with sample size N. 1 replicate/cell (loop carries a
#     replicate index so you can raise N_REP later for bands).
#
# Equivalence used throughout (conditional on covariates, Var(T)=1):
#     Sigma_latent = [[ lambda^2 + sigma^2 ,   -lambda^2        ],
#                     [   -lambda^2        ,  lambda^2 + sigma^2 ]]
#   => the MVN-equivalent off-diagonal is  s12_true = -lambda^2,
#      and the residual correlation is rho_true = -lambda^2/(lambda^2+sigma^2).
#
# Metrics (figure priority):
#   (1) MAIN : paired PSIS-LOO elpd_diff (latent - mvn) +/- SE, vs lambda, by N
#   (2)        rank1_share of residual covariance (your diagnostic)
#   (3)        recovery of lambda (latent) and s12 (MVN) vs matched truth
#   (4)        posterior residual correlation vs truth
#
# !! VALIDITY NOTE (read this) ----------------------------------------------
# With a PER-PERSON latent T_i, person-level LOO is contaminated unless the
# latent model's log_lik[i] integrates T_i OUT (marginal bivariate normal).
# If latent.stan emits the CONDITIONAL (given-T) log_lik, the elpd comparison
# is biased toward the latent model. Two options below:
#   A) Fix it in latent.stan: in generated quantities emit the marginal
#      log_lik[i] = multi_normal_lpdf({priv_i,pub_i} | mu_i, Sigma_latent).
#   B) Set RECOMPUTE_MARGINAL_LL <- TRUE to recompute both models' per-person
#      log_lik in R from posterior draws (airtight, ignores Stan's log_lik).
#      Requires the parameter names in the EXTRACTORS section to match.
# --------------------------------------------------------------------------
############################################################

suppressPackageStartupMessages({
  library(cmdstanr)
  library(rethinking)
  library(posterior)
  library(loo)
  library(MASS)
  library(dplyr)
})
have_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

############################################################
# 0. CONFIG  --  edit these to match your repo
############################################################

STAN_DIR <- "code/stan_models"
LATENT_STAN <- file.path(STAN_DIR, "latent.stan")
MVN_STAN    <- file.path(STAN_DIR, "base_model.stan")

OUT_DIR <- "si_output"

# --- sweep grid (WTP_TEST=1 keeps only the corner cells) ---
LAMBDA_GRID <- if (WTP_TEST) SWEEP_LAMBDA_TEST else SWEEP_LAMBDA_FULL
N_GRID      <- if (WTP_TEST) SWEEP_N_TEST else SWEEP_N_FULL
N_REP       <- 1                       # raise later for error bands
SIGMA_TRUE  <- 0.5                     # idiosyncratic SD in the DGP

# --- sampler ---
CHAINS <- SAMPLER_CHAINS
WARMUP <- if (WTP_TEST) 20 else 500
SAMPLING <- if (WTP_TEST) 20 else 500
ADAPT_DELTA <- 0.95; MAX_TREEDEPTH <- if (WTP_TEST) 5 else 12
BASE_SEED <- 20260601

# --- log-lik handling ---
LOGLIK_VAR <- "log_lik"        # per-person log_lik variable name in GQ
RECOMPUTE_MARGINAL_LL <- TRUE  # set TRUE to recompute marginal LL in R (option B)

# --- parameter / fitted-value names emitted by your Stan files ---
# (used for residuals, recovery, and optional marginal-LL recompute)
NAMES <- list(
  latent_mu_priv = "mu_priv_hat",   # fitted mean for priv, latent model
  latent_mu_pub  = "mu_pub_hat",    # fitted mean for pub,  latent model
  latent_lambda  = "lambda",        # scalar loading
  latent_sigma   = "sigma",         # length-2 idiosyncratic SDs (or scalar)
  latent_T       = "T",             # per-person latent (for recompute option)
  mvn_a   = "a", mvn_bW = "bW", mvn_bA = "bA", mvn_bA2 = "bA2", # MVN mean coefs
  mvn_Sigma = "L_Sigma"             # 2x2 residual covariance (or set Cholesky below)
)
MVN_SIGMA_IS_CHOLESKY <- TRUE       # TRUE if your MVN model stores L (lower-tri)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

############################################################
# 1. DATA-GENERATING PROCESS (single, matched, lambda-controlled)
############################################################
# One shared mean model; dependence enters ONLY through the latent factor.
# Both Stan models receive this same data, so any mean misspecification is
# common to both and cancels in the PAIRED elpd_diff.

simulate_dgp_latent <- function(N, lambda, sigma = SIGMA_TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  wealth <- rnorm(N)
  age    <- runif(N, 20, 60)
  sex    <- sample(1:2, N, replace = TRUE)
  edu    <- sample(1:5, N, replace = TRUE)
  Tlat   <- rnorm(N)                       # Var(T) = 1

  a   <- c(0.5, 0.5)
  bW  <- c(0.4, 0.4)
  bA  <- c(0.1, 0.1)
  bA2 <- c(-0.001, -0.001)
  bS  <- matrix(c(0.2, -0.2,
                  0.1,  0.3), nrow = 2, byrow = TRUE)

  base_priv <- a[1] + bS[1, sex] + bW[1]*wealth + bA[1]*age + bA2[1]*age^2
  base_pub  <- a[2] + bS[2, sex] + bW[2]*wealth + bA[2]*age + bA2[2]*age^2

  priv <- rnorm(N, base_priv - lambda*Tlat, sigma)
  pub  <- rnorm(N, base_pub  + lambda*Tlat, sigma)

  df <- data.frame(priv, pub, wealth, age, sex, edu, T = Tlat)

  # matched MVN truth (what the MVN model should recover)
  attr(df, "truth") <- list(
    lambda  = lambda,
    sigma   = sigma,
    s12     = -lambda^2,                                   # off-diagonal cov
    s_diag  = lambda^2 + sigma^2,                          # marginal var
    rho     = -lambda^2 / (lambda^2 + sigma^2)             # residual corr
  )
  df
}

make_standata <- function(df, N) {
  list(
    N = N, priv = df$priv, pub = df$pub,
    wealth = df$wealth, age = df$age, sex = df$sex,
    K_edu = length(unique(df$edu)), edu = df$edu,
    alpha = rep(1, length(unique(df$edu)) - 1)
  )
}

############################################################
# 2. COMPILE ONCE
############################################################
message("Compiling models ...")
latent_model <- cmdstan_model(LATENT_STAN)
mvn_model    <- cmdstan_model(MVN_STAN)

fit_one <- function(model, standata, seed) {
  t0 <- Sys.time()
  fit <- model$sample(
    data = standata, seed = seed,
    chains = CHAINS, parallel_chains = CHAINS,
    iter_warmup = WARMUP, iter_sampling = SAMPLING,
    adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH,
    refresh = 0, show_messages = FALSE
  )
  diag <- fit$diagnostic_summary()
  summ <- suppressWarnings(fit$summary())
  list(
    fit = fit,
    runtime_s = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    n_div     = sum(diag$num_divergent),
    max_rhat  = max(summ$rhat, na.rm = TRUE),
    min_ess   = min(summ$ess_bulk, na.rm = TRUE)
  )
}

############################################################
# 3. METRIC HELPERS
############################################################

# --- per-person log_lik -> loo object ---
get_loglik_array <- function(fit) fit$draws(LOGLIK_VAR)  # iter x chain x N

loo_from_fit <- function(fit) {
  ll <- get_loglik_array(fit)
  ll_mat <- posterior::as_draws_matrix(ll)               # (iter*chain) x N
  cid <- rep(seq_len(dim(ll)[2]), each = dim(ll)[1])
  r_eff <- loo::relative_eff(exp(ll_mat), chain_id = cid)
  list(loo = suppressWarnings(loo::loo(ll_mat, r_eff = r_eff)),
       pw  = NULL)
}

# Optional: recompute MARGINAL per-person log_lik in R (validity option B).
# Returns (iter*chain) x N matrix for each model.
recompute_marginal_ll <- function(fit, df, type = c("latent","mvn")) {
  type <- match.arg(type)
  post <- extract.samples(fit)
  Y <- cbind(df$priv, df$pub); N <- nrow(Y)
  if (type == "latent") {
    lam   <- post[[NAMES$latent_lambda]]      # S
    T_lat <- post[[NAMES$latent_T]]           # S x N
    sig   <- post[[NAMES$latent_sigma]]       # S (scalar) or S x 2
    # mu_priv_hat/mu_pub_hat include lambda*T[i]; recover T-excluded marginal means
    mu_p_raw <- post[[NAMES$latent_mu_priv]]  # S x N, T-included
    mu_c_raw <- post[[NAMES$latent_mu_pub]]   # S x N, T-included
    S <- nrow(mu_p_raw)
    lam_mat <- matrix(rep(lam, N), nrow = S, ncol = N)
    mu_p <- mu_p_raw + lam_mat * T_lat        # base_priv = mu_priv_hat + lambda*T
    mu_c <- mu_c_raw - lam_mat * T_lat        # base_pub  = mu_pub_hat  - lambda*T
    ll <- matrix(NA_real_, S, N)
    for (s in seq_len(S)) {
      sg <- if (is.null(dim(sig))) c(sig[s], sig[s]) else sig[s, ]
      Sigma <- matrix(c(lam[s]^2 + sg[1]^2, -lam[s]^2,
                        -lam[s]^2,           lam[s]^2 + sg[2]^2), 2, 2)
      ll[s, ] <- mvtnorm_logpdf(Y, cbind(mu_p[s, ], mu_c[s, ]), Sigma)
    }
    ll
  } else {
    mu_p <- reconstruct_mvn_mu(post, df, 1)
    mu_c <- reconstruct_mvn_mu(post, df, 2)
    Sig  <- get_mvn_sigma_draws(post)         # list of 2x2 per draw
    S <- nrow(mu_p); ll <- matrix(NA_real_, S, N)
    for (s in seq_len(S)) ll[s, ] <- mvtnorm_logpdf(Y, cbind(mu_p[s,], mu_c[s,]), Sig[[s]])
    ll
  }
}

mvtnorm_logpdf <- function(Y, MU, Sigma) {
  # vectorised bivariate-normal log density over rows of Y
  d  <- Y - MU
  iS <- solve(Sigma); ld <- determinant(Sigma, logarithm = TRUE)$modulus
  q  <- rowSums((d %*% iS) * d)
  as.numeric(-0.5 * (2*log(2*pi) + ld + q))
}

loo_from_matrix <- function(ll_mat, chains = CHAINS) {
  cid <- rep(seq_len(chains), each = nrow(ll_mat)/chains)
  r_eff <- loo::relative_eff(exp(ll_mat), chain_id = cid)
  suppressWarnings(loo::loo(ll_mat, r_eff = r_eff))
}

# --- signed paired elpd difference (latent - mvn) with paired SE ---
paired_elpd_diff <- function(loo_lat, loo_mvn) {
  d <- loo_lat$pointwise[, "elpd_loo"] - loo_mvn$pointwise[, "elpd_loo"]
  n <- length(d)
  list(elpd_diff = sum(d), se_diff = sqrt(n) * sd(d),
       n_pareto_bad = sum(pareto_k_values(loo_lat) > 0.7) +
                      sum(pareto_k_values(loo_mvn) > 0.7))
}

# --- residual extraction (fitted means) ---
reconstruct_mvn_mu <- function(post, df, which_outcome) {
  j   <- which_outcome
  N   <- nrow(df)
  a   <- post[[NAMES$mvn_a]];   bW  <- post[[NAMES$mvn_bW]]
  bA  <- post[[NAMES$mvn_bA]];  bA2 <- post[[NAMES$mvn_bA2]]
  bS_d  <- post[["bS"]]         # S x 2 x 2  (outcome x sex level)
  edu_d <- post[["edu_effect"]] # S x 2 x K_edu (outcome x edu level)
  # draws x N: intercept + covariate slopes + sex effect + education effect
  outer(a[, j], rep(1, N)) +
    outer(bW[, j], df$wealth) + outer(bA[, j], df$age) + outer(bA2[, j], df$age^2) +
    bS_d[, j, ][, df$sex] +     # S x N sex contribution
    edu_d[, j, ][, df$edu]      # S x N education contribution
}

extract_residuals <- function(fit, df, type = c("latent","mvn")) {
  type <- match.arg(type); post <- extract.samples(fit)
  if (type == "latent") {
    mu_p <- apply(post[[NAMES$latent_mu_priv]], 2, mean)
    mu_c <- apply(post[[NAMES$latent_mu_pub]],  2, mean)
  } else {
    mu_p <- apply(reconstruct_mvn_mu(post, df, 1), 2, mean)
    mu_c <- apply(reconstruct_mvn_mu(post, df, 2), 2, mean)
  }
  list(res_p = df$priv - mu_p, res_c = df$pub - mu_c)
}

rank1_share <- function(r1, r2) {
  sv <- svd(cov(cbind(r1, r2)))$d
  sv[1] / sum(sv)
}

# --- covariance / recovery extraction ---
get_mvn_sigma_draws <- function(post) {
  S <- post[[NAMES$mvn_Sigma]]
  # accept draws x 2 x 2 array, or Cholesky
  out <- vector("list", dim(S)[1])
  for (s in seq_len(dim(S)[1])) {
    M <- S[s, , ]
    if (MVN_SIGMA_IS_CHOLESKY) M <- M %*% t(M)
    out[[s]] <- M
  }
  out
}

recovery_row <- function(fit_lat, fit_mvn, truth) {
  post_l <- extract.samples(fit_lat)
  out <- list()

  # lambda (latent)
  lam <- post_l[[NAMES$latent_lambda]]
  qi  <- quantile(lam, c(.05,.95))
  out$lambda_mean <- mean(lam)
  out$lambda_lo   <- qi[[1]]; out$lambda_hi <- qi[[2]]
  out$lambda_cov  <- truth$lambda >= qi[[1]] && truth$lambda <= qi[[2]]

  # s12 (mvn) + residual corr from MVN Sigma draws
  Sig <- tryCatch(get_mvn_sigma_draws(extract.samples(fit_mvn)), error = function(e) NULL)
  if (!is.null(Sig)) {
    s12 <- sapply(Sig, function(M) M[1,2])
    rho <- sapply(Sig, function(M) M[1,2]/sqrt(M[1,1]*M[2,2]))
    qs <- quantile(s12, c(.05,.95)); qr <- quantile(rho, c(.05,.95))
    out$s12_mean <- mean(s12); out$s12_lo <- qs[[1]]; out$s12_hi <- qs[[2]]
    out$s12_cov  <- truth$s12 >= qs[[1]] && truth$s12 <= qs[[2]]
    out$rho_mean <- mean(rho); out$rho_lo <- qr[[1]]; out$rho_hi <- qr[[2]]
    out$rho_cov  <- truth$rho >= qr[[1]] && truth$rho <= qr[[2]]
  } else {
    warning("Could not read MVN Sigma draws (set NAMES$mvn_Sigma / MVN_SIGMA_IS_CHOLESKY).")
  }
  out
}

############################################################
# 4. SWEEP
############################################################
results <- list(); ri <- 0
total <- length(LAMBDA_GRID) * length(N_GRID) * N_REP
message(sprintf("Sweep: %d cells x 2 models = %d fits", total, 2*total))

for (lam in LAMBDA_GRID) for (N in N_GRID) for (rep in seq_len(N_REP)) {
  cell_seed <- BASE_SEED + ri
  df <- simulate_dgp_latent(N, lambda = lam, seed = cell_seed)
  truth <- attr(df, "truth")
  sdat  <- make_standata(df, N)

  fl <- fit_one(latent_model, sdat, seed = cell_seed + 1)
  fm <- fit_one(mvn_model,    sdat, seed = cell_seed + 2)

  # (1) ELPD via PSIS-LOO -------------------------------------------------
  if (RECOMPUTE_MARGINAL_LL) {
    loo_l <- loo_from_matrix(recompute_marginal_ll(fl$fit, df, "latent"))
    loo_m <- loo_from_matrix(recompute_marginal_ll(fm$fit, df, "mvn"))
  } else {
    loo_l <- loo_from_fit(fl$fit)$loo
    loo_m <- loo_from_fit(fm$fit)$loo
  }
  cmp <- paired_elpd_diff(loo_l, loo_m)

  # (2) rank1_share -------------------------------------------------------
  rl <- extract_residuals(fl$fit, df, "latent")
  rm <- extract_residuals(fm$fit, df, "mvn")
  r1_lat <- rank1_share(rl$res_p, rl$res_c)
  r1_mvn <- rank1_share(rm$res_p, rm$res_c)

  # (3)+(4) recovery ------------------------------------------------------
  rec <- recovery_row(fl$fit, fm$fit, truth)

  ri <- ri + 1
  results[[ri]] <- data.frame(
    lambda = lam, N = N, rep = rep, seed = cell_seed,
    s12_true = truth$s12, rho_true = truth$rho,
    elpd_diff = cmp$elpd_diff, se_diff = cmp$se_diff,
    pareto_bad = cmp$n_pareto_bad,
    elpd_latent = loo_l$estimates["elpd_loo","Estimate"],
    elpd_mvn    = loo_m$estimates["elpd_loo","Estimate"],
    rank1_latent = r1_lat, rank1_mvn = r1_mvn,
    lambda_mean = rec$lambda_mean %||% NA, lambda_cov = rec$lambda_cov %||% NA,
    s12_mean = rec$s12_mean %||% NA, s12_cov = rec$s12_cov %||% NA,
    rho_mean = rec$rho_mean %||% NA, rho_cov = rec$rho_cov %||% NA,
    div_lat = fl$n_div, div_mvn = fm$n_div,
    rhat_lat = fl$max_rhat, rhat_mvn = fm$max_rhat,
    rt_lat = fl$runtime_s, rt_mvn = fm$runtime_s
  )
  message(sprintf("[%d/%d] lambda=%.2f N=%d  elpd_diff(lat-mvn)=%.2f +/- %.2f  | div %d/%d",
                  ri, total, lam, N, cmp$elpd_diff, cmp$se_diff, fl$n_div, fm$n_div))
}

res <- dplyr::bind_rows(results)
res$winner <- with(res, ifelse(elpd_diff >  2*se_diff, "latent",
                        ifelse(elpd_diff < -2*se_diff, "mvn", "tie")))

saveRDS(res, file.path(OUT_DIR, "sweep_results.rds"))
write.csv(res, file.path(OUT_DIR, "sweep_results.csv"), row.names = FALSE)
message("Saved tables to ", OUT_DIR)

############################################################
# 5. FIGURES
############################################################
fig_path <- function(f) file.path(OUT_DIR, f)

if (have_ggplot) {
  library(ggplot2)
  # (1) MAIN: elpd_diff vs lambda, faceted by N
  p1 <- ggplot(res, aes(lambda, elpd_diff)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_ribbon(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff), alpha = .2) +
    geom_line() + geom_point(aes(color = winner), size = 2) +
    facet_wrap(~ N, labeller = label_both) +
    labs(x = expression(lambda~"(  "*s[12]==-lambda^2*"  )"),
         y = "elpd_diff  (latent - MVN)",
         title = "Held-out predictive advantage vs dependence strength") +
    theme_bw()
  ggsave(fig_path("fig1_elpd_diff.png"), p1, width = 9, height = 3.2, dpi = 200)

  # (2) rank1_share
  r1 <- res |>
    tidyr::pivot_longer(c(rank1_latent, rank1_mvn), names_to="model", values_to="rank1")
  p2 <- ggplot(r1, aes(lambda, rank1, color = model)) +
    geom_line() + geom_point() + facet_wrap(~ N, labeller = label_both) +
    labs(x = expression(lambda), y = "rank-1 share of residual cov",
         title = "Residual structure left by each model") + theme_bw()
  ggsave(fig_path("fig2_rank1_share.png"), p2, width = 9, height = 3.2, dpi = 200)

  # (3) recovery of lambda and s12
  p3 <- ggplot(res, aes(lambda)) +
    geom_abline(slope = 1, linetype = 2) +
    geom_point(aes(y = lambda_mean), color = "steelblue") +
    geom_point(aes(y = -s12_mean), color = "firebrick", shape = 17) +
    facet_wrap(~ N, labeller = label_both) +
    labs(y = "posterior mean", x = expression(lambda~"(true)"),
         title = "Recovery: lambda (blue) and -s12 (red) vs truth") + theme_bw()
  ggsave(fig_path("fig3_recovery.png"), p3, width = 9, height = 3.2, dpi = 200)
  message("Saved figures (ggplot) to ", OUT_DIR)
} else {
  png(fig_path("fig1_elpd_diff.png"), width = 1100, height = 380)
  par(mfrow = c(1, length(N_GRID)))
  for (nn in N_GRID) {
    s <- res[res$N == nn, ]
    plot(s$lambda, s$elpd_diff, type = "b", ylim = range(c(s$elpd_diff-s$se_diff, s$elpd_diff+s$se_diff)),
         xlab = expression(lambda), ylab = "elpd_diff (latent - MVN)", main = paste("N =", nn))
    arrows(s$lambda, s$elpd_diff-s$se_diff, s$lambda, s$elpd_diff+s$se_diff, code=3, angle=90, length=.03)
    abline(h = 0, lty = 2)
  }
  dev.off(); par(mfrow = c(1,1))
  message("Saved base-R figure to ", OUT_DIR, " (install ggplot2 + tidyr for full panels)")
}

############################################################
# 6. CALIBRATION (run separately, NOT in the per-cell loop)
############################################################
# SBC is expensive and ranked last. Run on 2-3 representative cells only,
# e.g. lambda in {0, 0.4, 0.8}, by repeatedly: draw params ~ prior, simulate,
# fit, store rank statistic of each truth within its posterior. See loo / SBC
# tooling. Stub left intentionally; wire in once the sweep above is trusted.

cat("\nDONE. Headline crossover (where winner flips):\n")
print(res[, c("lambda","N","elpd_diff","se_diff","winner")])
