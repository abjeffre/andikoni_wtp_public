###################################################################
############ 3. MODEL COMPARISON AND RESIDUAL STRUCTURE ###########
# Reads  data/wtp_data.csv, data/posterior_{base,latent,latent_pred}.RDS
# Writes data/model_comparison.csv, data/rank1_summary.csv,
#        figures/rank1_share.pdf, figures/rq3.pdf

library(rethinking)
library(posterior)
library(loo)
source("code/functions/utility.R")
source("code/prepare_data.R")

post_base <- readRDS("data/posterior_base.RDS")
post_latent <- readRDS("data/posterior_latent.RDS")
post_latent_pred <- readRDS("data/posterior_latent_pred.RDS")

###################################################################
############ PSIS-LOO AND WAIC ####################################
# NOTE ON VALIDITY: the latent model's log_lik is conditional on the per-person
# latent T_i, so this person-level comparison is not a clean out-of-sample
# criterion for it. The supplementary sweep (step 4) is the calibrated version
# of this contest; see README "Known issues".

ll <- list(
  base = post_base$log_lik,
  latent = post_latent$log_lik,
  latent_pred = post_latent_pred$log_lik
)

loo_fits <- lapply(ll, function(m) loo::loo(m))
waic_fits <- lapply(ll, function(m) loo::waic(m))

model_comparison <- as.data.frame(loo::loo_compare(loo_fits))
# loo >= 2.10 returns a data frame carrying the model names in a "model" column;
# older versions return a matrix that carries them in the row names.
if (!"model" %in% names(model_comparison)) {
  model_comparison <- cbind(model = rownames(model_comparison), model_comparison)
}
model_comparison$waic <- sapply(
  as.character(model_comparison$model),
  function(m) waic_fits[[m]]$estimates["waic", "Estimate"]
)
rownames(model_comparison) <- NULL

print(model_comparison[, c("model", "elpd_diff", "se_diff", "elpd_loo", "p_loo", "waic")])
write.csv(model_comparison, "data/model_comparison.csv", row.names = FALSE)

###################################################################
############ RESIDUAL COVARIANCE STRUCTURE ########################

cov_mvn <- get_cov_mvn(post_base)
cov_lat <- get_cov_latent(post_latent)

phi_mvn <- rank1_post(cov_mvn)
phi_lat <- rank1_post(cov_lat)

rank1_summary <- data.frame(
  model = c("MVN", "Latent"),
  mean = c(mean(phi_mvn), mean(phi_lat)),
  lower95 = c(quantile(phi_mvn, 0.025), quantile(phi_lat, 0.025)),
  upper95 = c(quantile(phi_mvn, 0.975), quantile(phi_lat, 0.975)),
  row.names = NULL
)

lat <- extract_mu_residuals(post_latent, stan_data, "latent")
mvn <- extract_mu_residuals(post_base, stan_data, "mvn")

# Rank-1 share of the FITTED residual pair, within each model. (The original
# script paired the MVN private residual with the latent community residual,
# mixing the two models; corrected here to compare like with like.)
rank1_summary$residual_rank1 <- c(
  rank1_share(mvn$res_p, mvn$res_c),
  rank1_share(lat$res_p, lat$res_c)
)
write.csv(rank1_summary, "data/rank1_summary.csv", row.names = FALSE)
print(rank1_summary)

###################################################################
############ FIGURE: RANK-1 SHARE #################################

pdf("figures/rank1_share.pdf", width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot(density(phi_mvn),
  xlim = c(0, 1), lwd = 2, col = "black", main = "Posterior Rank-1 Share",
  xlab = "Rank-1 share (dominant eigenvalue proportion)", ylab = "Density"
)
lines(density(phi_lat), lwd = 2, col = "grey50")
legend("topright", legend = c("MVN", "Latent"), col = c("black", "grey50"),
  lwd = 2, bty = "n")

plot(1:2, rank1_summary$mean,
  ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "Rank-1 share",
  pch = 19, cex = 1.5
)
axis(1, at = 1:2, labels = rank1_summary$model)
arrows(
  x0 = 1:2, y0 = rank1_summary$lower95, x1 = 1:2, y1 = rank1_summary$upper95,
  angle = 90, code = 3, length = 0.05
)

plot(density(1 - phi_mvn),
  xlim = c(0, 1), lwd = 2, col = "black",
  main = "Deviation from Rank-1 Structure",
  xlab = "1 - rank-1 share", ylab = "Density"
)
lines(density(1 - phi_lat), lwd = 2, col = "grey50")
legend("topright", legend = c("MVN", "Latent"), col = c("black", "grey50"),
  lwd = 2, bty = "n")
dev.off()

###################################################################
############ FIGURE: RQ3 - COVARIANCE ELLIPSES ####################

Sigma_mvn <- Reduce("+", cov_mvn) / length(cov_mvn)
Sigma_lat <- Reduce("+", cov_lat) / length(cov_lat)

pdf("figures/rq3.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot_cov_with_data(Sigma_mvn, mvn$res_p, mvn$res_c,
  main = "Unrestricted (MVN)", col_ellipse = "black")
plot_cov_with_data(Sigma_lat, lat$res_p, lat$res_c,
  main = "Latent-factor model", col_ellipse = "black")
dev.off()

cat("step 3 complete: model_comparison.csv, rank1_share.pdf, rq3.pdf\n")
