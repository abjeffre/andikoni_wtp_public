###################################################################
############ 2. PAPER FIGURES #####################################
# Reads  data/wtp_data.csv, data/posterior_base.RDS,
#        data/posterior_latent.RDS
# Writes figures/base.pdf, figures/predictors.pdf, figures/latent.pdf
#
# figures/rq3.pdf is written by step 3, which needs the residual covariances.

library(rethinking)
library(posterior)
source("code/functions/utility.R")
source("code/prepare_data.R")

post_base <- readRDS("data/posterior_base.RDS")
post_latent <- readRDS("data/posterior_latent.RDS")

# Panel helper: posterior mean with an 89% interval, one row per coefficient.
coef_panel <- function(draws, main, xlim = c(-.9, .9)) {
  draws <- rev(draws)
  mu <- colMeans(draws)
  errors <- apply(draws, 2, PI, .89)
  my_dotchart(mu, cex = .9, xlab = "Posterior estimate", pch = 16,
    xlim = xlim, main = main)
  for (i in 1:ncol(draws)) {
    segments(
      x0 = errors[1, i], x1 = errors[2, i],
      y0 = i, y1 = i, lwd = 2, col = "black"
    )
  }
  abline(v = 0, col = "black", lty = 2)
}

###################################################################
############ FIGURE: BASE RELATIONSHIP ############################

pdf("figures/base.pdf", width = 8, height = 4)
seq <- c(1, 5, 10, 25, 50, 100, 200, 300, 400)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(jitter(log10(df$priv_price / 2500), 2) ~ jitter(log10(df$pub_price / 2500), 2),
  xlab = "log(Community WTP)", ylab = "log(Private WTP)", pch = 16,
  yaxt = "n", xaxt = "n",
  xlim = c(0, log10(max(seq))), ylim = c(0, log10(max(seq))),
  col = col.alpha("black", .75)
)
axis(2, at = log10(seq), labels = seq)
axis(1, at = log10(seq), labels = seq)
lines(c(-10, 100), c(-10, 100), lty = 2)

dens2(post_latent$a[, 1] - post_latent$a[, 2],
  col_fill = col.alpha("grey9", .2), show.HPDI = .95, show.zero = TRUE,
  xlab = "Intercept Contrast (Comm - Priv)"
)
dev.off()

###################################################################
############ FIGURE: PREDICTORS ###################################

pdf("figures/predictors.pdf", width = 10)
post <- post_base
par(mfrow = c(1, 3), mar = c(5, 6, 3, 1))

coef_panel(data.frame(
  Female = post$bS[, 1, 1] - post$bS[, 1, 2],
  Age = post$bA[, 1],
  "Age^2" = post$bA2[, 1],
  Edu = post$bE[, 1],
  Wealth = post$bW[, 1]
), main = "Private")

coef_panel(data.frame(
  Female = post$bS[, 2, 1] - post$bS[, 2, 2],
  Age = post$bA[, 2],
  "Age^2" = post$bA2[, 2],
  Edu = post$bE[, 2],
  Wealth = post$bW[, 2]
), main = "Community")

# Every row is a community-minus-private contrast of the same coefficient.
# (The original script indexed post$bE[2] and post$bW[2] without a comma on the
# last two rows, recycling one draw instead of taking the second outcome's
# column; corrected here.)
coef_panel(data.frame(
  Female = post$bS[, 2, 1] - post$bS[, 2, 2] - (post$bS[, 1, 1] - post$bS[, 1, 2]),
  Age = post$bA[, 2] - post$bA[, 1],
  "Age^2" = post$bA2[, 2] - post$bA2[, 1],
  Edu = post$bE[, 2] - post$bE[, 1],
  Wealth = post$bW[, 2] - post$bW[, 1]
), main = "Contrasts (Comm - Priv)")
dev.off()

###################################################################
############ FIGURE: LATENT AXIS ##################################

pdf("figures/latent.pdf")
dens2(post_latent$lambda,
  col_fill = col.alpha("dodgerblue", .2), show.HPDI = .95, show.zero = TRUE,
  xlab = "Latent Pro-Social - Selfish Dimension (Lambda)"
)
dev.off()

cat("step 2 complete: base.pdf, predictors.pdf, latent.pdf\n")
