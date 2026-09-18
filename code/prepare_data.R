###################################################################
############ BUILD THE STAN INPUT LIST ############################
# Reads  data/wtp_data.csv
# Defines stan_data (and df, the analysis frame) in the calling environment.
# Sourced by 1_MAKE_fit_models.R, 2_MAKE_figures.R and 3_MAKE_model_comparison.R
# so that all three agree on exactly one construction of the model inputs.

df <- read.csv("data/wtp_data.csv", stringsAsFactors = FALSE)

if (nrow(df) != 71) {
  stop("expected 71 respondents in data/wtp_data.csv, found ", nrow(df))
}
if (anyNA(df)) {
  stop("data/wtp_data.csv contains missing values")
}

# Outcomes are bids in TSh; log1p keeps the zero bids, which are meaningful
# (a respondent who would give up no trees on community land).
df$priv <- log1p(df$priv_price)
df$pub <- log1p(df$pub_price)

edu_index <- as.numeric(as.factor(df$ed_rev))

stan_data <- list(
  N      = nrow(df),
  priv   = as.vector(df$priv),
  pub    = as.vector(df$pub),
  wealth = rethinking::standardize(log1p(df$wealth)),
  sex    = as.integer(as.numeric(as.factor(df$sex))),
  age    = rethinking::standardize(df$age),
  edu    = edu_index,
  K_edu  = max(edu_index),
  alpha  = rep(2, max(edu_index) - 1)
)

# The residual helpers in utility.R read $priv, $pub, $wealth and $age off this
# list, so they must be the standardized versions the models actually saw.
