data {
  int<lower=1> N;
  vector[N] priv;                       // outcome 1: private planting
  vector[N] pub;                        // outcome 2: public planting

  // covariates (demeaned continuous variables strongly recommended)
  vector[N] wealth;
  array[N] int<lower=1, upper=2> sex;
  vector[N] age;
  int<lower=2> K_edu;
  array[N] int<lower=1, upper=K_edu> edu;
  vector[K_edu-1] alpha;                // dirichlet prior hyperparameter
}

parameters {
  // residual SDs per outcome
  vector<lower=0>[2] sigma;             // [priv, pub]

  // outcome-specific intercepts only — covariates no longer enter here
  vector[2] a;

  // loading: single parameter governing trade-off magnitude
  real<lower=0> lambda;

  // latent factor regression coefficients (on the T scale)
  real gW;                              // wealth effect on T
  real gA;                              // age (linear) effect on T
  real gA2;                             // age (quadratic) effect on T
  real gE;                              // education slope on T
  vector[2] gS;                         // sex effects on T (indexed 1,2)

  // ordinal education: monotone steps on T scale
  simplex[K_edu - 1] delta_edu;

  // latent factor residuals — mean zero, scale 1
  // T_i = mu_T_i + eta_i, eta_i ~ Normal(0,1)
  // We sample the raw residuals and construct T in transformed parameters
  vector[N] eta;                        // individual residuals, ~ Normal(0,1)
}

transformed parameters {
  // Education cumulative effects on latent scale
  vector[K_edu] edu_effect_T;
  vector[K_edu - 1] inc;

  inc = (K_edu - 1) * delta_edu;
  edu_effect_T[1] = 0;
  for (k in 2:K_edu) {
    edu_effect_T[k] = edu_effect_T[k-1] + gE * inc[k-1];
  }

  // Latent factor: covariate-driven mean + individual residual
  // Identification notes:
  //   (1) Continuous covariates should be demeaned in data block so that
  //       E[mu_T] ~ 0 at covariate means, avoiding conflict with a[1], a[2]
  //   (2) eta ~ Normal(0,1) anchors the residual scale of T to 1
  //   (3) lambda is unrestricted in sign; the sign convention is:
  //       higher T  =>  lower priv WTP, higher pub WTP (pro-community tilt)
  vector[N] mu_T;
  vector[N] T;

  for (i in 1:N) {
    mu_T[i] = gS[sex[i]]
            + gW * wealth[i]
            + gA * age[i]
            + gA2 * square(age[i])
            + edu_effect_T[edu[i]];
    T[i] = mu_T[i] + eta[i];
  }
}

model {
  // ---- Priors ----

  // Outcome intercepts: relatively diffuse
  a ~ normal(0, 1);

  // Loading: moderately tight — without outcome-specific covariate effects
  // the model is more dependent on lambda carrying the right magnitude,
  // so we keep this tighter than a flat prior to aid geometry
  lambda ~ normal(0, 0.5);

  // Outcome residual SDs
  sigma ~ normal(0, 1);                 // half-normal via <lower=0>

  // Latent factor residuals — anchors scale of T
  eta ~ normal(0, 1);

  // Covariate effects on latent factor
  gS ~ normal(0, 1);
  gW ~ normal(0, 1);
  gA ~ normal(0, 1);
  gA2 ~ normal(0, 1);
  gE ~ normal(0, 1);

  // Ordinal education steps
  delta_edu ~ dirichlet(alpha);

  // ---- Likelihood ----
  // All covariate influence flows through T; outcomes get intercepts only
  priv ~ normal(a[1] - lambda * T, sigma[1]);
  pub  ~ normal(a[2] + lambda * T, sigma[2]);
}

generated quantities {
  vector[N] mu_priv_hat;
  vector[N] mu_pub_hat;
  vector[N] gap_mu;
  real delta_pop;
  vector[N] log_lik;

  for (i in 1:N) {
    mu_priv_hat[i] = a[1] - lambda * T[i];
    mu_pub_hat[i]  = a[2] + lambda * T[i];
    gap_mu[i]      = mu_pub_hat[i] - mu_priv_hat[i];

    log_lik[i] =
      normal_lpdf(priv[i] | mu_priv_hat[i], sigma[1]) +
      normal_lpdf(pub[i]  | mu_pub_hat[i],  sigma[2]);
  }

  // Population-average gap (averages out individual T variation)
  delta_pop = mean(gap_mu);
}
