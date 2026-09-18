data {
  int<lower=1> N;
  vector[N] priv;                       // outcome 1: private planting
  vector[N] pub;                        // outcome 2: public planting

  // outcome-side covariates (already on your chosen scale)
  vector[N] wealth;
  array[N] int<lower=1, upper=2> sex;   // 1 or 2
  vector[N] age;
  int<lower=2> K_edu;
  array[N] int<lower=1, upper=K_edu> edu;
  vector[K_edu-1] alpha;
}

parameters {
  // residual SDs per outcome
  vector<lower=0>[2] sigma;             // [priv, pub]

  // outcome-specific intercepts and slopes
  vector[2] a;                          // intercepts [priv, pub]
  matrix[2,2] bS;                       // sex effects by outcome (row=outcome, col=sex)
  vector[2] bW;                       // wealth slopes by outcome and sex
  vector[2] bA;                         // age slopes by outcome
  vector[2] bA2;                         // age slopes by outcome
  vector[2] bE;                         // education slopes by outcome
  //matrix[2,2] bI;                         // intervention slopes by outcome
  array[2] simplex[K_edu - 1] delta_edu;        // relative step sizes between levels
  // latent prosocial orientation: shared random effect
  vector[N] T;                          // 1D latent factor per person
  //real<lower=0> lambda;                         // log(lambda_pub)
  real<lower =0> lambda;                         // log(lambda_pub)
}

transformed parameters{
  array[2] vector[K_edu] edu_effect;            // effect for each level 1..K
  array[2] vector[K_edu - 1] inc;
  // rescale simplex to K-1 total so average step ~ beta_edu
  for(i in 1:2){
    inc[i] = (K_edu - 1) * delta_edu[i];
    edu_effect[i][1] = 0;
    for (k in 2:K_edu) {
      edu_effect[i][k] = edu_effect[i][k - 1] + bE[i] * inc[i][k - 1];
      }
  }
}


model {
  vector[N] mu_priv;
  vector[N] mu_pub;

  // ---- Priors ----
  a ~ normal(0, 1);
  to_vector(bS) ~ normal(0, 1);
  delta_edu ~ dirichlet(alpha);
  bW ~ normal(0, 1);
  bA ~ normal(0, 1);
  bA2 ~ normal(0, 1);
  bE ~ normal(0, 1);
  //to_vector(bI) ~ normal(0, 1);

  lambda  ~ normal(0, 0.5);            // moderately tight to stabilize loadings
  sigma    ~ normal(0, 1);              // half-normal via <lower=0>

  // latent trait as simple random effect, scale anchored at 1
  T ~ normal(0, 1);

  // ---- Linear predictors ----
  for (i in 1:N) {
    real base_priv = a[1]
                   + bS[1, sex[i]]
                   + bW[1] * (wealth[i])
                   + bA[1] * age[i]
                   + bA2[1] * square(age[i])
                   + edu_effect[1][edu[i]];
                   //+ bI[1, sex[i]] * log1p(age[i]);

    real base_pub  = a[2]
                   + bS[2, sex[i]]
                   + bW[2] * (wealth[i])
                   + bA[2] * age[i]
                   + bA2[2] * square(age[i])
                   + edu_effect[2][edu[i]];
                   //+ bI[2, sex[i]] * log1p(age[i]);

    // SAME latent, opposite loadings: tradeoff structure
    mu_priv[i] = base_priv - lambda * T[i];
    mu_pub[i]  = base_pub  + lambda  * T[i];
  }

  // ---- Likelihood ----
  priv ~ normal(mu_priv, sigma[1]);
  pub  ~ normal(mu_pub,  sigma[2]);
}

generated quantities {
  // Fitted means
  vector[N] mu_priv_hat;
  vector[N] mu_pub_hat;

  // Fitted gap per person and population-average gap
  vector[N] gap_mu;
  real delta_pop;

  // Pointwise log-likelihood for WAIC/LOO (joint per person)
  vector[N] log_lik;

  for (i in 1:N) {
    real base_priv = a[1]
                   + bS[1, sex[i]]
                   + bW[1] * (wealth[i])
                   + bA[1] * age[i]
                   + bA2[1] * square(age[i])
                   + edu_effect[1][edu[i]];
                   //+ bI[1, sex[i]] * log1p(age[i]);

    real base_pub  = a[2]
                   + bS[2, sex[i]]
                   + bW[2] * (wealth[i])
                   + bA[2] * age[i]
                   + bA2[2] * square(age[i])
                   + edu_effect[2][edu[i]];
                  // + bI[2, sex[i]] * log1p(age[i]);

    mu_priv_hat[i] = base_priv - lambda * T[i];
    mu_pub_hat[i]  = base_pub  + lambda  * T[i];

    gap_mu[i] = mu_pub_hat[i] - mu_priv_hat[i];

    // joint log-likelihood contribution of (priv_i, pub_i)
    log_lik[i] =
      normal_lpdf(priv[i] | mu_priv_hat[i], sigma[1]) +
      normal_lpdf(pub[i]  | mu_pub_hat[i],  sigma[2]);
  }
}
