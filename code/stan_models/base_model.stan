data {
  int<lower=1> N;
  vector[N] priv;
  vector[N] pub;
  vector[N] wealth;
  array[N] int<lower=1, upper=2> sex;
  vector[N] age;
  int<lower=2> K_edu;
  array[N] int<lower=1, upper=K_edu> edu;
  vector[K_edu-1] alpha;
}

parameters {
  // Outcome-specific intercepts and slopes (same as before)
  vector[2] a;                          // [priv, pub]
  matrix[2,2] bS;                       // sex effects by outcome (row=outcome, col=sex)
  vector[2] bW;                       // wealth slopes by outcome and sex
  vector[2] bA;                         // age slopes by outcome
  vector[2] bA2;                         // age slopes by outcome
  vector[2] bE;                         // education slopes by outcome
  //matrix[2,2] bI;                         // intervention slopes by outcome
  array[2] simplex[K_edu - 1] delta_edu;        // relative step sizes between levels
  // Residual scale and correlation instead of latent T
  vector<lower=0>[2] tau;              // marginal SDs for [priv, pub]
  cholesky_factor_corr[2] L_Omega;     // Cholesky of 2x2 correlation matrix
}

transformed parameters {
  matrix[2,2] L_Sigma;
  L_Sigma = diag_pre_multiply(tau, L_Omega); // Cholesky of Sigma
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
  // Priors (match your latent model scale roughly)
  a ~ normal(0, 1);
  to_vector(bS) ~ normal(0, 1);
  delta_edu ~ dirichlet(alpha);
  bW ~ normal(0, 1);
  bA ~ normal(0, 1);
  bA2 ~ normal(0, 1);
  bE ~ normal(0, 1);
  //to_vector(bI) ~ normal(0, 1);

  tau ~ normal(0, 1);  // half-normal via <lower=0>
  L_Omega ~ lkj_corr_cholesky(2);  // weak prior on correlation

  // Likelihood
  for (i in 1:N) {
    vector[2] mu;
    vector[2] y_i;

    real base_priv = a[1]
                   + bS[1, sex[i]]
                   + bW[1] * (wealth[i])
                   + bA[1] * age[i]
                   + bA2[1] * square(age[i])
                   + edu_effect[1][edu[i]];
//                   + bI[1, sex[i]] * log1p(age[i]);

    real base_pub  = a[2]
                   + bS[2, sex[i]]
                   + bW[2] * (wealth[i])
                   + bA[2] * age[i]
                   + bA2[2] * square(age[i])
                   + edu_effect[2][edu[i]];
                   //+ bI[2, sex[i]] * log1p(age[i]);

    mu[1] = base_priv;
    mu[2] = base_pub;

    y_i[1] = priv[i];
    y_i[2] = pub[i];

    y_i ~ multi_normal_cholesky(mu, L_Sigma);
  }
}


generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    vector[2] mu;
    vector[2] y;

    real base_priv = a[1]
                   + bS[1, sex[i]]
                   + bW[1] * (wealth[i])
                   + bA[1] * age[i]
                   + bA2[1] * square(age[i])
                   + edu_effect[1][edu[i]]
;

    real base_pub  = a[2]
                   + bS[2, sex[i]]
                   + bW[2] * (wealth[i])
                   + bA[2] * age[i]
                   + bA2[2] * square(age[i])
                   + edu_effect[2][edu[i]];

    mu[1] = base_priv;
    mu[2] = base_pub;

    y[1] = priv[i];
    y[2] = pub[i];

    log_lik[i] = multi_normal_cholesky_lpdf(y | mu, L_Sigma);
  }
}
