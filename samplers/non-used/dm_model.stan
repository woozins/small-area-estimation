data {
  int<lower=1> m;                 // # observations
  int<lower=1> p;                 // # predictors
  vector[m] y;                    // response
  vector<lower=0>[m] D;           // known variances per obs
  matrix[m, p] X;                 // design matrix
}

transformed data {
  vector[m] sd_e = sqrt(D);
  real<lower=0> d_bar = mean(D);
  real<lower=0> d2_bar = mean(square(D));
}


parameters {
  vector[p] beta;                 // JAGS: dunif(-1e8, 1e8)
  real<lower=0, upper=1> pi;      // mixture weight (renamed from p)
  real<lower=0> inv_sigma_v2;     // JAGS: dgamma(3, 1/(2*d_bar))  [rate]
}


transformed parameters {
  vector[m] mu = X * beta;
  real<lower=0> sigma_v2 = 1.0 / inv_sigma_v2;
}

model {
  // --- Priors (explicit) ---
  for (j in 1:p)  beta[j] ~ uniform(-1e8, 1e8);
  pi ~ beta(1, 4);
  inv_sigma_v2 ~ gamma(d_bar, d2_bar);   // shape=3, rate=1/(2*d_bar)

  // --- Hierarchical prior for theta ---
  theta ~ normal()

    // --- Likelihood ---
  y ~ normal(theta, sd_e);
}

