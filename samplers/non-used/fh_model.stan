data {
  int<lower=1> m;
  int<lower=1> p;
  vector[m] Y;                 // sample means
  vector<lower=0>[m] D;         // known sampling variances
  matrix[m, p] X;               // design matrix
}

transformed data {
  vector[m] sd_e = sqrt(D);
}

parameters {
  vector[p] beta;                       // JAGS: dunif(-1e8, 1e8)
  real<lower=0> sigma2_fh;              // JAGS: dunif(0, 1e8)
  vector[m] theta;                      // latent true values
}

transformed parameters {
  vector[m] mu = X * beta;
}

model {
  // --- Explicit priors (JAGS와 동일 의미) ---
  for (j in 1:p)  beta[j]   ~ uniform(-1e8, 1e8);
  sigma2_fh       ~ uniform(0, 1e8);

  // --- Hierarchical prior for theta ---
  
  theta ~ normal(mu, sqrt(sigma2_fh));

  // --- Likelihood ---
  Y ~ normal(theta, sd_e);
}
