data { int<lower=1> m; // 소지역 수 
  int<lower=1> p; // 공변량 수 
  vector[m] Y; // direct estimator
  matrix[m, p] X; // area-level covariates 
  vector<lower=0>[m] D; // known sampling variances
}

parameters { // random area effects 
  vector[p] beta; // regression coefficients 
  vector<lower=0>[m] lambda; // local scales ~ laplacian 
  real<lower=0> tau; // global scale ~ half-Cauchy(0,1)
  vector[m] z;
}

transformed parameters {
  vector[m] u;
  u = z .* (lambda * tau);
}

model { 
  // priors
  lambda ~ double_exponential(0, 1); // laplacian 
  tau ~ cauchy(0, 1); // half-Cauchy via <lower=0> 
  beta ~ uniform(-1e+8, 1e+8);
  // hierarchical prior on beta 
  z ~ normal(0, 1); 
  
  // likelihood (hierarchical form)
  for (i in 1:m) {
    Y[i] ~ normal( dot_product(row(X, i), beta) + u[i], sqrt(D[i])); 
  }
}

generated quantities {
  vector[m] theta; // area means
  for (i in 1:m)
  theta[i] = dot_product(row(X, i), beta) + u[i]; 
}  
