set.seed(123)
library(coda)
library(cmdstanr)

options(mc.cores = parallel::detectCores())

##################### GENERATE SAMPLE DATA ###############################################################
m <- 30                            # small areas
p <- 2                             # number of covariates
X <- cbind(1, rnorm(m))            # design matrix (intercept + one covariate)
true.beta <- c(5, 2)               # true beta
true.sigma2 <- 10               # true area-level variance
true.v <- rnorm(m, sd = sqrt(true.sigma2))
true.theta <- X %*% true.beta + true.v

D <- runif(m, 0.5, 2.0)          # known sampling variances
Y <- rnorm(m, mean = true.theta, sd = sqrt(D))  # direct estimators

test.data <- list(X = X, Y = Y, D = D, m = m, p = p, u = true.v, theta = true.theta, beta = true.beta)



##################### STAN code and function ###############################################################

# n_iter = 10000; burn_in = 1000; data = test.data; seed = 123

func.glLA <- function(n_iter, burn_in, data = data, seed = 123){
  
  model <- cmdstan_model("stan_samplers/gl_model.stan")
  
  fit <- model$sample(
    data = data,
    seed = seed,
    chains = 3,
    parallel_chains = 1,
    iter_warmup = burn_in,
    iter_sampling = n_iter
  )
  
  result <- fit$summary()
  
  theta_idx <- which(result$variable == 'theta[1]')
  
  theta_mean <- result[theta_idx:(theta_idx + data$m - 1), 2]
  theta_ci <- result[theta_idx:(theta_idx + data$m - 1), 6:7]
  
  if (max(result$rhat) <= 1.1){
    converge <- 'TRUE'
  } else {converge <- 'FALSE'}
  
  samples <- list(theta = theta_mean, coverage = theta_ci, convergence = converge)
  return(samples)
}


##################### STAN code and function ###############################################################
# results <- func.glLA(10000, 1000, data)
# 
# results
# 
# 
# df <- data.frame(true = data$theta,
#                  model = results$theta,
#                  reg = data$X %*% data$beta)
