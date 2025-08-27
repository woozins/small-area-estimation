library(MASS)
library(bayesplot)
library(progress)
library(statmod)


##################### GENERATE SAMPLE DATA ###############################################################
# m <- 30                            # small areas
# p <- 2                             # number of covariates
# X <- cbind(1, rnorm(m))            # design matrix (intercept + one covariate)
# true.beta <- c(5, 2)               # true beta
# true.sigma2 <- 1.0                 # true area-level variance
# true.v <- rnorm(m, sd = sqrt(true.sigma2))
# true.theta <- X %*% true.beta + true.v
# 
# D <- runif(m, 0.5, 2.0)          # known sampling variances
# Y <- rnorm(m, mean = true.theta, sd = sqrt(D))  # direct estimators
# 
# test.data <- list(X = X, Y = Y, D = D, m = m, p = p, u = true.v, theta = true.theta, beta = true.beta)
# 


# 
# ##################### Partially collapsed Gibbs sampler !!!  ###############################################################
# n_iter = 10000; burn_in = 1000; data = test.data; seed = 123;index = 1;thin = 5
func.FH <- function(n_iter, burn_in, data, thin = 10)
{
  # parameters
  X <- data$X
  Y <- data$Y
  D <- data$D
  n <- nrow(X)
  p <- ncol(X)
  D.mat <- diag(D)
  
  # initial values
  Theta <- rep(1, n)
  Beta <- rep(1, p)
  Sigma2 <- 1
  
  # chain containers
  Theta.chain <- array(0, dim=c(n, n_iter))
  Beta.chain <- array(0, dim=c(p,n_iter))
  Sigma2.chain <- rep(0, n_iter)
  
  for (index in 1:n_iter)
  {
    if (index%%10000==0) cat(index,"\n")
    
    # update Sigma2
    Sigma2 <- 1/rgamma(1,shape=n/2-1, rate=sum((Theta-X%*%Beta)^2)/2)
    
    # update Beta
    mean.Beta <- solve(t(X)%*%X)%*%t(X)%*%Theta
    var.Beta <- Sigma2*solve(t(X)%*%X)
    Beta <- mvrnorm(1, mu=mean.Beta, Sigma=var.Beta)
    
    # update Theta
    var.Theta <- 1/(1/Sigma2+1/D)
    mean.Theta <- var.Theta*(Y/D+X%*%Beta/Sigma2)
    Theta <- rnorm(n, mean=mean.Theta, sd=sqrt(var.Theta))
    

    Sigma2.chain[index] <- Sigma2
    Theta.chain[,index] <- Theta
    Beta.chain[,index] <- Beta		
  
  }
  
  idx <- seq(burn_in, n_iter, thin) #1-5 thinn
  
  beta.thinned <- Beta.chain[, idx]
  sigma2.thinned <- Sigma2.chain[idx]
  theta.thinned <- Theta.chain[, idx]
  
  ci_95 <- function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  
  theta.hat <- rowMeans(theta.thinned)
  coverage <- t(apply(theta.thinned, 1,  ci_95))
  
  samples <- list(beta = beta.thinned, sigma2 = sigma2.thinned, theta = theta.hat, coverage = coverage)
  # samples <- list(theta = theta.hat, coverage = coverage)
  
  return(samples)
  
}

