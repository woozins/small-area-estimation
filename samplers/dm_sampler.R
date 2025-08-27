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
# 
# # 
# ##################### Partially collapsed Gibbs sampler !!!  ###############################################################
# n_iter = 10000; burn_in = 1000; data = test.data; seed = 123;index = 1;thin = 5

func.DM <- function(n_iter, burn_in, data, thin = 10)
{
  
  # parameters
  
  X <- data$X
  Y <- data$Y
  D <- data$D
  m <- nrow(X)
  np <- ncol(X)
  D.mat <- diag(D)
  
  a <- 2*mean(D)
  b <- 3
  c <- 1
  d <- 4
  
  # initial values
  V <- rep(1,m)
  Sigma2 <- 1
  Delta <- rep(1,m)
  Beta <- rep(1,np)
  p <- 0.5
  
  # chain containers
  V.chain <- array(0, dim=c(m,n_iter))
  Sigma2.chain <- rep(1, n_iter)
  Delta.chain <- array(1,dim=c(m,n_iter))
  Beta.chain <- array(0, dim=c(np, n_iter))
  p.chain <- rep(0, n_iter)
  
  for (index in 1:(n_iter))
  {
    if (index%%10000==0) cat(index,"\n")
    # update Sigma2
    
    Sigma2 <- 1/rgamma(1,shape=b+sum(Delta)/2,rate=a+sum(Delta*V*V)/2)
    
    # update p
    p <- rbeta(1, shape1=c+sum(Delta), shape2=d+m-sum(Delta))
    
    # update Delta
    p.delta <- p/(p+(1-p)*sqrt((Sigma2+D)/D)*exp(-0.5*(Y-X%*%Beta)^2*Sigma2/(D+Sigma2)/D))
    for (i in 1:m) Delta[i] <- rbinom(1,1,p.delta[i])
    
    # update Beta
    Xstd <- X/sqrt(D)
    sigma <- solve(t(Xstd)%*%Xstd)
    mean <- apply(X*(Y-Delta*V)/D,2,sum)
    mean <- sigma%*%mean
    Beta <- mvrnorm(1, mu=mean,Sigma=sigma)
    
    # update V
    mean.v <- Sigma2*(Y-X%*%Beta)/(Sigma2+D)
    var.v <- Sigma2*D/(Sigma2+D)
    V <- rnorm(m, mean=mean.v, sd=sqrt(var.v))*Delta
    
    
    
    Sigma2.chain[index] <- Sigma2
    Delta.chain[,index] <- Delta
    p.chain[index] <- p
    V.chain[,index] <- V
    Beta.chain[,index] <- Beta		
  }
  
  idx <- seq(burn_in, n_iter, thin) #1-5 thinn
  
  beta.thinned <- Beta.chain[, idx]
  sigma2.thinned <- Sigma2.chain[idx]
  v.thinned <- V.chain[,idx]
  
  theta.thinned <- X%*%beta.thinned + v.thinned 
  
  ci_95 <- function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  
  theta.hat <- rowMeans(theta.thinned)
  coverage <- t(apply(theta.thinned, 1,  ci_95))
  
  samples <- list(beta = beta.thinned, sigma2 = sigma2.thinned, theta = theta.hat, coverage = coverage)
  # samples <- list(theta = theta.hat, coverage = coverage)
  
  return(samples)
}

