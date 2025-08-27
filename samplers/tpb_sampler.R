library(MASS)
library(bayesplot)
library(progress)
library(statmod)
library(GIGrvg)

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
# ##################### ###############################################################
# n_iter = 10000; burn_in = 1000; data = test.data; seed = 123;index = 1;thin = 5;p=1;q=1
func.tpb <- function(n_iter, burn_in, data,  p, q, thin = 10)
{
  # parameters
  
  X <- data$X
  Y <- data$Y
  D <- data$D
  m <- nrow(X)
  np <- ncol(X)
  D.mat <- diag(D)
  
  a <- 10^-10
  b <- 10^-10
  phi <- 1
  
  # initial values
  Lambda2 = rep(1,m)
  Tau2 = 1
  Z = rep(1,m)
  Beta = rep(1,np)
  U = rep(1,m)
  
  # chain containers
  Lambda2.chain = array(0, dim=c(m,n_iter))
  Tau2.chain = rep(0, n_iter)
  Z.chain = array(0,dim=c(m,n_iter))
  Beta.chain = array(0, dim=c(np,n_iter))
  U.chain = array(0, dim=c(m,n_iter))
  
  for (index in 1:(n_iter))
  {
    if (index%%10000==0) cat(index,"\n")
    # update Tau2
    Tau2 <- 1/rgamma(1,shape=(a+m)/2,rate=sum(U^2/Lambda2)/2+b/2)
    
    # update Z
    Z <- rgamma(m, shape=p+q, rate=phi+Lambda2)
    
    # update Lambda2
    for (i in 1:m) Lambda2[i] = rgig(1,p-0.5,U[i]^2/Tau2,2*Z[i])
    
    # update Beta
    Xstd <- X/sqrt(D)
    sigma <- solve(t(Xstd)%*%Xstd)
    mean <- apply(X*(Y-U)/D,2,sum)
    mean <- sigma%*%mean
    Beta <- mvrnorm(1, mu=mean,Sigma=sigma)
    
    # update U
    sigma2 <- 1/(1/D+1/Lambda2/Tau2)
    mean <- sigma2 * (Y-X%*%Beta)/D
    U <- rnorm(m, mean=mean, sd=sqrt(sigma2))
    

    Tau2.chain[index] <- Tau2
    Lambda2.chain[,index] <- Lambda2
    Z.chain[,index] <- Z
    U.chain[,index] <- U
    Beta.chain[,index] <- Beta		
  }
  
  
  idx <- seq(burn_in, n_iter, thin) #1-5 thinn
  
  beta.thinned <- Beta.chain[, idx]
  u.thinned <- U.chain[, idx]
  tau2.thinned <- Tau2.chain[idx]
  
  theta.thinned <- X%*%beta.thinned + u.thinned
  
  ci_95 <- function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  
  theta.hat <- rowMeans(theta.thinned)
  coverage <- t(apply(theta.thinned, 1,  ci_95))
  
  samples <- list(beta = beta.thinned, theta = theta.hat, coverage = coverage)
  # samples <- list(theta = theta.hat, coverage = coverage)
  
  return(samples)
}
