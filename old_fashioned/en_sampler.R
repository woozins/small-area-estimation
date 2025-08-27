rm(list = ls())

library(MASS)
library(bayesplot)
library(progress)

##################### define functions ###############################################################
func.g <- function(l, rho) { # checked
  return (2*l / ((1 - rho)**2 + 2*l*rho))
}

func.c <- function(rho){ # 
  return((sqrt(rho)/(2*sqrt(2*pi)))*exp(-(1-rho)**2 / (2*rho))*(pnorm(-(1-rho)/sqrt(rho)))**(-1))
}

func.pi <- function(l, rho){ # 
  return(2*func.c(rho)*exp(-l)/sqrt(2*l*rho + (1-rho)**2))
}

# ll function for lambda sampling
func.ll <- function(l, sigma2, vs, rho, m){ 
  return ((-m/2)*log(l) + (m/2 - 1/2)*log(2*l*rho + (1-rho)**2) - l - (((1-rho)**2)*sum(vs**2)/(4*l*sigma2)))
}


##################### Gibbs sampler <- suppose we know rho  ###############################################################
func.gibbs <- function(n_iter, burn_in, X, Y, hyper.s, hyper.t, tune, rho){
  
  n <- nrow(X)
  p <- ncol(X)
  #list : row <- iterations
  beta <- matrix(0, nrow = n_iter, ncol = p)
  sigma2 <- numeric(n_iter)
  lambda <- numeric(n_iter)
  v <- matrix(0, nrow = n_iter, ncol = n) 
  
  # hyperparameter for prior
  # hyper.s <- 6; hyper.t <- 6   # sigma2 ~ Inv-Gamma(s/2, t/2)
  
  #tuning parameter
  # tune <- 0.0005
  
  # initial values
  beta[1, ] <- 2
  sigma2[1] <- 10
  lambda[1] <- 2
  v[1, ] <- rnorm(n, 1, sd= 3)  
  
  pb <- progress_bar$new(
    total = n_iter,
    format = "  [:bar] :percent eta: :eta"
  )
  
  
  for (t in 2:n_iter) {
    #1. v | Y, beta, lambda, sigma2
    B_en <- D / (D + sigma2[t-1]*func.g(lambda[t-1], rho))
    v[t,] <- rnorm(n, (1-B_en)*(Y-(X%*%beta[t-1, ])), sd = sqrt((1-B_en)*D))
    
    #2. Beta | v, lambda ,Y, sigma2
    temp1 <- t(X)%*%solve(diag(D))%*%X
    temp2 <- t((Y-v[t,])*(1/D))%*%(X)
    beta[t,] <- mvrnorm(1, solve(temp1)%*%t(temp2), solve(temp1))
    
    #3. lambda | v, beta, Y, sigma2 <- use MH with Gibbs sampler
    # note! lambda must be positive <- automatically reject it.
    
    lambda.cand <- rnorm(1, lambda[t-1], sd = tune) 
    if (lambda.cand <= 0){ # reject
      lambda[t] <- lambda[t-1]
    } else{
    r <- exp(func.ll(lambda.cand, sigma2[t-1], v[t,], rho, n) - func.ll(lambda[t-1], sigma2[t-1], v[t,], rho, n)) 
    if(runif(1) < r){
      lambda[t] <- lambda.cand
    } else {
      lambda[t] <- lambda[t-1]
    }
    }
    
    #4. sigma2 | v, beta, y, lambda
    rate.temp <- 0.5*(((((1-rho)**2 + (2*lambda[t]*rho))*sum(v[t,]**2))/(2*lambda[t])) + hyper.t)
    sigma.inv <- rgamma(1, shape = (n+hyper.s)/2,  rate = rate.temp)
    sigma2[t] <- 1/sigma.inv
    #update progress bar
    pb$tick()
  }
  
  #1-10 thinning
  
  idx <- seq(burn_in, n_iter, 10)
  v.thinned <- v[idx,]
  beta.thinned <- beta[idx,]
  sigma2.thinned <- sigma2[idx]
  lambda.thinned <- lambda[idx]
  
  samples <- cbind(v.thinned, beta.thinned, sigma2.thinned, lambda.thinned)
  return(samples)
  }

##################### GENERATE SAMPLE DATA ###############################################################
set.seed(123)

# y : observed value
# X : covariate
# v : random effects
# XB + v : true theta.
# 100 country, 3 covariate(constant + 2 variable)

n <- 5
p <- 3
rho <- 0.5 # suppose set
true.lambda <- 2
X <- cbind(1, matrix(rnorm(n * (p - 1), mean = 10, sd = 5), n, p - 1))
D <- rep(1,n) # suppose D = 1 for all i
true.beta <- c(1, 2, -1)
true.sigma2 <- 8
true.v <- rnorm(n, 0, sd = sqrt(true.sigma2*func.g(true.lambda, rho)))
true.e <- rnorm(n, 0, sd = sqrt(D))
Y <- X %*% true.beta + true.v + true.e

##################### play sampler ###############################################################
results <- func.gibbs(n_iter = 10000, burn_in = 1000, X = X,Y = Y, hyper.s = 1, hyper.t = 1, tune =0.005, rho = 0.5)
##################### check convergence ###############################################################
burn_in = 1000
dim(results)
plot(results[,n+1], type="l", ylab="beta_0")
plot(results[,n+2], type="l", ylab="beta_1")
plot(results[,n+3], type="l", ylab="beta_2")
plot(results[,n+4], type="l", ylab="sigma2")
plot(results[,n+5], type="l", ylab="lambda")

##################### rho selection via MCEM ###############################################################

# NOTICE 
# CONCERN rho converging to local minima

## Empirical bayes with gibbs sampling(MCEM)

#initial setting / define functions
init.rho <- 0.5
rhos <- c(0)
rhos <- append(rhos, init.rho)
tol <- 10e-4



# calculating q estimate
func.q <- function(rho, m, sigma2, l, v){
  return(-mean(-((m/2)*log(sigma2*func.g(l, rho))) - (sum(v**2)/(2*sigma2*func.g(l, rho))) + log(func.pi(l, rho))))
}


#iteration
t <- 2
while (abs(rhos[t] - rhos[t-1]) > tol){
  
  # E step <- approximate via monte carlo.
  samples <- func.gibbs(n_iter = 10000, burn_in = 1000, X = X, Y = Y, hyper.s = 6, hyper.t = 6, tune = 0.005, rho = rhos[t])
  
  # M step
  M_step <- optimize(function(x) func.q(rho, ncol(X), sigma2 = samples[, nrow(X)+ncol(X)+1],
                                        l = samples[,nrow(X)+ncol(X)+2], 
                                        v = samples[,1:nrow(X)]), 
                                        interval = c(0.001,1e+15))
  
  rhos <- append(rhos, M_step$minimum)
  
  print(paste0(t-1, 'th iteration : rho = ', M_step$minimum))
  
  t <- t+1
}


