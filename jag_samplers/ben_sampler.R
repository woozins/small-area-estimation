library(MASS)
library(bayesplot)
library(progress)

##################### define functions ###############################################################
func_g <- function(l, rho) { # checked
  return( 2*l / ((1 - rho)^2 + 2*l*rho) )
}

func_c <- function(rho){ # 
  return( (sqrt(rho)/(2*sqrt(2*pi)))*exp(-(1-rho)^2 / (2*rho)) * (pnorm(-(1-rho)/sqrt(rho)))^(-1) )
}

func_pi <- function(l, rho){ # 
  return( 2*func_c(rho)*exp(-l)/sqrt(2*l*rho + (1-rho)^2) )
}

# ll function for lambda sampling
func_ll <- function(l, sigma2, v, rho, m){ 
  return ( -m/2*log(l) + (m-1)/2*log(2*l*rho + (1-rho)^2) - l - (1-rho)^2*sum(v^2)/(4*l*sigma2)) 
}


##################### Gibbs sampler <- suppose we know rho  ###############################################################
BEN_fun <- function(n_iter, burn_in, data, rho, hyper.s=6, hyper.t=6, tune=0.005){
  
  X <- data$X
  Y <- data$Y
  n <- nrow(X)
  p <- ncol(X)
  D <- data$D
  
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
    B_en <- D / (D + sigma2[t-1]*func_g(lambda[t-1], rho))
    v[t,] <- rnorm(n, (1-B_en)*(Y-X%*%beta[t-1, ]), sd = sqrt((1-B_en)*D))
    
    #2. Beta | v, lambda ,Y, sigma2
    temp1 <- t(X)%*%solve(diag(D))%*%X
    temp2 <- t((Y-v[t,]))%*%solve(diag(D))%*%X
    beta[t,] <- mvrnorm(1, solve(temp1)%*%t(temp2), solve(temp1))
    
    #3. lambda | v, beta, Y, sigma2 <- use MH with Gibbs sampler
    # note! lambda must be positive <- automatically reject it.
    
    lambda.cand <- rnorm(1, lambda[t-1], sd = tune) 
    if (lambda.cand <= 0){ # reject
      lambda[t] <- lambda[t-1]
    } else{
      r <- exp(func_ll(lambda.cand, sigma2[t-1], v[t,], rho, n) - func_ll(lambda[t-1], sigma2[t-1], v[t,], rho, n)) 
      if(runif(1) < r){
        lambda[t] <- lambda.cand
      } else {
        lambda[t] <- lambda[t-1]
      }
    }
    
    #4. sigma2 | v, beta, y, lambda
    rate.temp <- 0.5*( ( ( ((1-rho)^2 + (2*lambda[t]*rho))*sum(v[t,]^2)) / (2*lambda[t]) ) + hyper.t )
    sigma.inv <- rgamma(1, shape = (n+hyper.s)/2,  rate = rate.temp)
    sigma2[t] <- 1/sigma.inv
    #update progress bar
    pb$tick()
  }
  
  # thinning
  idx <- seq(burn_in, n_iter, 5)
  
  v.thinned <- v[idx,]
  beta.thinned <- beta[idx,]
  sigma2.thinned <- sigma2[idx]
  lambda.thinned <- lambda[idx]
  
  theta.hat <- X%*%colMeans(beta.thinned) + colMeans(v.thinned)
  
  ci_95 <- function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  
  coverage <- t(apply(X%*%t(beta.thinned) + t(v.thinned), 1,  ci_95))  
  
  df <- cbind(v.thinned, beta.thinned, sigma2.thinned, lambda.thinned)
  
  samples <- list(v = v.thinned, lambda = lambda.thinned, beta = beta.thinned, 
                  sigma2 = sigma2.thinned, theta = theta.hat, coverage = coverage, samples_df = df)
  
  
  return(samples)
}

##################### GENERATE SAMPLE DATA ###############################################################
# set.seed(123)
# n <- 100; p <- 2
# X <- cbind(1, matrix(rnorm(n * (p - 1), mean = 10, sd = sqrt(2)), n, p - 1))
# true_beta <- c(20, 1)
# D <- rep(c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5),each=10)
# rho <- 0.5 
# true_lambda <- 2
# true_sigma2 <- 1
# true_v <- rnorm(n, 0, sd = sqrt(true_sigma2*func_g(true_lambda, rho)))
# true_e <- rnorm(n,mean=0,sd=sqrt(D))
# true_theta <- X %*% true_beta + true_v
# Y <- true_theta + true_e
# BEN_sample <- list(
#   Y = Y,
#   X = X,
#   D = D,
#   true_theta = true_theta,
#   true_v = true_v,
#   true_beta = true_beta,
#   true_e = true_e,
#   true_sigma2 = true_sigma2,
#   rho = rho,
#   true_lambda = true_lambda,
#   n=n,
#   p=p)

##################### play sampler ###############################################################
# BEN_results <- BEN_fun(n_iter = 10000, burn_in = 10, data=BEN_sample, rho = 0.5)

##################### check convergence ###############################################################
# results <- BEN_results
# plot(results[,1], type="l", ylab="v1")
# plot(results[,5], type="l", ylab="v5")
# plot(results[,10], type="l", ylab="v10")
# plot(results[,30], type="l", ylab="v30")
# plot(results[,50], type="l", ylab="v50")
# plot(results[,100], type="l", ylab="v100")
# plot(results[,101], type="l", ylab="beta_0")
# plot(results[,102], type="l", ylab="beta_1")
# plot(results[,103], type="l", ylab="sigma2")
# plot(results[,104], type="l", ylab="lambda")
# 
# colMeans(results)
# 
# results <- round(results,2)
# colnames(results)[101:102]   <- c("beta_0", "beta_1")
# 
# estimate_mode <- function(x) {
#   tb <- table(x)
#   as.numeric(names(tb)[which.max(tb)])
# }
# 
# res1 <- apply(results[,101:104],2,estimate_mode)
# res2 <- t(apply(results[,101:104],2,quantile, probs=c(0.025,0.975)))
# cbind(mode=res1,res2)
# 
# beta_hat <- round(apply(results[,101:103],2,mean),2)
# v_hat <- round(apply(results[,1:100],2,mean),2)
# 
# true_theta <- X%*%true.beta+true.v
# theta_hat <- X%*%beta_hat+v_hat
# 
# sum((true_theta-theta_hat)^2) #101
# sum(abs(true_theta-theta_hat)) #80


##################### rho selection via MCEM ###############################################################
mcem <- function(data, mc_iter, init_rho=0.5, tol = 1e-2){
  #initial setting / define functions
  rho <- init_rho
  
  #iteration
  t <- 1
  repeat{
    
    # E step <- approximate via monte carlo.
    samples <- BEN_fun(n_iter = mc_iter, burn_in = mc_iter/10, data=data, rho = rho[t])$samples_df # how to choose n_iter, burn_in?
    X = data$X
    n = nrow(X)
    p = ncol(X)
    v = samples[,1:n]
    sigma2 = samples[, n+p+1]
    l = samples[,n+p+2]
    w = rowSums(v^2)/(2*sigma2)
    func_q <- function(rho){
      return(mean(-n/2*log(sigma2*func_g(l, rho)) - w/func_g(l, rho) + log(func_pi(l, rho))))
    }
    
    # M step
    M_step <- optimize(func_q, 
                       interval = c(0,1),
                       maximum = TRUE)
    
    rho <- append(rho, M_step$maximum)
    cat("\n", t, "th iteration : rho = ", M_step$maximum, "\n", sep = "")
    if (abs(rho[t+1]-rho[t]) < tol) break
    t <- t+1
  }
  return(rho)
}
# rho <- mcem(BEN_sample, init_rho = 0.5)
# 
# rho=tail(rho,1) 
# 
# rhos.grid <- round(seq(1e-3, 0.999, length=200),4)
# Qvals    <- sapply(rhos.grid, func_q)
# 
# plot(rhos.grid, Qvals, type="l",
#      xlab="rho", ylab="Q(rho)",
#      main="MCEM Q(rho)")
# abline(v=rhos.grid[which.max(Qvals)], col="red", lty=2)

##################### fun_BEN + MCEM ###############################################################
func.ben <- function(n_iter, burn_in, data, mc_iter, tol = 1e-2, hyper.s=6, hyper.t=6, tune=0.005, init_rho = 0.5){
  rho <- tail(mcem(data, mc_iter, init_rho, tol = tol), 1)
  return(BEN_fun(n_iter, burn_in, data, rho=rho, hyper.s, hyper.t, tune))
}

