library(MASS)
library(bayesplot)
library(progress)
library(statmod)
##################### define functions ###############################################################

func.Lambda <- function(lam.matrix , r){ # generate Lambda using lambda and r
  return(-lam.matrix + diag(diag(lam.matrix)) + diag(rowSums(lam.matrix) + r))
}


##################### GENERATE SAMPLE DATA ###############################################################
# set.seed(123)
# # y : observed value
# # X : covariate
# # v : random effects
# # XB + v : true theta.
# # 100 country, 3 covariate(constant + 2 variable)
# 
# #observed quantities
# n <- 10
# p <- 3
# X <- cbind(1, matrix(rnorm(n * (p - 1), mean = 10, sd = 5), n, p - 1))
# D <- rep(1,n) # suppose D = 1 for all i
# 
# #hyperparameters
# ga = 0.01;gb = 0.01;gr = 0.01;hr = 0.01;hyper.s = 6; hyper.t = 6
# 
# 
# #parameters
# true.a <- rexp(1, ga)
# true.a
# true.b <- rexp(1, gb)
# true.r <- 0.1
# 
# A <- matrix(0, nrow = n, ncol = n)
# for(i in 1:n-1){
#   for(j in i:n){
#     if(j == i){
#       A[i,i] <- rinvgauss(1, mean = 1/2, shape = true.a**2/2)
#     } else {
#       A[i,j] <- rinvgauss(1, mean = 1/2, shape = true.b**2/2)
#       A[j,i] <- A[i,j]
#     }
#   }
# }
# true.Lambda <- func.Lambda(A, r = true.r)
# 
# temp <- rgamma(1, hyper.s/2, rate = hyper.t/2 )
# true.sigma2 <- 1/temp
# true.v <- mvrnorm(1, rep(0,n) , Sigma = true.sigma2*solve(true.Lambda))
# true.beta <- c(1, 2, -1)
# true.e <- rnorm(n, 0, sd = sqrt(D))
# 
# Y <- X %*% true.beta + true.v + true.e
# true.theta <- true.v + X%*%true.beta
# 
# ##################### Partially collapsed Gibbs sampler !!!  ###############################################################
# n_iter = 50000

func.GLlasso <- function(n_iter, burn_in, data, hyper.s, hyper.t, ga = 0.01, gb = 0.01, gr = 0.01, hr = 0.01, thin = 10){

  # define lists  
  X <- data$X
  Y <- data$Y
  D <- data$D
  
  n <- nrow(X)
  p <- ncol(X)
  D.mat <- diag(D)
  beta <- matrix(0, nrow = n_iter, ncol = p) #n_iter x p matrix
  sigma2 <- numeric(n_iter)
  lambda <-  array(0, dim = c(n, n, n_iter)) # note must be symmetric : only update upper triangle.
  v <- matrix(0, nrow = n_iter, ncol = n) #n_iter x n matrix 
  a <- numeric(n_iter)
  b <- numeric(n_iter)
  r <- numeric(n_iter)
  
  # initial values
  beta[1, ] <- rnorm(p, 1, 0)
  sigma2[1] <- 10
  lambda[,,1] <- matrix(rbeta(n**2, 2,2), nrow = n)
  v[1, ] <- rnorm(n, 1, sd = 3)
  a[1] <- 1
  b[1] <- 1
  r[1] <- 1
  
  pb <- progress_bar$new(
    total = n_iter,
    format = "  [:bar] :percent eta: :eta"
  )
  
  for (t in 2:n_iter) {
    #1. full conditional of v -> checked
    
    Lambda <- func.Lambda(lambda[,,t-1], r[t-1])
    
    mu <- solve(diag(1, n) + (D.mat%*%Lambda/sigma2[t-1])) %*% (Y- X%*%beta[t-1, ])
    cov <- solve(diag(1, n) + (D.mat%*%Lambda / sigma2[t-1])) %*% D.mat
        
    v[t,] <- mvrnorm(1, mu, cov)
    

    #2. full conditional of beta
    
    temp1 <- t(X)%*%solve(D.mat)%*%X
    temp2 <- t(Y-v[t,])%*%solve(D.mat)%*%(X)
     
    mu <- solve(temp1)%*%t(temp2)
    cov <- solve(temp1)
    
    beta[t,] <- mvrnorm(1, mu, cov)
    # print('beta sample done')
    
    #3. full conditional of lambda_ii

    mu <- a[t-1]*sqrt(sigma2[t-1])*(1/abs(v[t,]))
    lambda[,,t] <- diag(rinvgauss(n = n, mean = mu, shape =  a[t-1]**2))
    # print('lambda ii sample done')
    
    
    #4. full conditional of lambda_ij (i neq j)
    
    for (i in 1:n-1){
      for (j in (i+1):n){
        mu <- b[t-1]*sqrt(sigma2[t-1])*(1/abs(v[t, i] - v[t, j]))
        lambda[i,j,t] <- rinvgauss(n = 1, mean = mu, shape = b[t-1]**2)
        lambda[j,i,t] <- lambda[i,j,t]
        
      } 
    }
    
    # print('lambda ij sample done')
    
    #5. full conditional of sigma2
    
    temp.v <- matrix(0, ncol = n, nrow= n)
    for (i in 1:n-1){
      for (j in (i+1):n){
        temp.v[i, j] <- abs(v[t,i] - v[t,j])
      }
    }
    temp.v2 <- temp.v ** 2
    
    rate.temp <- 0.5*(sum(diag(lambda[,,t] + r[t-1])*(v[t,]**2)) + sum(temp.v2 * lambda[,,t]) + hyper.t)
    
    sigma.inv <- rgamma(1, shape = (n+hyper.s)/2,  rate = rate.temp)
    sigma2[t] <- 1/sigma.inv
    
    # print('sigma sample done')
    
    
    #6 full conditional of r
    
    r[t] <- rgamma(1, hr, (sum(v[t,]**2) / (2*sigma2[t])) + gr)
    # print('r sample done')
    
    #7 full conditional of a
    
    a[t] <- rexp(1, rate = ga + (sum(abs(v[t,]))/(2*sqrt(sigma2[t]))))
    # lam.inv <- sum(1/diag(lambda[,,t]))
    # a[t] <- rnorm(1, -ga/(2*lam.inv), sqrt(1/lam.inv))
    # print('a sample done')
    #8 full conditional of b ! lets use temp.v in #5.
    
    b[t] <- rexp(1, rate = gb + (sum(temp.v)/(2*sqrt(sigma2[t]))))
    # lam.inv <- (sum(1/(lambda[,,t][lambda[,,t] > 0])) - sum(1/diag(lambda[,,t])))
    # b[t] <- rnorm(1, -gb/(2*lam.inv), sqrt(1/lam.inv))
    # print('b sample done')
    
    
    
    #update progress bar
    pb$tick()
  }
  
  #1-10 thinning / outputs
  
  idx <- seq(burn_in, n_iter, 10)
  v.thinned <- v[idx,]
  beta.thinned <- beta[idx,]
  sigma2.thinned <- sigma2[idx]
  lambda.thinned <- lambda[,,idx]
  a.thinned <- a[idx]
  b.thinned <- b[idx]
  r.thinned <- r[idx]
  
  theta.hat <- X%*%colMeans(beta.thinned) + colMeans(v.thinned)
  
  ci_95 <- function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  
  coverage <- t(apply(X%*%t(beta.thinned) + t(v.thinned), 1,  ci_95))
  
  samples <- list(v = v.thinned, lambda = lambda.thinned, beta = beta.thinned,
                  sigma2 = sigma2.thinned, a = a.thinned, b= b.thinned, r = r.thinned, theta = theta.hat, coverage = coverage)
# 
  # samples <- list(theta = theta.hat, coverage = coverage)

  return(samples)
}
# # 
# ##################### play sampler ###############################################################
# results <- func.GLlasso(n_iter = 30000, burn_in = 5000, X = X,Y = Y, hyper.s = 1, hyper.t = 1,ga = 0.01, gb = 0.01, gr = 0.01, hr = 0.01)
# 
# ##################### check convergence ###############################################################
# # check v
# plot(results$v[,1], type="l", ylab="v1") #v1
# summary(results$v[,1])
# true.v[1]
# 
# plot(results$v[,2], type="l", ylab="v2") #v2
# summary(results$v[,2])
# true.v[2]
# 
# plot(results$v[,10], type="l", ylab="v50") #v2
# summary(results$v[,10])
# true.v[10]
# 
# plot(results$beta[,1], type="l", ylab="beta0") #beta = (1,-2, 1)
# summary(results$beta[,1])
# plot(results$beta[,2], type="l", ylab="beta0") #beta = (1,-2, 1)
# summary(results$beta[,2])
# plot(results$beta[,3], type="l", ylab="beta0") #beta = (1,-2, 1)
# summary(results$beta[,3])
# 
# plot(results$sigma2, type="l", ylab="sigma2") #sigma2 = 2
# summary(results$sigma2)
# true.sigma2
# 
# plot(results$a, type="l", ylab="a") #0.03
# summary(results$a)
# true.a
# 
# plot(results$b, type="l", ylab="b") #0.005
# summary(results$b)
# true.b
# 
# plot(results$r, type="l", ylab="r") #1
# summary(results$r)
# true.r
# 
# 
# true.theta <- X%*%true.beta + true.v
# results$theta
# true.theta
# 
# colMeans(results$v)
# true.v
# 
# X%*%colMeans(results$beta)
# X%*%true.beta
# # ############# toy example section
# 
# # 옵션 1) θ를 prior에서 랜덤하게 뽑기 -----------------------
# # 1) 하이퍼파라미터
# n = 5;p=2
# sigma2_0 <- 0.5
# r_0      <- 1.2
# a_0      <- 0.8
# b_0      <- 1.1
# 
# # 2) θ ~ prior
# beta_true  <- rnorm(p, 0, 3)      # β는 넓은 Gaussian prior 가정
# lam_true   <- matrix(0, n, n)
# diag(lam_true) <- statmod::rinvgauss(n,
#                                      mean  = a_0 * sqrt(sigma2_0),
#                                      shape = a_0^2)
# for(i in 1:(n-1)) for(j in (i+1):n){
#   lam_true[i,j] <- lam_true[j,i] <-
#     statmod::rinvgauss(1,
#                        mean  = b_0 * sqrt(sigma2_0) / abs(i-j),
#                        shape = b_0^2)
# }
# 
# r_true      <- r_0
# sigma2_true <- sigma2_0
# 
# # 3) v ~ p(v | λ, r, σ²)
# A_true   <- func.Lambda(lam_true, r_true)
# D        <- diag(1, n)
# Sigma_v  <- solve(diag(1,n) + (D %*% A_true)/sigma2_true) %*% D
# v_true   <- as.numeric(MASS::mvrnorm(1, rep(0,n), Sigma_v))
# 
# # 4) Y ~ p(Y | v, β, σ²)
# X         <- matrix(rnorm(n*p), n, p)
# Y         <- X %*% beta_true + v_true + rnorm(n, 0, sqrt(sigma2_true))
# 
# 
# #–– run Gibbs sampler ––#
# samples1 <- func.GLlasso(
#   n_iter   = 50000,
#   burn_in  = 10000,
#   X        = X,
#   Y        = Y,
#   hyper.s  = 1,    hyper.t = 1,
#   D        = diag(1, n),
#   ga       = 1,    gb      = 1,
#   gr       = 1,    hr      = 1
# )
# # 
# # samples1
# # sigma2_0 <- 0.5
# # r_0      <- 1.2
# # a_0      <- 0.8
# # b_0      <- 1.1
# # check v
# plot(samples1$v[,1], type="l", ylab="v1") #v1
# summary(samples1$v[,1])
# v_true[1]
# plot(samples1$v[,2], type="l", ylab="v2") #v2
# summary(samples1$v[,2])
# v_true[2]
# plot(samples1$v[,5], type="l", ylab="v50") #v2
# summary(samples1$v[,5])
# v_true[5]
# 
# 
# beta_true
# plot(samples1$beta[,1], type="l", ylab="beta0") #beta = (1,-2, 1)
# summary(samples1$beta[,1])
# plot(samples1$beta[,2], type="l", ylab="beta0") #beta = (1,-2, 1)
# summary(samples1$beta[,2])
# 
# plot(samples1$sigma2, type="l", ylab="sigma2") #sigma2 = 2
# summary(samples1$sigma2)
# # sigma2_0 <- 0.5
# 
# 
# 
# plot(samples1$a, type="l", ylab="a") #0.03
# summary(samples1$a)
# # a_0      <- 0.8
# 
# plot(samples1$b, type="l", ylab="b") #0.005
# summary(samples1$b)
# # b_0      <- 1.1
# 
# plot(samples1$r, type="l", ylab="r") #1
# summary(samples1$r)
# # r_0      <- 1.2




