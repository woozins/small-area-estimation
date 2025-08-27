rm(list = ls())

library(MASS)
library(bayesplot)
library(progress)
library(statmod)

##################### define functions ###############################################################

func.Lambda <- function(lam , r){ # generate Lambda using lambda and r. : lambda <- only upper + diagonal elements.
  temp <- (t(lam) + lam) - diag(diag(lam))
  return(diag(rowSums(temp) + r) - temp + diag(diag(lam)))
}


n_iter = 5000

##################### Gibbs sampler <- suppose we know rho  ###############################################################
func.GLlasso.test <- function(n_iter, burn_in, X, Y, hyper.s, hyper.t, D, ga = 0.01, gb = 0.01, gr = 0.01, hr = 0.01){
  
  # define lists  
  n <- nrow(X)
  if (missing(D)) D <- diag(1, n)
  p <- ncol(X)
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
    
    # mu <- solve(diag(1, n) + (D%*%func.Lambda(lambda[,,t-1], r[t-1])/sigma2[t-1])) %*% (Y- X%*%beta[t-1, ])
    # cov <- solve(diag(1, n) + (D%*%func.Lambda(lambda[,,t-1], r[t-1]) / sigma2[t-1])) %*% D
    # 
    # v[t,] <- mvrnorm(1, mu, cov)
    
    v[t,] <- true.v
    
    
    #2. full conditional of beta
    temp1 <- t(X)%*%solve(D)%*%X
    temp2 <- t(Y-v[t,])%*%solve(D)%*%(X)
    
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
    # r[t] <- 1
    # print('r sample done')
    
    #7 full conditional of a
    
    a[t] <- rexp(1, ga + (sum(abs(v[t,]))/(2*sqrt(sigma2[t]))))
    # a[t] <- 0.03
    # print('a sample done')
    
    #8 full conditional of b ! lets use temp.v in #5.
    
    b[t] <- rexp(1, gb + (sum(temp.v)/(2*sqrt(sigma2[t]))))
    # b[t] <- 0.005
    # print('b sample done')
    
    #update progress bar
    pb$tick()
  }
  
  #1-5 thinning
  
  idx <- seq(burn_in, n_iter, 5)
  v.thinned <- v[idx,]
  beta.thinned <- beta[idx,]
  sigma2.thinned <- sigma2[idx]
  lambda.thinned <- lambda[,,idx]
  a.thinned <- a[idx]
  b.thinned <- b[idx]
  r.thinned <- r[idx]
  
  
  samples <- list(v = v.thinned, lambda = lambda.thinned, beta = beta.thinned, 
                  sigma2 = sigma2.thinned, a = a.thinned, b= b.thinned, r = r.thinned)
  return(samples)
}

##################### GENERATE SAMPLE DATA ###############################################################
set.seed(123)
# y : observed value
# X : covariate
# v : random effects
# XB + v : true theta.
# 100 country, 3 covariate(constant + 2 variable)

#observed quantities
n <- 10
p <- 3
X <- cbind(1, matrix(rnorm(n * (p - 1), mean = 10, sd = 5), n, p - 1))
D <- rep(1,n) # suppose D = 1 for all i

#hyperparameters
ga = 0.01;gb = 0.01;gr = 0.01;hr = 0.01


#parameters
true.a <- rexp(1, ga)
true.b <- rexp(1, gb)
true.r <- rgamma(1, hr, gr)

A <- matrix(0, nrow = n, ncol = n)
for(i in 1:n-1){
  for(j in i:n){
    if(j == i){
      A[i,i] <- rinvgauss(1, mean = 1/2, shape = true.a**2/2)
    } else {
    A[i,j] <- rinvgauss(1, mean = 1/2, shape = true.b**2/2)
    }
  }
}
true.Lambda <- A + t(A) - diag(diag(A))

true.sigma2 <- r
true.v <- mvrnorm(1, rep(0,n) , Sigma = true.sigma2*solve(true.Lambda))

true.beta <- c(1, 2, -1)
true.e <- rnorm(n, 0, sd = sqrt(D))

Y <- X %*% true.beta + true.v + true.e


##################### play sampler ###############################################################
results <- func.GLlasso.test(n_iter = 50000, burn_in = 1000, X = X,Y = Y, hyper.s = 6, hyper.t = 6, ga = 0.01, gb = 0.01, gr = 0.01, hr = 0.01)

##################### check convergence ###############################################################
# check v
plot(results$v[,1], type="l", ylab="v1") #v1
summary(results$v[,1])
true.v[1]
plot(results$v[,2], type="l", ylab="v2") #v2
summary(results$v[,2])
true.v[2]
plot(results$v[,50], type="l", ylab="v50") #v2
summary(results$v[,50])
true.v[50]

plot(results$beta[,1], type="l", ylab="beta0") #beta = (1,-2, 1)
summary(results$beta[,1])
plot(results$beta[,2], type="l", ylab="beta0") #beta = (1,-2, 1)
summary(results$beta[,2])
plot(results$beta[,3], type="l", ylab="beta0") #beta = (1,-2, 1)
summary(results$beta[,3])

plot(results$sigma2, type="l", ylab="sigma2") #sigma2 = 2
summary(results$sigma2)

plot(results$a, type="l", ylab="a") #0.03
summary(results$a)
plot(results$b, type="l", ylab="b") #0.005
summary(results$b)
plot(results$r, type="l", ylab="r") #1
summary(results$r)
