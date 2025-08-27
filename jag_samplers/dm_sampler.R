set.seed(123)
library(coda)
library(rjags)

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
# data <- list(X = X, Y = Y, D = D, u = true.v, theta = true.theta, beta = true.beta)
##################### JAGS scripts ###############################################################

dm_model_string <- "
model {
  # 1) 관측식
  for (i in 1:m) {
    y[i] ~ dnorm(theta[i], prec.D[i])
    prec.D[i] <- 1 / D[i]
    mu[i] <- inprod(X[i,], beta[])  
    theta[i] <- mu[i] + delta[i] * v[i]
    
    # 2) 혼합 랜덤효과: 지표 S, 슬랩 v
    delta[i] ~ dbern(p)
    v[i] ~ dnorm(0, prec.v)             
  }

  prec.v <- inv_sigma_v2

  # 3) 사전분포 (proper priors)
  for (j in 1:q) {
    beta[j] ~ dunif(-1e+8, 1e+8)         # 약한 정규 사전 (또는 dflat())
  }
  p ~ dbeta(1, 4)          # c,d > 0
  inv_sigma_v2 ~ dgamma(3, 1/(2*d_bar))     # a,b > 0  (Inv-Gamma on sigma_v^2)

  # 4) 파생량
  sigma_v2 <- 1.0 / inv_sigma_v2
}

"
#모델 파일로 저장
writeLines(dm_model_string, con = "samplers/dm_model.jags")
##################### JAGS code and function ###############################################################
func.DM <- function(n_iter, burn_in , data){
  
  data_list <- list(
    y = as.numeric(data$Y),
    D = data$D,
    d_bar = mean(data$D),
    X  = data$X,
    m  = length(data$Y),
    q  = ncol(data$X)
  )
  
  # 1) 모델 파일로 저장
  
  # 2) 초기값 함수 (체인별로 약간씩 다르게 주면 수렴 진단에 유리)
  inits <- function() {
    list(
      beta = rep(0, data_list$q)
      )
  }
  
  # 3) JAGS 모델 객체 생성
  jags_mod <- jags.model(
    file    = "samplers/dm_model.jags",
    data    = data_list,
    inits   = inits,
    n.chains= 3,
    quiet   = TRUE
  )
  
  # 4) 예열 (버닝) 단계
  update(jags_mod, n.iter = burn_in)
  
  # 5) 사후표본 추출
  params <- c("beta", "sigma_v2", "theta")
  coda_samples <- coda.samples(
    model      = jags_mod,
    variable.names = params,
    n.iter     = n_iter,
    thin       = 5
  )
  
  
  # 4) 소지역별 사후 추정값
  theta_post <- do.call(rbind, lapply(coda_samples, function(x) as.matrix(x)[, grep("^theta\\[", colnames(x))]))
  theta_mean <- apply(theta_post, 2, mean)
  theta_ci   <- apply(theta_post, 2, quantile, probs = c(0.025, 0.975))
  
  samples <- list(theta = theta_mean, coverage = t(theta_ci), data = coda_samples)
  return(samples)
}

# results <- func.DM(50000, 1000, data)

############### checking convergence #############################################################################
# # 1) 요약 통계
# summary(results$data)
# # 
# # # 2) 눈으로 확인: 트레이스플롯, 밀도플롯
# plot(results$data)
# # 
# # # 3) Gelman–Rubin 진단
# gelman.diag(results$data)
# 
# data.frame(pred = results$theta,
#            true = true.theta)
########################################################################################################################
