set.seed(123)
library(coda)
# library(rjags)

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

fh_model_string <- "
model {
  # 1) 관측식: 표본평균 yi와 오차분산 vi
  for (i in 1:m) {
    yi[i] ~ dnorm(theta[i], prec.e[i])
    prec.e[i] <- 1/D[i]
  }

  # 2) 계층식: 잠재진실 theta[i]
  for (i in 1:m) {
    theta[i] ~ dnorm(mu[i], 1/sigma2.fh)
    mu[i] <- inprod(X[i, ], beta[])
  }

  # 3) 회귀계수 베타에 비정보적 사전
  for (j in 1:p) {
    beta[j] ~ dunif(-1e+8, 1e+8)
  }

  # 4) 영역간 분산(σ_u^2)의 사전: 감마→정밀도 τ_u
  sigma2.fh ~ dunif(0, 1e+8)
}
"
#모델 파일로 저장
writeLines(fh_model_string, con = "samplers/fh_model.jags")
##################### JAGS code and function ###############################################################
func.FH <- function(n_iter, burn_in , data){
    
  data_list <- list(
    yi = as.numeric(data$Y),
    D = data$D,
    X  = data$X,
    m  = length(data$Y),
    p  = ncol(data$X)
  )
  
  # 1) 모델 파일로 저장
  
  # 2) 초기값 함수 (체인별로 약간씩 다르게 주면 수렴 진단에 유리)
  inits <- function() {
    list(
      beta = rep(0, data_list$p),
      sigma2.fh = 1
    )
  }
  
  # 3) JAGS 모델 객체 생성
  jags_mod <- jags.model(
    file    = "samplers/fh_model.jags",
    data    = data_list,
    inits   = inits,
    n.chains= 3,
    quiet   = TRUE
  )
  
  # 4) 예열 (버닝) 단계
  update(jags_mod, n.iter = burn_in)
  
  # 5) 사후표본 추출
  params <- c("beta", "sigma2.fh", "theta")
  coda_samples <- coda.samples(
    model      = jags_mod,
    variable.names = params,
    n.iter     = n_iter,
    thin       = 5
  )
  
  
  
  # 1) 요약 통계
  # summary(coda_samples)
  # 
  # # 2) 눈으로 확인: 트레이스플롯, 밀도플롯
  # plot(coda_samples)
  # 
  # # 3) Gelman–Rubin 진단
  # gelman.diag(coda_samples)
  
  # 4) 소지역별 사후 추정값
  theta_post <- do.call(rbind, lapply(coda_samples, function(x) as.matrix(x)[, grep("^theta\\[", colnames(x))]))
  theta_mean <- apply(theta_post, 2, mean)
  theta_ci   <- apply(theta_post, 2, quantile, probs = c(0.025, 0.975))

  samples <- list(theta = theta_mean, coverage = t(theta_ci))
  return(samples)
}

