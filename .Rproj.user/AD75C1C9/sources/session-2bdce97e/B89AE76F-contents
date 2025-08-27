rm(list = ls())


# library(kableExtra)
library(dplyr)

# for parallel computing
library(future)
library(furrr)

# for mcmc
source('samplers/fh_sampler.R')
source('samplers/ben_sampler.R')
source('samplers/lasso_sampler.R')
source('samplers/dm_sampler.R')
source('samplers/tpb_sampler.R')
source('samplers/ng_sampler.R')


source('evalftns.R')

library(coda)
library(cmdstanr)

# for progress bar
library(progressr)

# 오버서브스크립션 방지
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

###########################Initial setting ########################### ########################### ########################### 

m.list <- c(1000, 500, 100) # Note it is reverse!
V.list <- seq(0.5, 5, 0.5)
re.list <- c('FH', 'DM', 'a', 'b','c')
models <- c('GLlasso', 'BEN', 'FH', 'DM', 'tpb' ,'ng')
eval.list <- c('AAD','ASD','ARB','ASRB','Coverage')
param.grid <- data.frame(expand.grid(m.list, re.list))
colnames(param.grid) <- c('m', 're')

########################### data generation. ########################### ########################### ########################### 
sim.data <- function(m, re , D){
  X <- cbind(rep(1, m), rnorm(m, 10, sd = sqrt(2)))
  beta <- c(20, 1)
  if (re == 'FH'){
    u <- rnorm(m, 0, 1)   
  } else if (re == 'DM') {
    delta <- rbinom(m, 1, 0.2)
    v <- rnorm(m, 0, 5)
    u <- delta*v
  } else if (re == 'a') {
    vec <- cbind(rep(c(1,1,1,1,0), m/5), rep(1-c(1,1,1,1,0), m/5)) 
    v <- cbind(rnorm(m, 0, 5),rnorm(m, 0, 1))
    u <- rowSums(v*vec)
  } else if (re == 'b') {
    sigma.list <- seq(0.5, 5, 0.5)
    sigma2 <- sample(sigma.list, 1)
    u <- rnorm(m, 0, sqrt(sigma2))
  } else if (re == 'c') {
    u <- rt(m, df = 3)  
  }

  e <- rnorm(m, 0, sd = sqrt(diag(D)))
  Y <- X%*%beta + u + e
  theta <- X%*%beta + u
  result <- list(Y = as.vector(Y), X = X, m = nrow(X), p = ncol(X), theta = theta, beta = beta, u = u, D = D)
  return(result)
}

######## Model Fitting & evaluation ################################################################################################
#Compare 5 models : FH model, DM model, GL model, our EN, LASSO model.

#setting
n_iter <- 10000
burn_in <- 1000

n.models <- length(models) # 테스트 모델 수
n.ms <- 3                  # 
n.setting <- 5             # 실험 셋팅 수
n.datas <- 10              # 각 셋팅 별 데이터 생성 수


# Define result array
result <- array(NA, dim = c(length(m.list), length(re.list), length(eval.list),  n.datas, n.models)) # m x re x eval x each dataset x models  
dimnames(result) <- list(
  m = paste0("m", m.list),
  re = paste0("re", re.list),
  eval = paste0('eval', eval.list),
  rep = paste0("rep", 1:n.datas),
  model = paste0("model", models)
)


# 작업목록
tasks <- tidyr::crossing(
  i    = seq_len(nrow(param.grid)),
  iter = seq_len(n.datas)
)


# 병렬 작업
plan(multisession, workers = 15)

## Progress Bar style
handlers(global = TRUE)

library(progressr)

handlers(handler_progress(
  format = ":current/:total [:bar] :percent eta: :eta",
  clear = FALSE,
  width = 150
))

#병렬처리 코드

with_progress({
  p <- progressor(along = seq_len(nrow(tasks)))
  
  res_list <- furrr::future_map(
    seq_len(nrow(tasks)), # i,iter 별
    .options = furrr::furrr_options(seed = TRUE),  # L'Ecuyer-CMRG 스트림 자동 관리
    .f = function(k) { # 적용함수
      p()  # progress
      
      i    <- tasks$i[k]
      iter <- tasks$iter[k]
      
      m  <- param.grid$m[i]
      re <- param.grid$re[i]
      
      m.idx  <- which(m.list  == m)
      re.idx <- which(re.list == re)
      
      # D 벡터 구성 (원 코드 로직 유지)
      D <- rep(V.list, m / length(V.list))
      
      # --- 데이터 생성 ---
      data <- sim.data(m = m, re = re, D = D)
      
      # --- 모델 적합(한 작업 안에서 순차 실행: 데이터 재사용) ---
      #model fitting
      GLlasso <- func.GLlasso(n_iter = n_iter, burn_in = burn_in, data = data, hyper.s = 6, hyper.t = 6,ga = 0.01, gb = 0.01, gr = 0.01, hr = 0.01, thin = 10)
      FH <- func.FH(n_iter = n_iter, burn_in = burn_in, data= data)
      BEN <- func.ben(n_iter = n_iter, burn_in = burn_in, data = data, mc_iter = 30,  hyper.s=6, hyper.t=6, tune=0.005, init_rho = 0.5, thin = 10)
      DM <- func.DM(n_iter = n_iter, burn_in = burn_in, data = data)
      TPB <- func.tpb(n_iter = n_iter, burn_in = burn_in, p = 1, q = 1, data = data) #p and q?
      NG <- func.ng(n_iter = n_iter, burn_in = burn_in, data=data, a = 0.5)  # a?
      
      models <- list(GLlasso = GLlasso, FH = FH, BEN = BEN, DM = DM, TPB = TPB, NG = NG)
      
      true.theta <- data$theta
      
      # 5개 지표 x 3모델 매트릭스 반환
      metrics_mat <- matrix(NA_real_, nrow = 5L, ncol = length(models),
                            dimnames = list(c("AAD","ASD","ARB","ASRB","Coverage"),
                                            c("GLlasso","FH","BEN",'DM', 'TPB','NG')))
      j <- 1L
      for (nm in names(models)) {
        fit <- models[[nm]]
        metrics_mat[1, j] <- AAD(true.theta, fit$theta, m)
        metrics_mat[2, j] <- ASD(true.theta, fit$theta, m)
        metrics_mat[3, j] <- ARB(true.theta, fit$theta, m)
        metrics_mat[4, j] <- ASRB(true.theta, fit$theta, m)
        metrics_mat[5, j] <- Coverage(true.theta, fit$coverage)
        j <- j + 1L
      }
      
      # 인덱스와 결과를 함께 반환 (공유 메모리 안 쓰기)
      list(m.idx = m.idx, re.idx = re.idx, iter = iter, metrics = metrics_mat)
    }
  )
})

## 5) 결과 배열에 채우기 (단일 스레드)
for (r in res_list) {
  for (model.idx in 1:n.models) {
    for (metric.idx in 1:length(eval.list)) {
      result[r$m.idx, r$re.idx, metric.idx, r$iter, model.idx] <-
        r$metrics[metric.idx, model.idx]
    }
  }
}

####### Visualization ################################################################################################

















