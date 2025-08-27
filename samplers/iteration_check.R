# for parallel computing
library(parallel)
source('samplers/lasso_sampler.R')


mydata <- read.table(file="analysis/data1.txt", header=TRUE)

m <- nrow(mydata)
X <- cbind(rep(1, m), mydata$X1, mydata$X2, mydata$X3)
Y <- mydata$Y
D <- mydata$d

n_iter <- nburn+nsim
data <- list(Y = Y, X = X, m = nrow(X), p = ncol(X), theta = NA , beta = NA , u = NA , D = D)

GLlasso_results <- func.GLlasso(n_iter = n_iter, burn_in = nburn, data = data, hyper.s = 6, hyper.t = 6,ga = 0.01, gb = 0.01, gr = 0.01, hr = 0.01, thin = 10)


task_list <- list(
  list(fun = function() func.GLlasso(n_iter = n_iter, burn_in = nburn, data = data, 
                                     hyper.s = 6, hyper.t = 6,ga = 0.01, gb = 0.01, gr = 0.01, hr = 0.01, thin = 10))
)
task_list <- unlist(lapply(task_list, function(x) rep(list(x),4)), recursive = FALSE)

cl <- makeCluster(20)
clusterExport(cl, varlist = c('data','func.GLlasso', 'nburn','n_iter'))

results <- parLapply(cl, task_list, function(task) {
  tryCatch({
    task$fun()
  }, error = function(e) {
    message("Error in task: ", conditionMessage(e))
    return(NA)  # 또는 NULL
  })
})
stopCluster(cl)


results