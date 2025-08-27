library(coda)
library(Cairo)

draw_plots_GLlasso <- function(result){
pdf("traceplots.pdf", width=5, height=5)   # 변수마다 PNG 생성
  for (i in 1:dim(result$beta)[2]){
    plot(result$beta[,i], type="l", ylab = paste0("beta", i)) 
  }

  for (i in 1:dim(result$v)[2]){
    plot(result$v[,i], type="l", ylab = paste0("v", i)) 
  }

  # for (i in 1:dim(result$lambda)[2]){
  #   plot(result$lambda[,i], type="l", ylab = paste0("lambda", i)) 
  # }

  plot(result$sigma, type="l", ylab = 'sigma') 
  # only for GLlasso
  plot(result$a, type="l", ylab = 'a') 
  plot(result$b, type="l", ylab = 'b') 
  plot(result$r, type="l", ylab = 'r') 
  
  dev.off()
}
result <- GLlasso_results
draw_plots_GLlasso(result)

draw_plots_BEN <- function(result){
  pdf("traceplots.pdf", width=5, height=5)   # 변수마다 PNG 생성
  for (i in 1:dim(result$beta)[2]){
    plot(result$beta[,i], type="l", ylab = paste0("beta", i)) 
  }
  
  for (i in 1:dim(result$v)[2]){
    plot(result$v[,i], type="l", ylab = paste0("v", i)) 
  }
  
  # for (i in 1:dim(result$lambda)[2]){
  #   plot(result$lambda[,i], type="l", ylab = paste0("lambda", i)) 
  # }
  
  plot(result$sigma, type="l", ylab = 'sigma') 
  dev.off()
}
result <- BEN_results
draw_plots_BEN(result)

draw_plots_FH <- function(result){
  pdf("traceplots.pdf", width=5, height=5)   # 변수마다 PNG 생성
  for (i in 1:dim(result$Theta.chain)[1]){
    plot(result$Theta.chain[i,], type="l", ylab = paste0("Theta", i)) 
  }
    plot(result$Sigma2.chain, type = 'l' , ylab = 'sigma')
    
  dev.off()
}
result <- FH_results
draw_plots_FH(result)

draw_plots_DM <- function(result){
  pdf("traceplots.pdf", width=5, height=5)   # 변수마다 PNG 생성
  
  for (i in 1:dim(result$Beta.chain)[1]){
    plot(result$Beta.chain[i,], type="l", ylab = paste0("Beta", i)) 
  }
  
  for (i in 1:dim(result$V.chain)[1]){
    plot(result$V.chain[i,], type="l", ylab = paste0("V", i)) 
  }
  
  for (i in 1:dim(result$Delta.chain)[1]){
    plot(result$Delta.chain[i,], type="l", ylab = paste0("Delta", i)) 
  }
  
  plot(result$Sigma2.chain, type = 'l' , ylab = 'sigma')
  
  dev.off()
}
result <- DM_results
draw_plots_DM(result)

draw_plots_HS <- function(result){
  pdf("traceplots.pdf", width=5, height=5)   # 변수마다 PNG 생성
  
  for (i in 1:dim(result$Beta.chain)[1]){
    plot(result$Beta.chain[i,], type="l", ylab = paste0("Beta", i)) 
  }
  
  for (i in 1:dim(result$U.chain)[1]){
    plot(result$U.chain[i,], type="l", ylab = paste0("U", i))
  }

  for (i in 1:dim(result$Lambda2.chain)[1]){
    plot(result$Lambda2.chain[i,], type="l", ylab = paste0("Lambda2", i))
  }

  plot(result$Tau2.chain, type = 'l' , ylab = 'tau2')

  dev.off()
}
result <- HS_results
result <- SB_results
result <- NEG_results
result <- LA_results
result <- NG_results
draw_plots_HS(result)
