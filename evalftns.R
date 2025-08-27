#true, pred : mx1 vectpr

AAD <- function(true, pred, m){
  return(sum(abs(true-pred))/m)
}


ASD <- function(true, pred, m){
  return(sum((true-pred)**2)/m)
}


ARB <- function(true, pred, m){
  return(sum(abs((pred-true)/true))/m)
}


ASRB <- function(true, pred, m){
  return(sum(((true-pred)/true)**2)/m)
}

Coverage <- function(true, interval){
  return(mean((true > interval[,1])*(true < interval[,2])))
}

