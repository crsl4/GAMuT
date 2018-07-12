library(CompQuadForm)

proj_GAMuT_pheno <- function(X){
  out = vector("list", 2)
  names(out) = c("Kc", "ev_Kc")
  out$Kc = X %*% solve(t(X) %*% X) %*% t(X) # Projection matrix
  out$ev_Kc = rep(1, ncol(X)) # Find eigenvalues
  return(out)	
}

linear_GAMuT_geno <- function(X){
  out = vector("list", 2)
  names(out) = c("Lc", "ev_Lc")
  
  Lc = t(X) %*% X # transposed kernel to find eigenvalues: runs faster
  ev_Lc = eigen(Lc, symmetric=TRUE, only.values=TRUE)$values  
  
  out$Lc = X %*% t(X) # similiarity kernel for test
  out$ev_Lc = ev_Lc[ev_Lc > 1e-08]
  return(out)	 
}

TestGAMuT <- function(Yc, lambda_Y, Xc, lambda_X) {
  
  ## Test statistic
  m = nrow(Yc) # Number of subjects in study
  GAMuT = (1/m) * sum(sum(t(Yc) %*% Xc))  
  
  ## Derive p-value of GAMuT statistic
  ## Form vector of eigenvalue products
  Z = (as.matrix(lambda_Y)) %*% t(as.matrix(lambda_X))
  Zsort <- sort(Z, decreasing=TRUE)
  scoredavies = GAMuT*m^2
  results_score <- davies(scoredavies, Zsort)
  davies_pvalue <- (results_score$Qq)
  return(davies_pvalue)
}
