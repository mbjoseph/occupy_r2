# R^2 in logit psi
r2 <- function(X, beta){
  # X is a design matrix
  # beta is a matrix (n_iter X p) of covariate posteriors
  pop_var <- function(x){
    sum((x - mean(x))^2) / length(x)
  }
  nsite <- nrow(X)
  ndraws <- nrow(beta)
  XB <- array(dim=c(nsite, ndraws))
  for (i in 1:ndraws){
    XB[, i] <- X %*% beta[i, ]
  }
  varXB <- apply(XB, 2, pop_var)
  r_squared <- varXB / (varXB + pi^2 / 3)
  return(r_squared)
}