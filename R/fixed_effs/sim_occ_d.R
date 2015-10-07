# Function to simulate occupancy data
# - normally distributed site-varying occupancy probabilities
# - normally distributed site-varying detection probabilities
# - normally distributed survey-varying detection probabilities
sim_occ_d <- function(nsite, nsurvey, 
                      npsicov=1, npcov=1){
  require(testit)
  x_psi <- matrix(rnorm(nsite * npsicov), ncol=npsicov)
  beta_psi <- runif(npsicov, -8, 8)
  # occupancy parameters
  mu_psi <- 0
  lpsi <- mu_psi + x_psi %*% beta_psi
  psi <- plogis(lpsi)
  z <- rbinom(nsite, 1, psi)
  
  # detection parameters
  site <- rep(1:nsite, nsurvey)
  survey <- rep(1:nsurvey, nsite)
  lmup <- 0
  x_p <- matrix(rnorm(nsite * npcov), ncol=npcov)
  beta_p <- runif(npcov, -8, 8)
  lp <- lmup + x_p %*% beta_p
  p <- plogis(lp)
  
  y <- rbinom(nsite*nsurvey, 1, z[site] * p[site])
  
  assert("No detections without occurence", !any(y==1 & rep(z, nsurvey) == 0))
  res <- data.frame(y, site, survey)
  res <- res[order(res$site, res$survey), ]
  return(list(df=res, x_psi=x_psi, x_p = x_p, lpsi=lpsi, z=z, 
              beta_psi=beta_psi, beta_p=beta_p))
}
