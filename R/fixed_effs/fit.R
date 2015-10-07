# Fitting the occupancy model with simple structure
library(rstan)
library(scales)
source('R/fixed_effs/sim_occ_d.R')
source('R/fixed_effs/r2.R')
nsite <- 100
nsurvey <- 10

occ_d <- sim_occ_d(nsite, nsurvey, npsicov = 1, npcov = 1)
occ_d$beta_psi
ymat <- matrix(occ_d$df$y, nrow=nsite, ncol=nsurvey, byrow=TRUE)
stan_d <- list(nsite=nsite, 
               n_psicov=ncol(occ_d$x_psi), x_psi=occ_d$x_psi, 
               n_pcov = ncol(occ_d$x_p), x_p = occ_d$x_p,
               nsurvey=nsurvey, y=ymat)
pars <- c('mu_psi', 'beta_psi', 'mu_p', 'beta_p', 'psi', 'rsq_psi', 'rsq_p')
m_init <- stan('R/fixed_effs/occ.stan', data=stan_d, chains=1, iter=10)
m_fit <- stan(fit=m_init, data=stan_d, pars=pars, chains=2, cores=2, iter=3000)
traceplot(m_fit)
post <- rstan::extract(m_fit)

par(mfrow=c(1, 3), bty='n')
rs_psi <- r2(X=stan_d$x_psi, beta = post$beta_psi)
br <- seq(0, 1, .02)
hist(rs_psi, breaks=br, freq=F,col='blac k', border=F,
     main=paste('beta_psi = ', round(occ_d$beta_psi, 2)))
lines(density(post$rsq_psi), col='red', lwd=2)

# R^2 in p
rs_p <- r2(X=stan_d$x_p, beta = post$beta_p)
hist(rs_p, breaks=br, freq=F, col='black', border=F,
     main=paste('beta_p = ', round(occ_d$beta_p, 2)))
lines(density(post$rsq_p), col='red', lwd=2)

plot(rs_psi, rs_p, xlim=c(0, 1), ylim=c(0, 1), col=alpha(1, .3))
abline(0, 1, lty=2)
