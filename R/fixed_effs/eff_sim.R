# Exploring the relationship between cov. effect size and R^2
# Fitting the occupancy model with simple structure
library(rstan)
library(scales)
library(modeest)
source('R/fixed_effs/sim_occ_d.R')
source('R/HDI.R')
nsite <- 50
nsurvey <- 2
iter <- 1200
n_sim <- 10

# initialize model
occ_d <- sim_occ_d(nsite, nsurvey, npsicov = 1, npcov = 1)
ymat <- matrix(occ_d$df$y, nrow=nsite, ncol=nsurvey, byrow=TRUE)
stan_d <- list(nsite=nsite, 
               n_psicov=ncol(occ_d$x_psi), x_psi=occ_d$x_psi, 
               n_pcov = ncol(occ_d$x_p), x_p = occ_d$x_p,
               nsurvey=nsurvey, y=ymat)
pars <- c('mu_psi', 'beta_psi', 'mu_p', 'beta_p', 'psi', 'rsq_psi', 'rsq_p')
m_init <- stan('R/fixed_effs/occ.stan', data=stan_d, chains=1, iter=10)

# iteratively simulate datasets and store covariate effects with R2
psi_beta <- array(dim=c(n_sim))
p_beta <- array(dim=c(n_sim))
psi_r2 <- array(dim=c(n_sim))
p_r2 <- array(dim=c(n_sim))
psi_lo <- array(dim=c(n_sim))
psi_hi <- array(dim=c(n_sim))
p_lo <- array(dim=c(n_sim))
p_hi <- array(dim=c(n_sim))

for (i in 1:n_sim){
  occ_d <- sim_occ_d(nsite, nsurvey, npsicov = 1, npcov = 1)
  ymat <- matrix(occ_d$df$y, nrow=nsite, ncol=nsurvey, byrow=TRUE)
  stan_d <- list(nsite=nsite, 
                 n_psicov=ncol(occ_d$x_psi), x_psi=occ_d$x_psi, 
                 n_pcov = ncol(occ_d$x_p), x_p = occ_d$x_p,
                 nsurvey=nsurvey, y=ymat)
  m_fit <- stan('R/fixed_effs/occ.stan', 
                data=stan_d, pars=pars, refresh=0, save_dso=FALSE,
                chains=2, cores=1, iter=iter, open_progress=FALSE, 
                verbose=FALSE)
  if (summary(m_fit)$summary['lp__', 'Rhat'] > 1.01) next
  post <- extract(m_fit)
  psi_beta[i] <- occ_d$beta_psi
  p_beta[i] <- occ_d$beta_p
  
  psi_r2[i] <- mlv(post$rsq_psi, method='shorth')$M
  hdi_psi <- HDI(post$rsq_psi)
  psi_lo[i] <- hdi_psi[1]
  psi_hi[i] <- hdi_psi[2]

  p_r2[i] <- mlv(post$rsq_p, method='shorth')$M
  hdi_p <- HDI(post$rsq_psi)
  p_lo[i] <- hdi_p[1]
  p_hi[i] <- hdi_p[2]

  closeAllConnections()
}

d <- data.frame(beta = c(psi_beta, p_beta), 
                r2 = c(psi_r2, p_r2), 
                parameter = factor(rep(c('psi', 'p'), each=n_sim), 
                                   labels = c('psi', 'p')), 
                lo = c(psi_lo, p_lo), 
                hi = c(psi_hi, p_hi))
d <- d[complete.cases(d), ]
# saveRDS(d, 'R/sim_d.rds')
# d <- readRDS('R/fixed_effs/sim_d.rds')
ggplot(d, aes(x=beta)) + 
  #geom_segment(aes(x=beta, xend=beta, y=lo, yend=hi), alpha=.3) +
  geom_point(aes(y=r2), alpha=.3) + 
  facet_grid(.~parameter, labeller= label_parsed) + 
  xlab(expression(paste('Coefficient: ', beta))) + 
  ylab(expression(paste('Posterior mode: ', R^2)))

plot(psi_beta, p_r2, cex=.1)
plot(p_beta, psi_r2, cex=.1)
