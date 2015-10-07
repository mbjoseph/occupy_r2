data {
  int<lower=0> nsite;
  int<lower=0> nsurvey;
  int<lower=0,upper=1> y[nsite,nsurvey];
  int<lower=0> n_psicov;
  matrix[nsite, n_psicov] x_psi;
  int<lower=0> n_pcov;
  matrix[nsite, n_pcov] x_p;
}

parameters {
  real mu_psi;
  vector[n_psicov] beta_psi;
  real mu_p;
  vector[n_pcov] beta_p;
}

transformed parameters {
  vector[nsite] logit_psi;
  vector[nsite] psi;
  vector[nsite] logit_p;

  logit_psi <- mu_psi + x_psi * beta_psi;
  
  for (i in 1:nsite){
    psi[i] <- inv_logit(logit_psi[i]);
  }
  
  logit_p <- mu_p + x_p * beta_p;
}

model {
  // priors
  mu_psi ~ normal(0, 3);
  beta_psi ~ normal(0, 3);
  mu_p ~ normal(0, 3);
  beta_p ~ normal(0, 3);

  // likelihood
  {
  // local variables to avoid recomputing log(psi1) and log(1 - psi1)
  vector[nsite] log_psi;
  real log1m_psi;

  log_psi <- log(psi);
  
  for (i in 1:nsite) {
    log1m_psi <- log1m(psi[i]);
    if (sum(y[i]) > 0)
      increment_log_prob(log_psi[i] + bernoulli_logit_log(y[i], logit_p[i]));
    else
      increment_log_prob(
        log_sum_exp(log_psi[i] + bernoulli_logit_log(y[i], logit_p[i]),
                    log1m_psi));
  }
  }
}

generated quantities {
  real<lower=0, upper=1> rsq_psi;
  real<lower=0, upper=1> rsq_p;
  vector[nsite] XBpsi;
  vector[nsite] XBp;
  real<lower=0> var_XBpsi;
  real<lower=0> var_XBp;
  
  XBpsi <- x_psi * beta_psi;
  XBp <- x_p * beta_p;
  
  var_XBpsi <- sum((XBpsi - mean(XBpsi)) .* (XBpsi - mean(XBpsi))) / nsite;
  var_XBp <- sum((XBp - mean(XBp)) .* (XBp - mean(XBp))) / nsite;
  
  rsq_psi <- var_XBpsi / (var_XBpsi + pi()^2 / 3);
  rsq_p <- var_XBp / (var_XBp + pi()^2 / 3);
}