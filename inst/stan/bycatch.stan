data {
  int<lower=0> n_row;
  vector[n_row] effort; // covariates
  vector[n_row] new_effort; // covariate for unobserved sets
  int yint[n_row]; // vector[n_pos] y;
  vector[n_row] yreal;
  int time[n_row]; // time variable
  int<lower=0> n_year; // number of unique years
  int<lower=0> K;
  matrix[n_row, K] x; // covariates
  int family; // 1 = poisson, 2 = negbin, 3 = hurdle-poisson, 4 = hurdle-negbin, 5 = lognormal, 6 = gamma, 7 = hurdle-lognormal, 8 = hurdle-gamma, 9 = normal, 10 = hurdle-normal
  int time_varying; // whether to treat model as dlm
}
transformed data {
  int est_phi;
  int est_theta;
  int est_sigma;
  int est_cv;
  int is_discrete;
  // Count zeros and non-zeros for hurdle models
  int n_zero = 0;
  int n_nonzero = 0;
  int zero_idx[n_row];
  int nonzero_idx[n_row];

  // initialize
  est_phi = 0;
  est_theta = 0;
  est_sigma = 0;
  est_cv = 0;
  is_discrete = 0;

  if(family == 2) est_phi = 1;
  if(family == 3) est_theta = 1;
  if(family == 4) {
    est_phi = 1;
    est_theta = 1;
  }
  if(family < 5) is_discrete = 1;

  // estimate lognormal variance
  if(family == 5) est_sigma = 1;
  if(family == 7) {
    est_sigma = 1;
    est_theta = 1;
  }

  // estimate gamma cv
  if(family == 6) est_cv = 1;
  if(family == 8) {
    est_cv = 1;
    est_theta = 1;
  }

  // estimate normal variance
  if(family == 9) est_sigma = 1;
  if(family == 10) {
    est_sigma = 1;
    est_theta = 1;
  }

  // Pre-compute indices for hurdle models (vectorization)
  if(est_theta == 1) {
    for(i in 1:n_row) {
      if(yint[i] == 0) {
        n_zero += 1;
        zero_idx[n_zero] = i;
      } else {
        n_nonzero += 1;
        nonzero_idx[n_nonzero] = i;
      }
    }
  }
}
parameters {
  vector[K] beta;
  vector[time_varying*(n_year-1)] est_time_dev;
  real<lower=0> sigma_rw[time_varying];
  real<lower=0> sigma_logn[est_sigma];
  real<lower=0> cv_gamma[est_cv];
  real<lower=0> nb2_phi[est_phi];
  real<lower=0,upper=1> theta[est_theta];
}
transformed parameters {
  vector[n_row] log_lambda;
  vector[n_row] lambda;
  vector[n_row] pred;
  real<lower=0> gammaA[est_cv];
  vector[time_varying*n_year] time_dev;
  pred = x * beta;

  if(time_varying == 1) {
    time_dev[1] = 0;
    for(i in 2:n_year) {
      time_dev[i] = est_time_dev[i-1];
    }
  }

  for(i in 1:n_row) {
    if(time_varying == 1) {
      pred[i] = pred[i] + time_varying*time_dev[time[i]];
    }
    log_lambda[i] = pred[i] + log(effort[i]); // exp(pred) = theta
    lambda[i] = exp(log_lambda[i]);
  }
  // gamma model
  if(est_cv == 1) gammaA[1] = inv(pow(cv_gamma[1], 2.0));
}
model {
  beta ~ student_t(3, 0, 2);

  if(time_varying == 1) {
    sigma_rw ~ student_t(3, 0, 1);
    est_time_dev[1] ~ student_t(3, 0, 2);
    for(i in 2:(n_year-1)) {
      // model evolution of temporal deviations as random walk
      est_time_dev[i] ~ normal(est_time_dev[i-1], sigma_rw[1]);
    }
  }

  if(est_theta == 1) {
    theta ~ beta(1, 1);
  }

  // Family 1: Poisson
  if(family == 1) {
    yint ~ poisson_log(log_lambda);
  }
  // Family 2: Negative Binomial
  else if(family == 2) {
    nb2_phi ~ student_t(3, 0, 2);
    yint ~ neg_binomial_2_log(log_lambda, nb2_phi[1]);
  }
  // Family 3: Poisson-Hurdle (vectorized)
  else if(family == 3) {
    // Use sufficient statistics for efficiency
    target += n_zero * log(theta[1]);
    target += n_nonzero * log1m(theta[1]);

    // Handle positive counts - vectorize where possible
    for(i in 1:n_nonzero) {
      int idx = nonzero_idx[i];
      target += poisson_log_lpmf(yint[idx] | log_lambda[idx]) -
                poisson_lccdf(0 | lambda[idx]);
    }
  }
  // Family 4: Negative Binomial-Hurdle (vectorized)
  else if(family == 4) {
    nb2_phi ~ student_t(3, 0, 2);

    // Use sufficient statistics for efficiency
    target += n_zero * log(theta[1]);
    target += n_nonzero * log1m(theta[1]);

    // Handle positive counts
    for(i in 1:n_nonzero) {
      int idx = nonzero_idx[i];
      target += neg_binomial_2_log_lpmf(yint[idx] | log_lambda[idx], nb2_phi[1]) -
                neg_binomial_2_lccdf(0 | lambda[idx], nb2_phi[1]);
    }
  }
  // Family 5: Lognormal
  else if(family == 5) {
    sigma_logn ~ student_t(3, 0, 2);
    yreal ~ lognormal(log_lambda, sigma_logn[1]);
  }
  // Family 6: Gamma
  else if(family == 6) {
    cv_gamma[1] ~ student_t(3, 0, 2);
    yreal ~ gamma(gammaA[1], gammaA[1] ./ lambda);
  }
  // Family 7: Lognormal-Hurdle (vectorized)
  else if(family == 7) {
    sigma_logn ~ student_t(3, 0, 2);

    // Use sufficient statistics for efficiency
    target += n_zero * log(theta[1]);
    target += n_nonzero * log1m(theta[1]);

    // Handle positive values
    for(i in 1:n_nonzero) {
      int idx = nonzero_idx[i];
      target += lognormal_lpdf(yreal[idx] | log_lambda[idx], sigma_logn[1]);
    }
  }
  // Family 8: Gamma-Hurdle (vectorized)
  else if(family == 8) {
    cv_gamma[1] ~ student_t(3, 0, 2);

    // Use sufficient statistics for efficiency
    target += n_zero * log(theta[1]);
    target += n_nonzero * log1m(theta[1]);

    // Handle positive values
    for(i in 1:n_nonzero) {
      int idx = nonzero_idx[i];
      target += gamma_lpdf(yreal[idx] | gammaA[1], gammaA[1] / lambda[idx]);
    }
  }
  // Family 9: Normal
  else if(family == 9) {
    sigma_logn ~ student_t(3, 0, 2);
    yreal ~ normal(lambda, sigma_logn[1]);
  }
  // Family 10: Normal-Hurdle (vectorized)
  else if(family == 10) {
    sigma_logn ~ student_t(3, 0, 2);

    // Use sufficient statistics for efficiency
    target += n_zero * log(theta[1]);
    target += n_nonzero * log1m(theta[1]);

    // Handle positive values
    for(i in 1:n_nonzero) {
      int idx = nonzero_idx[i];
      target += normal_lpdf(yreal[idx] | lambda[idx], sigma_logn[1]);
    }
  }
}
generated quantities {
  vector[n_row] log_lik;
  int<lower = 0> y_new[n_row*is_discrete];
  vector[n_row*(1-is_discrete)] y_new_real;

  // Pre-compute log_new_lambda
  vector[n_row] log_new_lambda;
  for(n in 1:n_row) {
    log_new_lambda[n] = pred[n] + log(new_effort[n]);
  }

  // Family 1: Poisson
  if(family == 1) {
    for(n in 1:n_row) {
      log_lik[n] = poisson_log_lpmf(yint[n] | log_lambda[n]);

      // sample posterior predictive distribution
      y_new[n] = 0;
      if(new_effort[n] > 0) {
        y_new[n] = poisson_log_rng(log_new_lambda[n]);
      }
    }
  }
  // Family 2: Negative Binomial
  else if(family == 2) {
    for(n in 1:n_row) {
      log_lik[n] = neg_binomial_2_log_lpmf(yint[n] | log_lambda[n], nb2_phi[1]);

      // sample posterior predictive distribution
      y_new[n] = 0;
      if(new_effort[n] > 0) {
        y_new[n] = neg_binomial_2_log_rng(log_new_lambda[n], nb2_phi[1]);
      }
    }
  }
  // Family 3: Poisson-Hurdle
  else if(family == 3) {
    for(n in 1:n_row) {
      if(yint[n] == 0) {
        log_lik[n] = log(theta[1]);
      } else {
        // (1 - theta) * Pr(Pois(y|lambda)) / (1 - PoissonCDF(0|lambda))
        log_lik[n] = log1m(theta[1]) + poisson_log_lpmf(yint[n] | log_lambda[n]) -
                     poisson_lccdf(0 | lambda[n]);
      }

      // sample posterior predictive distribution
      y_new[n] = 0;
      if(new_effort[n] > 0) {
        y_new[n] = (1 - bernoulli_rng(theta[1])) * poisson_log_rng(log_new_lambda[n]);
      }
    }
  }
  // Family 4: Negative Binomial-Hurdle
  else if(family == 4) {
    for(n in 1:n_row) {
      if(yint[n] == 0) {
        log_lik[n] = log(theta[1]);
      } else {
        // (1 - theta) * Pr(NB2(y|lambda)) / (1 - NB2CDF(0|lambda))
        log_lik[n] = log1m(theta[1]) + neg_binomial_2_log_lpmf(yint[n] | log_lambda[n], nb2_phi[1]) -
                     neg_binomial_2_lccdf(0 | lambda[n], nb2_phi[1]);
      }

      // sample posterior predictive distribution
      y_new[n] = 0;
      if(new_effort[n] > 0) {
        y_new[n] = (1 - bernoulli_rng(theta[1])) * neg_binomial_2_log_rng(log_new_lambda[n], nb2_phi[1]);
      }
    }
  }
  // Family 5: Lognormal
  else if(family == 5) {
    for(n in 1:n_row) {
      log_lik[n] = lognormal_lpdf(yreal[n] | log_lambda[n], sigma_logn[1]);

      // sample posterior predictive distribution
      y_new_real[n] = 0;
      if(new_effort[n] > 0) {
        y_new_real[n] = lognormal_rng(log_new_lambda[n], sigma_logn[1]);
      }
    }
  }
  // Family 6: Gamma
  else if(family == 6) {
    for(n in 1:n_row) {
      log_lik[n] = gamma_lpdf(yreal[n] | gammaA[1], gammaA[1] / lambda[n]);

      // sample posterior predictive distribution
      y_new_real[n] = 0;
      if(new_effort[n] > 0) {
        y_new_real[n] = gamma_rng(gammaA[1], gammaA[1] / exp(log_new_lambda[n]));
      }
    }
  }
  // Family 7: Lognormal-Hurdle
  else if(family == 7) {
    for(n in 1:n_row) {
      if(yint[n] == 0) {
        log_lik[n] = log(theta[1]);
      } else {
        log_lik[n] = log1m(theta[1]) + lognormal_lpdf(yreal[n] | log_lambda[n], sigma_logn[1]);
      }

      // sample posterior predictive distribution
      y_new_real[n] = 0;
      if(new_effort[n] > 0) {
        y_new_real[n] = (1 - bernoulli_rng(theta[1])) * lognormal_rng(log_new_lambda[n], sigma_logn[1]);
      }
    }
  }
  // Family 8: Gamma-Hurdle
  else if(family == 8) {
    for(n in 1:n_row) {
      if(yint[n] == 0) {
        log_lik[n] = log(theta[1]);
      } else {
        log_lik[n] = log1m(theta[1]) + gamma_lpdf(yreal[n] | gammaA[1], gammaA[1] / lambda[n]);
      }

      // sample posterior predictive distribution
      y_new_real[n] = 0;
      if(new_effort[n] > 0) {
        // Fixed: removed double semicolon (improvement 1.2)
        y_new_real[n] = (1 - bernoulli_rng(theta[1])) * gamma_rng(gammaA[1], gammaA[1] / exp(log_new_lambda[n]));
      }
    }
  }
  // Family 9: Normal
  else if(family == 9) {
    for(n in 1:n_row) {
      log_lik[n] = normal_lpdf(yreal[n] | lambda[n], sigma_logn[1]);

      // sample posterior predictive distribution
      y_new_real[n] = 0;
      if(new_effort[n] > 0) {
        y_new_real[n] = normal_rng(exp(log_new_lambda[n]), sigma_logn[1]);
      }
    }
  }
  // Family 10: Normal-Hurdle
  else if(family == 10) {
    for(n in 1:n_row) {
      if(yint[n] == 0) {
        log_lik[n] = log(theta[1]);
      } else {
        log_lik[n] = log1m(theta[1]) + normal_lpdf(yreal[n] | lambda[n], sigma_logn[1]);
      }

      // sample posterior predictive distribution
      y_new_real[n] = 0;
      if(new_effort[n] > 0) {
        y_new_real[n] = (1 - bernoulli_rng(theta[1])) * normal_rng(exp(log_new_lambda[n]), sigma_logn[1]);
      }
    }
  }
}
