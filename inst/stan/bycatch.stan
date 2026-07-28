data {
  int<lower=0> n_row;
  vector[n_row] effort; // For backwards compatibility (single stream)
  vector[n_row] new_effort; // covariate for unobserved sets
  int yint[n_row]; // For backwards compatibility (single stream)
  vector[n_row] yreal; // For backwards compatibility (single stream)
  int time[n_row]; // time variable
  int<lower=0> n_year; // number of unique years
  int<lower=0> K;
  matrix[n_year, K] x; // covariates (year-level)
  int family; // 1 = poisson, 2 = negbin, 3 = hurdle-poisson, 4 = hurdle-negbin, 5 = lognormal, 6 = gamma, 7 = hurdle-lognormal, 8 = hurdle-gamma, 9 = normal, 10 = hurdle-normal
  int time_varying; // whether to treat model as dlm

  // Multi-stream data
  int<lower=0,upper=1> multi_stream;  // multi-stream mode 0/1
  int<lower=0> n_obs;     // Number of obs-only observations
  int<lower=0> n_em;      // Number of EM-only observations
  int<lower=0> n_both;    // Number of both observations
  int<lower=0> n_selfreport; // Number of self-report observations

  // Separate data streams (only used if multi_stream == 1)
  array[n_obs] int yint_obs;
  array[n_em] int yint_em;
  array[n_both] int yint_both;
  array[n_selfreport] int yint_selfreport;
  vector[n_obs] effort_obs;
  vector[n_em] effort_em;
  vector[n_both] effort_both;
  vector[n_selfreport] effort_selfreport;
  vector[n_obs] yreal_obs;
  vector[n_em] yreal_em;
  vector[n_both] yreal_both;
  vector[n_selfreport] yreal_selfreport;

  // Indices mapping stream observations to time periods
  array[n_obs] int time_idx_obs;
  array[n_em] int time_idx_em;
  array[n_both] int time_idx_both;
  array[n_selfreport] int time_idx_selfreport;

  // Effort for unobserved expansion (by time period)
  vector[n_year] new_effort_by_year;

  // ---- NEW: self-reporting probability sub-model ----
  // Ground truth comes from OBS/EM/BOTH-covered vessels only, where we know
  // both monitoring status and self-report status. Pooled (constant) across years.
  int<lower=0,upper=1> estimate_self_report;   // turn this sub-model on/off
  array[n_year] int<lower=0> n_vessels_observed;      // denominator: all observed (obs/em/both) vessels
  array[n_year] int<lower=0> n_vessels_self_report;   // numerator: subset that also self-report
}
transformed data {
  int est_phi;
  int est_theta;
  int est_sigma;
  int est_cv;
  int is_discrete;

  // Total monitored effort (obs+em+both) by year, and total self-report stream
  // effort by year. Used in generated quantities to figure out how much of the
  // self-report stream's effort is genuinely NEW coverage vs. overlap with
  // vessels already captured by obs/em/both.
  vector[n_year] effort_monitored_by_year;
  vector[n_year] effort_selfreport_by_year;

  // initial
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

  // ---- NEW: sum monitored (obs+em+both) and self-report effort by year ----
  for(t in 1:n_year) {
    effort_monitored_by_year[t] = 0;
    effort_selfreport_by_year[t] = 0;
  }
  for(i in 1:n_obs)  effort_monitored_by_year[time_idx_obs[i]]  += effort_obs[i];
  for(i in 1:n_em)   effort_monitored_by_year[time_idx_em[i]]   += effort_em[i];
  for(i in 1:n_both) effort_monitored_by_year[time_idx_both[i]] += effort_both[i];
  for(i in 1:n_selfreport) effort_selfreport_by_year[time_idx_selfreport[i]] += effort_selfreport[i];
}
parameters {
  vector[K] beta;
  vector[time_varying*(n_year-1)] est_time_dev;
  real<lower=0> sigma_rw[time_varying];
  real<lower=0> sigma_logn[est_sigma];
  real<lower=0> cv_gamma[est_cv];
  real<lower=0> nb2_phi[est_phi];
  real<lower=0,upper=1> theta[est_theta];

  // ---- NEW: self-reporting sub-model parameters ----
  // Single pooled logit probability - self-reporting compliance assumed constant over time.
  vector[estimate_self_report] logit_p_report;
}
transformed parameters {
  vector[n_year] log_lambda_base;  // Base lambda for each year
  vector[n_year] lambda_base;
  vector[n_year] pred;
  real<lower=0> gammaA[est_cv];
  vector[time_varying*n_year] time_dev;

  // Single stream transformed parameters (backwards compatibility)
  vector[n_row] log_lambda;
  vector[n_row] lambda;

  // ---- NEW: self-reporting probability by year ----
  vector[n_year] p_report; // defaults to 0 (unused) when estimate_self_report == 0

  // base prediction for each year
  pred = x * beta;

  if(time_varying == 1) {
    time_dev[1] = 0;
    for(i in 2:n_year) {
      time_dev[i] = est_time_dev[i-1];
    }
    for(t in 1:n_year) {
      pred[t] = pred[t] + time_dev[t];
    }
  }

  // lambda for each year (without effort multiplier)
  for(t in 1:n_year) {
    log_lambda_base[t] = pred[t];
    lambda_base[t] = exp(log_lambda_base[t]);
  }

  // Single stream lambda (backwards compatibility)
  if(multi_stream == 0) {
    for(i in 1:n_row) {
      log_lambda[i] = pred[time[i]] + log(effort[i]);
      lambda[i] = exp(log_lambda[i]);
    }
  }

  // gamma model
  if(est_cv == 1) gammaA[1] = inv(pow(cv_gamma[1], 2.0));

  // ---- NEW: build p_report[t] (pooled/constant across years) ----
  for(t in 1:n_year) p_report[t] = 0;
  if(estimate_self_report == 1) {
    real p_pooled = inv_logit(logit_p_report[1]);
    for(t in 1:n_year) {
      p_report[t] = p_pooled;
    }
  }
}
model {
  beta ~ student_t(3, 0, 2);

  if(time_varying == 1) {
    sigma_rw ~ student_t(3, 0, 1);
    est_time_dev[1] ~ student_t(3, 0, 2);
    for(i in 2:(n_year-1)) {
      est_time_dev[i] ~ normal(est_time_dev[i-1], sigma_rw[1]);
    }
  }

  if(est_theta == 1) {
    theta ~ beta(1, 1);
  }

  // Single stream
  if(multi_stream == 0) {

    if(family == 1) {
      yint ~ poisson_log(log_lambda);
    }
    else if(family == 2) {
      nb2_phi ~ student_t(3, 0, 2);
      yint ~ neg_binomial_2_log(log_lambda, nb2_phi[1]);
    }
    else if(family == 3) {
      for(i in 1:n_row) {
        if (yint[i] == 0)
          1 ~ bernoulli(theta);
        else {
          0 ~ bernoulli(theta);
          yint[i] ~ poisson(lambda[i]) T[1, ];
        }
      }
    }
    else if(family == 4) {
      nb2_phi ~ student_t(3, 0, 2);
      for(i in 1:n_row) {
        if (yint[i] == 0)
          1 ~ bernoulli(theta);
        else {
          0 ~ bernoulli(theta);
          yint[i] ~ neg_binomial_2(lambda[i], nb2_phi[1]) T[1, ];
        }
      }
    }
    else if(family == 5) {
      sigma_logn ~ student_t(3, 0, 2);
      yreal ~ lognormal(log_lambda, sigma_logn[1]);
    }
    else if(family == 6) {
      cv_gamma[1] ~ student_t(3, 0, 2);
      yreal ~ gamma(gammaA[1], gammaA[1] ./ lambda);
    }
    else if(family == 7) {
      sigma_logn ~ student_t(3, 0, 2);
      for(i in 1:n_row) {
        if (yint[i] == 0)
          1 ~ bernoulli(theta);
        else {
          0 ~ bernoulli(theta);
          yreal[i] ~ lognormal(log_lambda[i], sigma_logn[1]);
        }
      }
    }
    else if(family == 8) {
      cv_gamma[1] ~ student_t(3, 0, 2);
      for(i in 1:n_row) {
        if (yint[i] == 0)
          1 ~ bernoulli(theta);
        else {
          0 ~ bernoulli(theta);
          yreal[i] ~ gamma(gammaA[1], gammaA[1] / lambda[i]);
        }
      }
    }
    else if(family == 9) {
      sigma_logn ~ student_t(3, 0, 2);
      yreal ~ normal(lambda, sigma_logn[1]);
    }
    else if(family == 10) {
      sigma_logn ~ student_t(3, 0, 2);
      for(i in 1:n_row) {
        if (yint[i] == 0)
          1 ~ bernoulli(theta);
        else {
          0 ~ bernoulli(theta);
          yreal[i] ~ normal(lambda[i], sigma_logn[1]);
        }
      }
    }
  }

  // Multi-stream (separate likelihoods, detection = 1)
  else {

    // Family 1: Poisson
    if(family == 1) {
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            yint_obs[i] ~ poisson(lambda_obs_i);
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            yint_em[i] ~ poisson(lambda_em_i);
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            yint_both[i] ~ poisson(lambda_both_i);
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            yint_selfreport[i] ~ poisson(lambda_selfreport_i);
          }
        }
      }
    }

    // Family 2: Negative Binomial
    else if(family == 2) {
      nb2_phi ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            yint_obs[i] ~ neg_binomial_2(lambda_obs_i, nb2_phi[1]);
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            yint_em[i] ~ neg_binomial_2(lambda_em_i, nb2_phi[1]);
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            yint_both[i] ~ neg_binomial_2(lambda_both_i, nb2_phi[1]);
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            yint_selfreport[i] ~ neg_binomial_2(lambda_selfreport_i, nb2_phi[1]);
          }
        }
      }
    }

    // Family 3: Poisson-Hurdle
    else if(family == 3) {
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            if (yint_obs[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_obs[i] ~ poisson(lambda_obs_i) T[1, ];
            }
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            if (yint_em[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_em[i] ~ poisson(lambda_em_i) T[1, ];
            }
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            if (yint_both[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_both[i] ~ poisson(lambda_both_i) T[1, ];
            }
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            if (yint_selfreport[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_selfreport[i] ~ poisson(lambda_selfreport_i) T[1, ];
            }
          }
        }
      }
    }

    // Family 4: Negative Binomial-Hurdle
    else if(family == 4) {
      nb2_phi ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            if (yint_obs[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_obs[i] ~ neg_binomial_2(lambda_obs_i, nb2_phi[1]) T[1, ];
            }
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            if (yint_em[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_em[i] ~ neg_binomial_2(lambda_em_i, nb2_phi[1]) T[1, ];
            }
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            if (yint_both[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_both[i] ~ neg_binomial_2(lambda_both_i, nb2_phi[1]) T[1, ];
            }
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            if (yint_selfreport[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yint_selfreport[i] ~ neg_binomial_2(lambda_selfreport_i, nb2_phi[1]) T[1, ];
            }
          }
        }
      }
    }

    // Family 5: Lognormal
    else if(family == 5) {
      sigma_logn ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0 && yreal_obs[i] > 0) {
            real log_lambda_obs_i = log_lambda_base[time_idx_obs[i]] + log(effort_obs[i]);
            yreal_obs[i] ~ lognormal(log_lambda_obs_i, sigma_logn[1]);
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0 && yreal_em[i] > 0) {
            real log_lambda_em_i = log_lambda_base[time_idx_em[i]] + log(effort_em[i]);
            yreal_em[i] ~ lognormal(log_lambda_em_i, sigma_logn[1]);
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0 && yreal_both[i] > 0) {
            real log_lambda_both_i = log_lambda_base[time_idx_both[i]] + log(effort_both[i]);
            yreal_both[i] ~ lognormal(log_lambda_both_i, sigma_logn[1]);
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0 && yreal_selfreport[i] > 0) {
            real log_lambda_selfreport_i = log_lambda_base[time_idx_selfreport[i]] + log(effort_selfreport[i]);
            yreal_selfreport[i] ~ lognormal(log_lambda_selfreport_i, sigma_logn[1]);
          }
        }
      }
    }

    // Family 6: Gamma
    else if(family == 6) {
      cv_gamma[1] ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0 && yreal_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            yreal_obs[i] ~ gamma(gammaA[1], gammaA[1] / lambda_obs_i);
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0 && yreal_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            yreal_em[i] ~ gamma(gammaA[1], gammaA[1] / lambda_em_i);
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0 && yreal_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            yreal_both[i] ~ gamma(gammaA[1], gammaA[1] / lambda_both_i);
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0 && yreal_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            yreal_selfreport[i] ~ gamma(gammaA[1], gammaA[1] / lambda_selfreport_i);
          }
        }
      }
    }

    // Family 7: Lognormal-Hurdle
    else if(family == 7) {
      sigma_logn ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real log_lambda_obs_i = log_lambda_base[time_idx_obs[i]] + log(effort_obs[i]);
            if (yint_obs[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_obs[i] > 0) {
                yreal_obs[i] ~ lognormal(log_lambda_obs_i, sigma_logn[1]);
              }
            }
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real log_lambda_em_i = log_lambda_base[time_idx_em[i]] + log(effort_em[i]);
            if (yint_em[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_em[i] > 0) {
                yreal_em[i] ~ lognormal(log_lambda_em_i, sigma_logn[1]);
              }
            }
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real log_lambda_both_i = log_lambda_base[time_idx_both[i]] + log(effort_both[i]);
            if (yint_both[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_both[i] > 0) {
                yreal_both[i] ~ lognormal(log_lambda_both_i, sigma_logn[1]);
              }
            }
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real log_lambda_selfreport_i = log_lambda_base[time_idx_selfreport[i]] + log(effort_selfreport[i]);
            if (yint_selfreport[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_selfreport[i] > 0) {
                yreal_selfreport[i] ~ lognormal(log_lambda_selfreport_i, sigma_logn[1]);
              }
            }
          }
        }
      }
    }

    // Family 8: Gamma-Hurdle
    else if(family == 8) {
      cv_gamma[1] ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            if (yint_obs[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_obs[i] > 0) {
                yreal_obs[i] ~ gamma(gammaA[1], gammaA[1] / lambda_obs_i);
              }
            }
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            if (yint_em[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_em[i] > 0) {
                yreal_em[i] ~ gamma(gammaA[1], gammaA[1] / lambda_em_i);
              }
            }
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            if (yint_both[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_both[i] > 0) {
                yreal_both[i] ~ gamma(gammaA[1], gammaA[1] / lambda_both_i);
              }
            }
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            if (yint_selfreport[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              if(yreal_selfreport[i] > 0) {
                yreal_selfreport[i] ~ gamma(gammaA[1], gammaA[1] / lambda_selfreport_i);
              }
            }
          }
        }
      }
    }

    // Family 9: Normal
    else if(family == 9) {
      sigma_logn ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            yreal_obs[i] ~ normal(lambda_obs_i, sigma_logn[1]);
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            yreal_em[i] ~ normal(lambda_em_i, sigma_logn[1]);
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            yreal_both[i] ~ normal(lambda_both_i, sigma_logn[1]);
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            yreal_selfreport[i] ~ normal(lambda_selfreport_i, sigma_logn[1]);
          }
        }
      }
    }

    // Family 10: Normal-Hurdle
    else if(family == 10) {
      sigma_logn ~ student_t(3, 0, 2);
      if(n_obs > 0) {
        for(i in 1:n_obs) {
          if(effort_obs[i] > 0) {
            real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
            if (yint_obs[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yreal_obs[i] ~ normal(lambda_obs_i, sigma_logn[1]);
            }
          }
        }
      }
      if(n_em > 0) {
        for(i in 1:n_em) {
          if(effort_em[i] > 0) {
            real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
            if (yint_em[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yreal_em[i] ~ normal(lambda_em_i, sigma_logn[1]);
            }
          }
        }
      }
      if(n_both > 0) {
        for(i in 1:n_both) {
          if(effort_both[i] > 0) {
            real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
            if (yint_both[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yreal_both[i] ~ normal(lambda_both_i, sigma_logn[1]);
            }
          }
        }
      }
      if(n_selfreport > 0) {
        for(i in 1:n_selfreport) {
          if(effort_selfreport[i] > 0) {
            real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
            if (yint_selfreport[i] == 0)
              1 ~ bernoulli(theta);
            else {
              0 ~ bernoulli(theta);
              yreal_selfreport[i] ~ normal(lambda_selfreport_i, sigma_logn[1]);
            }
          }
        }
      }
    }

    // ---- NEW: self-reporting probability sub-model likelihood ----
    // (Pooled/constant across years; ground truth from obs/em/both-covered vessels)
    if(estimate_self_report == 1) {
      logit_p_report ~ normal(0, 1.5); // weakly informative on the probability scale
      for(t in 1:n_year) {
        if(n_vessels_observed[t] > 0) {
          n_vessels_self_report[t] ~ binomial(n_vessels_observed[t], p_report[t]);
        }
      }
    }
  }
}
generated quantities {
  // Log-likelihood needs different sizes for pooled vs multi-stream models
  vector[multi_stream == 0 ? n_row : (n_obs + n_em + n_both + n_selfreport)] log_lik;
  int<lower = 0> y_new[n_year*is_discrete];
  vector[n_year*(1-is_discrete)] y_new_real;

  // ---- NEW: self-report-adjusted expansion effort ----
  // The self-report stream's effort may overlap with vessels already covered by
  // obs/em/both (a vessel can be both monitored AND self-report). p_report tells us,
  // per posterior draw, what fraction of MONITORED effort we'd expect to be
  // double-counted in the self-report stream. What's left over
  // (effort_selfreport_by_year - expected overlap) is genuinely NEW coverage, and
  // that portion is netted out of the truly-unaccounted effort before generating
  // y_new. Because this happens per draw, the width of p_report's posterior shows
  // up directly in the width of the expanded estimates.
  vector[n_year] effort_selfreport_new_coverage;
  vector[n_year] new_effort_adjusted;

  for(t in 1:n_year) {
    real expected_overlap = 0;
    if(estimate_self_report == 1) {
      expected_overlap = p_report[t] * effort_monitored_by_year[t];
    }
    effort_selfreport_new_coverage[t] = fmax(effort_selfreport_by_year[t] - expected_overlap, 0);
    new_effort_adjusted[t] = fmax(new_effort_by_year[t] - effort_selfreport_new_coverage[t], 0);
  }

  // POOLED MODEL (single stream)
  if(multi_stream == 0) {
    // Calculate pointwise log-likelihood for each observation
    for(n in 1:n_row) {
      real lambda_n = lambda_base[time[n]] * effort[n];
      real log_lambda_n = log_lambda_base[time[n]] + log(effort[n]);
      
      if(family == 1) {  // Poisson
        log_lik[n] = poisson_log_lpmf(yint[n] | log_lambda_n);
      } 
      else if(family == 2) {  // Negative Binomial
        log_lik[n] = neg_binomial_2_log_lpmf(yint[n] | log_lambda_n, nb2_phi[1]);
      } 
      else if(family == 3) {  // Poisson-Hurdle
        if(yint[n] == 0) {
          log_lik[n] = log(theta[1]);
        } else {
          log_lik[n] = log1m(theta[1]) + poisson_log_lpmf(yint[n] | log_lambda_n);
        }
      } 
      else if(family == 4) {  // NB-Hurdle
        if(yint[n] == 0) {
          log_lik[n] = log(theta[1]);
        } else {
          log_lik[n] = log1m(theta[1]) + neg_binomial_2_log_lpmf(yint[n] | log_lambda_n, nb2_phi[1]);
        }
      }
      else if(family == 5) {  // Lognormal
        log_lik[n] = lognormal_lpdf(yreal[n] | log_lambda_n, sigma_logn[1]);
      }
      else if(family == 6) {  // Gamma
        log_lik[n] = gamma_lpdf(yreal[n] | gammaA[1], gammaA[1] / lambda_n);
      }
      else if(family == 7) {  // Lognormal-Hurdle
        if(yreal[n] == 0) {
          log_lik[n] = log(theta[1]);
        } else {
          log_lik[n] = log1m(theta[1]) + lognormal_lpdf(yreal[n] | log_lambda_n, sigma_logn[1]);
        }
      }
      else if(family == 8) {  // Gamma-Hurdle
        if(yreal[n] == 0) {
          log_lik[n] = log(theta[1]);
        } else {
          log_lik[n] = log1m(theta[1]) + gamma_lpdf(yreal[n] | gammaA[1], gammaA[1] / lambda_n);
        }
      }
      else if(family == 9) {  // Normal
        log_lik[n] = normal_lpdf(yreal[n] | lambda_n, sigma_logn[1]);
      }
      else if(family == 10) {  // Normal-Hurdle
        if(yreal[n] == 0) {
          log_lik[n] = log(theta[1]);
        } else {
          log_lik[n] = log1m(theta[1]) + normal_lpdf(yreal[n] | lambda_n, sigma_logn[1]);
        }
      }
    }
  } 
  // MULTI-STREAM MODEL
  else {
    int idx = 1;  // Index for filling log_lik
    
    // OBS sector observations
    for(i in 1:n_obs) {
      real lambda_obs_i = lambda_base[time_idx_obs[i]] * effort_obs[i];
      real log_lambda_obs_i = log_lambda_base[time_idx_obs[i]] + log(effort_obs[i]);
      
      if(family == 1) {  // Poisson
        log_lik[idx] = poisson_log_lpmf(yint_obs[i] | log_lambda_obs_i);
      } else if(family == 2) {  // NB
        log_lik[idx] = neg_binomial_2_log_lpmf(yint_obs[i] | log_lambda_obs_i, nb2_phi[1]);
      } else if(family == 3) {  // Poisson-Hurdle
        if(yint_obs[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + poisson_log_lpmf(yint_obs[i] | log_lambda_obs_i);
        }
      } else if(family == 4) {  // NB-Hurdle
        if(yint_obs[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + neg_binomial_2_log_lpmf(yint_obs[i] | log_lambda_obs_i, nb2_phi[1]);
        }
      }
      idx += 1;
    }
    
    // EM sector observations
    for(i in 1:n_em) {
      real lambda_em_i = lambda_base[time_idx_em[i]] * effort_em[i];
      real log_lambda_em_i = log_lambda_base[time_idx_em[i]] + log(effort_em[i]);
      
      if(family == 1) {  // Poisson
        log_lik[idx] = poisson_log_lpmf(yint_em[i] | log_lambda_em_i);
      } else if(family == 2) {  // NB
        log_lik[idx] = neg_binomial_2_log_lpmf(yint_em[i] | log_lambda_em_i, nb2_phi[1]);
      } else if(family == 3) {  // Poisson-Hurdle
        if(yint_em[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + poisson_log_lpmf(yint_em[i] | log_lambda_em_i);
        }
      } else if(family == 4) {  // NB-Hurdle
        if(yint_em[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + neg_binomial_2_log_lpmf(yint_em[i] | log_lambda_em_i, nb2_phi[1]);
        }
      }
      idx += 1;
    }
    
    // BOTH sector observations
    for(i in 1:n_both) {
      real lambda_both_i = lambda_base[time_idx_both[i]] * effort_both[i];
      real log_lambda_both_i = log_lambda_base[time_idx_both[i]] + log(effort_both[i]);
      
      if(family == 1) {  // Poisson
        log_lik[idx] = poisson_log_lpmf(yint_both[i] | log_lambda_both_i);
      } else if(family == 2) {  // NB
        log_lik[idx] = neg_binomial_2_log_lpmf(yint_both[i] | log_lambda_both_i, nb2_phi[1]);
      } else if(family == 3) {  // Poisson-Hurdle
        if(yint_both[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + poisson_log_lpmf(yint_both[i] | log_lambda_both_i);
        }
      } else if(family == 4) {  // NB-Hurdle
        if(yint_both[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + neg_binomial_2_log_lpmf(yint_both[i] | log_lambda_both_i, nb2_phi[1]);
        }
      }
      idx += 1;
    }

    // SELFREPORT sector observations
    for(i in 1:n_selfreport) {
      real lambda_selfreport_i = lambda_base[time_idx_selfreport[i]] * effort_selfreport[i];
      real log_lambda_selfreport_i = log_lambda_base[time_idx_selfreport[i]] + log(effort_selfreport[i]);
      
      if(family == 1) {  // Poisson
        log_lik[idx] = poisson_log_lpmf(yint_selfreport[i] | log_lambda_selfreport_i);
      } else if(family == 2) {  // NB
        log_lik[idx] = neg_binomial_2_log_lpmf(yint_selfreport[i] | log_lambda_selfreport_i, nb2_phi[1]);
      } else if(family == 3) {  // Poisson-Hurdle
        if(yint_selfreport[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + poisson_log_lpmf(yint_selfreport[i] | log_lambda_selfreport_i);
        }
      } else if(family == 4) {  // NB-Hurdle
        if(yint_selfreport[i] == 0) {
          log_lik[idx] = log(theta[1]);
        } else {
          log_lik[idx] = log1m(theta[1]) + neg_binomial_2_log_lpmf(yint_selfreport[i] | log_lambda_selfreport_i, nb2_phi[1]);
        }
      }
      idx += 1;
    }
  }
  
  // Generate posterior predictive samples for unobserved effort (by year)
  // NOTE: uses new_effort_adjusted (nets out self-report-inferred coverage)
  // instead of the raw new_effort_by_year.
  for(t in 1:n_year) {
    if(is_discrete == 1) {
      y_new[t] = 0;
      if(new_effort_adjusted[t] > 0) {
        real lambda_new_t = lambda_base[t] * new_effort_adjusted[t];
        
        if(family == 1) {
          y_new[t] = poisson_rng(lambda_new_t);
        } else if(family == 2) {
          y_new[t] = neg_binomial_2_rng(lambda_new_t, nb2_phi[1]);
        } else if(family == 3) {
          y_new[t] = (1 - bernoulli_rng(theta[1])) * poisson_rng(lambda_new_t);
        } else if(family == 4) {
          y_new[t] = (1 - bernoulli_rng(theta[1])) * neg_binomial_2_rng(lambda_new_t, nb2_phi[1]);
        }
      }
    } else {
      y_new_real[t] = 0;
      if(new_effort_adjusted[t] > 0) {
        real log_lambda_new_t = log_lambda_base[t] + log(new_effort_adjusted[t]);
        real lambda_new_t = lambda_base[t] * new_effort_adjusted[t];
        
        if(family == 5) {
          y_new_real[t] = lognormal_rng(log_lambda_new_t, sigma_logn[1]);
        } else if(family == 6) {
          y_new_real[t] = gamma_rng(gammaA[1], gammaA[1] / lambda_new_t);
        } else if(family == 7) {
          y_new_real[t] = (1 - bernoulli_rng(theta[1])) * lognormal_rng(log_lambda_new_t, sigma_logn[1]);
        } else if(family == 8) {
          y_new_real[t] = (1 - bernoulli_rng(theta[1])) * gamma_rng(gammaA[1], gammaA[1] / lambda_new_t);
        } else if(family == 9) {
          y_new_real[t] = normal_rng(lambda_new_t, sigma_logn[1]);
        } else if(family == 10) {
          y_new_real[t] = (1 - bernoulli_rng(theta[1])) * normal_rng(lambda_new_t, sigma_logn[1]);
        }
      }
    }
  }
}
