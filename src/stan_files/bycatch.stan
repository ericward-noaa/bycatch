data {
  int<lower=0> n_row;
  vector[n_row] effort; // covariates
  vector[n_row] new_effort; // covariate for unobserved sets
  int yint[n_row]; // vector[n_pos] y; # data
  int time[n_row]; // time variable
  int<lower=0> n_year; // number of unique years
  int<lower=0> K;
  matrix[n_row, K] x; // covariates
  int family; // 1 = poisson, 2 = negbin, 3 = hurdle-poisson, 4 = hurdle-negbin
  int time_varying; // whether to treat model as dlm
}
transformed data {
  int est_phi;
  int est_theta;
  est_phi = 0;
  est_theta = 0;
  if(family == 2) est_phi=1;
  if(family == 3) est_theta=1;
  if(family == 4) {
    est_phi=1;
    est_theta=1;
  }
}
parameters {
  vector[K] beta;
  vector[time_varying*(n_year-1)] est_time_dev;
  real<lower=0> sigma_rw[time_varying];
  real<lower=0> nb2_phi[est_phi];
  real<lower=0,upper=1> theta[est_theta];
}
transformed parameters {
  vector[n_row] log_lambda;
  vector[n_row] lambda;
  vector[n_row] pred;
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
    theta ~ beta(1,1);
  }

  if(family == 1) {
    yint ~ poisson_log(log_lambda);
  }
  if(family == 2) {
    nb2_phi ~ student_t(3, 0, 2);
    yint ~ neg_binomial_2_log(log_lambda, nb2_phi[1]);
  }
  if(family == 3) {
    for(i in 1:n_row) {
      // this is a slow way, but ok for small N. Otherwise we'd want to switch to sufficient stats
      // https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
      if (yint[i] == 0)
        1 ~ bernoulli(theta);
      else {
        0 ~ bernoulli(theta);
        yint[i] ~ poisson(lambda[i]) T[1, ];
      }
    }
  }
  if(family == 4) {
    for(i in 1:n_row) {
      if (yint[i] == 0)
        1 ~ bernoulli(theta);
      else {
        0 ~ bernoulli(theta);
        yint[i] ~ neg_binomial_2(lambda[i], nb2_phi[1]) T[1, ];
      }
    }
  }
}
generated quantities {
  vector[n_row] log_lik;
  int<lower = 0> y_new[n_row];

  if(family==1) {
    for(n in 1:n_row) {
      log_lik[n] = poisson_log_lpmf(yint[n] | log_lambda[n]);

      // sample posterior predictive distribution
      y_new[n] = poisson_log_rng(pred[n] + log(new_effort[n]));
    }
  }
  if(family==2) {
    for(n in 1:n_row) {
      log_lik[n] = neg_binomial_2_log_lpmf(yint[n] | log_lambda[n], nb2_phi[1]);

      // sample posterior predictive distribution
      y_new[n] = neg_binomial_2_log_rng(pred[n] + log(new_effort[n]), nb2_phi[1]);
    }
  }
  if(family==3) {
    for(n in 1:n_row) {
      if(yint[n]==0) {
        // 1-theta because we want to calculate pr(0)
        log_lik[n] = log(theta[1]);
      } else {
        // (1 - theta) * Pr(Pois(y|lambda)) / (1 - PoissonCDF(0|lambda))
        log_lik[n] = log1m(theta[1]) + poisson_log_lpmf(yint[n] | log_lambda[n]) - poisson_lccdf(0 | lambda[n]);
      }

      // sample posterior predictive distribution
      y_new[n] = (1 - bernoulli_rng(theta[1])) * poisson_log_rng(pred[n] + log(new_effort[n]));
    }
  }
  if(family==4) {
    for(n in 1:n_row) {
      if(yint[n]==0) {
        // 1-theta because we want to calculate pr(0)
        log_lik[n] = log(theta[1]);
      } else {
        // (1 - theta) * Pr(NB2(y|lambda)) / (1 - NB2CDF(0|lambda))
        log_lik[n] = log1m(theta[1]) + neg_binomial_2_log_lpmf(yint[n] | log_lambda[n], nb2_phi[1]) - neg_binomial_2_lccdf(0 | lambda[n], nb2_phi[1]);
      }

      // sample posterior predictive distribution
      y_new[n] = (1 - bernoulli_rng(theta[1])) * neg_binomial_2_log_rng(pred[n] + log(new_effort[n]), nb2_phi[1]);
    }
  }
}
