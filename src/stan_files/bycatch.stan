data {
  int<lower=0> n_row;
  vector[n_row] effort; // covariates
  int yint[n_row]; // vector[n_pos] y; # data
  int time[n_row]; // time variable
  int<lower=0> n_year; // number of unique years
  int<lower=0> K;
  matrix[n_row, K] x; // covariates
  int family; // 0 = poisson, 1 = negbin
  int time_varying; // whether to treat model as dlm
}
parameters {
  vector[K] beta;
  vector[time_varying*(n_year-1)] est_time_dev;
  real<lower=0> sigma_rw[time_varying];
  real<lower=0> nb2_phi[family];
}
transformed parameters {
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
    lambda[i] = exp(pred[i] + log(effort[i])); // exp(pred) = theta
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

  if(family == 0) {
    for(i in 1:n_row) {
      yint[i] ~ poisson(lambda[i]);
    }
  } else {
    nb2_phi ~ student_t(3, 0, 2);
    for(i in 1:n_row) {
      yint[i] ~ neg_binomial_2(lambda[i], nb2_phi[1]);
    }
  }
}
generated quantities {
  vector[n_row] log_lik;
  if(family==0) {
    for(n in 1:n_row) {
      log_lik[n] = poisson_lpmf(yint[n] | lambda[n]);
    }
  } else {
    for(n in 1:n_row) {
      log_lik[n] = neg_binomial_2_lpmf(yint[n] | lambda[n], nb2_phi[1]);
    }
  }
}
