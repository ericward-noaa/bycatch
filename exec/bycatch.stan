data {
  int<lower=0> n_year;
  int<lower=0> K;
  matrix[n_year, K] x; # covariates
  vector[n_year] effort; # covariates
  int events[n_year]; # vector[n_pos] y; # data
  int family;
  int time_varying; # whether to treat model as dlm
}
parameters {
  vector[K] beta;
  vector[(n_year-1)] est_time_dev;
  real<lower=0> sigma_rw;
}
transformed parameters {
  vector[n_year] lambda;
  vector[n_year] pred;
  vector[n_year] time_dev;
  pred = x * beta;
  time_dev[1] = 0;
  for(i in 2:n_year) {
    time_dev[i] = est_time_dev[i-1];
  }
  for(i in 1:n_year) {
    pred[i] = pred[i] + time_varying*time_dev[i];
    lambda[i] = exp(pred[i]) * effort[i];
  }
}
model {
  beta ~ student_t(3, 0, 3);
  sigma_rw ~ cauchy(0,5);
  est_time_dev[i] ~ student_t(3, 0, 3);
  for(i in 2:(n_year-1)) {
    # model evolution of temporal deviations as random walk
    est_time_dev[i] ~ normal(est_time_dev[i], sigma_rw);
  }
  if(family == 1) {
    events ~ poisson(lambda);
  }

}
