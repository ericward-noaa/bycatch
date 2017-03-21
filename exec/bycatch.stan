data {
  int<lower=0> n_year;
  int<lower=0> K;
  matrix[n_year, K] x; # covariates
  vector[n_year] effort; # covariates
  int events[n_year]; # vector[n_pos] y; # data
  int family;
}
parameters {
  vector[K] beta;
}
transformed parameters {
  vector[n_year] lambda;
  vector[n_year] pred;
  pred = x * beta;
  for(i in 1:n_year) {
    lambda[i] = exp(pred[i]) * effort[i];
  }
}
model {
  beta ~ student_t(3, 0, 3);

  if(family == 1) {
    events ~ poisson(lambda);
  }

}
