data {
  int<lower=0> n_year;
  real effort[n_year]; # vector[n_pos] y; # data
  int events[n_year]; # vector[n_pos] y; # data
  int family;
}
parameters {
  real<lower=0> theta; # estimated coefficient for effort
}
transformed parameters {
  vector<lower=0>[n_year] lambda;
  for(i in 1:n_year) {
    lambda[i] = theta * effort[i];
  }
}
model {
  theta ~ cauchy(0,5);

  if(family == 1) {
    events ~ poisson(lambda);
  }

}
