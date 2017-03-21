#' bycatch_expansion is the primary function for fitting bycatch models to time series of takes and effort
#'
#' @param time Numeric, e.g. years
#' @param events Integer vector of takes (events)
#' @param covar Optional matrix of covariates
#' @param effort Metric of fishing effort to be used in estimation of mean bycatch rate
#' @param coverage Observer coverage (0 - 100) used for binomial expansion
#' @param family Observation error distribution, defaults to Poisson
#'
#' @return list of the data used to fit the model, the matrix of covariates, the expanded bycatch generated via the fit and simulations, and the fitted stan model
#'
#' @export
bycatch_expansion <- function(time = NULL, events = NULL, covar = NULL, effort = NULL, coverage = NULL, family = c("poisson")) {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  df = data.frame(time=time, events=events, effort=effort, coverage=coverage)

  stan_dir = find.package("bycatch")
  model = paste0(stan_dir, "/exec/bycatch.stan")

  if(is.null(covar)) {
    # covariate matrix =
    covar = matrix(1, nrow(df), ncol=1)
  }

  pars = c("beta", "lambda")
  datalist = list(n_year = nrow(df),
    effort = df$effort,
    events = df$events,
    x = covar,
    K = ncol(covar),
    family = 1)

  # Currently each year as modeled as independent
  if(family == "poisson") {
    # intercept dropped, because lambda = theta * effort
    stan_model = stan(file=model, data = datalist, pars = pars, chains=4, iter=2000)
  }

  pars = rstan::extract(stan_model)
  df$lambda = apply(pars$lambda, 2, mean)

  # do binomial expansion of estimated takes
  # estimated takes can be thought of as the observed draw from binomial distribution where p (but not N) is known
  binom_p = df$coverage/100
  # point estimates may be unstable - see Raftery (1988) - so we can simulate distribution

  # sample ~ 1000 random draws from posterior
  sigfig_multiplier = 1000
  mcmc_samples = 1000
  samples = sample(1:dim(pars$lambda)[1], size=mcmc_samples, replace=F)
  expanded_estimates = matrix(0, length(samples), nrow(df))

  for(y in 1:nrow(df)) {
    for(mcmc in 1:length(samples)) {
      # pmf of N | p
      X = round(pars$lambda[samples[mcmc],y] * sigfig_multiplier)
      probs = dbinom(x=X,
        size = seq(X, 20000), prob = binom_p[y])
      expanded_estimates[mcmc,y] = sample(seq(X, 20000),size=1, prob=probs)
    }
  }
  expanded_estimates = expanded_estimates / sigfig_multiplier

  return(list("data" = df, "covar"=covar, "expanded_estimates" = expanded_estimates, "fitted_model" = stan_model))
}
