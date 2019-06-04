#' bycatch_expansion is the primary function for fitting bycatch models to time series of takes and effort
#' @param formula The model formula.
#' @param data A data frame.
#' @param time Named column of the 'data' data frame with the label for the time (e.g. year) variable
#' @param effort Named column of the 'effort' variable in the data frame with the label for the fishing effort to be used in estimation of mean bycatch rate
#' @param coverage Named column of the 'observer coverage' variable in the data frame (0 - 100) used for binomial expansion
#' @param family "poisson" or "nbinom2", defaults to "poisson". Both default to using a log link function
#' @param time_varying boolean TRUE/FALSE, whether to include time varying component (this is a random walk, analogous to making this a Dynamic linear model)
#' @param control List to pass to [rstan::sampling()]. For example,
#'   increase \code{adapt_delta} if there are warnings about divergent
#'   transitions: \code{control = list(adapt_delta = 0.99)}. By default,
#'   \pkg{glmmfields} sets \code{adapt_delta = 0.9}.
#'
#' @return list of the data used to fit the model, the matrix of covariates, the expanded bycatch generated via the fit and simulations, and the fitted stan model
#'
#' @export
#' @importFrom rstan sampling vb
#' @import Rcpp
#' @importFrom stats model.frame model.matrix model.response
#'
fit_bycatch <- function(formula, data, time = "year", effort = "effort",
  family = c("poisson","nbinom2"), time_varying = FALSE) {

  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf, "numeric")

  if(time %in% colnames(data) == FALSE) {
    stop("The time variable needs to be specified as a named column in the data frame")
  }
  if(effort %in% colnames(data) == FALSE) {
    stop("The effort variable needs to be specified as a named column in the data frame")
  }
  if(coverage %in% colnames(data) == FALSE) {
    stop("The coverage variable needs to be specified as a named column in the data frame")
  }
  if(family %in% c("poisson","nbinom2") == FALSE) {
    stop("The family must be specified as poisson or nbinom2")
  }
  family = ifelse(family=="poisson",1,2)

  pars = c("beta", "lambda", "log_lik")
  datalist = list(n_year = nrow(df),
    effort = data[,effort],
    events = y,
    time = data[,time],
    x = X,
    K = ncol(X),
    family = 1,
    time_varying = as.numeric(time_varying))

  sampling_args <- list(
    object = stanmodels$bycatch,
    data = stan_data,
    pars = pars,
    control = control, ...
  )

  pars = rstan::extract(stan_model)
  df$lambda = apply(pars$lambda, 2, mean)

  # do binomial expansion of estimated takes
  # estimated takes can be thought of as the observed draw from binomial distribution where p (but not N) is known
  binom_p = df$coverage/100
  # point estimates may be unstable - see Raftery (1988) - so we can simulate distribution

  # sample ~ 1000 random draws from posterior
  sigfig_multiplier = control$sigfig_multiplier # default to 100
  mcmc_samples = control$mcmc_samples # default to 1000
  maxX = control$maxX # default to 20000

  samples = sample(1:dim(pars$lambda)[1], size=mcmc_samples, replace=F)
  expanded_estimates = matrix(0, length(samples), nrow(df))

  for(y in 1:nrow(df)) {
    for(mcmc in 1:length(samples)) {
      # We observe takes, df$events[y]
      # pmf of N | p
      # In other words, we need to sample from the density of Binomial N
      # given Binomial p and X

      # X needs to be an integer, so expand to a series of large numbers as approximation
      X = round(pars$lambda[samples[mcmc],y] * sigfig_multiplier)
      # calculate the probabilites of these Ns given the mean observed takes and observer coverage
      prob_N = dbinom(x=X, size = seq(X, maxX), prob = binom_p[y])
      # sample from distribution of N, representing total mean takes
      N = sample(seq(X, maxX),size=1, prob=prob_N)
      # sample from poisson to convert mean -> observed data with observation model
      N = rpois(n = 1, lambda = N)
      # use the expansion for unobserved sets, use observed as known perfectly for observed
      N = (1-binom_p[y]) * N / sigfig_multiplier + df$events[y]
      expanded_estimates[mcmc,y] = N
    }
  }

  return(list("data" = data, "effort"=effort, "time"=time))
}
