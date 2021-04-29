#' fit_bycatch is the primary function for fitting bycatch models to time series of takes and effort
#' @param formula The model formula.
#' @param data A data frame.
#' @param time Named column of the 'data' data frame with the label for the time (e.g. year) variable
#' @param effort Named column of the 'effort' variable in the data frame with the label for the fishing effort to be used in estimation of mean bycatch rate
#' @param family Family for response distribution ("poisson", "nbinom2", "poisson-hurdle","nbinom2-hurdle"), defaults to "poisson". The hurdle variants estimate
#' the probability of zeros (theta) separately from the other models and use truncated distribution to model positive counts. All use a log link function
#' @param time_varying boolean TRUE/FALSE, whether to include time varying component (this is a random walk, analogous to making this a Dynamic linear model)
#' @param iter the number of mcmc iterations, defaults to 1000
#' @param chains the number of mcmc chains, defaults to 3
#' @param control List to pass to [rstan::sampling()]. For example,
#'   increase \code{adapt_delta} if there are warnings about divergent
#'   transitions: \code{control = list(adapt_delta = 0.99)}. By default,
#'   \pkg{glmmfields} sets \code{adapt_delta = 0.9}.
#' @param ... Any other arguments to pass to [rstan::sampling()].
#'
#' @return list of the data used to fit the model, the matrix of covariates, the expanded bycatch generated via the fit and simulations, and the fitted stan model
#'
#' @export
#'
#' @importFrom rstan sampling vb
#' @import Rcpp
#' @importFrom stats model.frame model.matrix model.response
#'
#' @examples
#' \donttest{
#' d <- data.frame(
#'   "Year" = 2002:2014,
#'   "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
#'   "expansionRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
#'   "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
#' )
#' fit <- fit_bycatch(Takes ~ 1,
#'   data = d, time = "Year",
#'   effort = "Sets", family = "poisson", time_varying = FALSE
#' )
#' loo::loo(fit$fitted_model)$estimates
#' }
fit_bycatch <- function(formula, data, time = "year", effort = "effort",
                        family = c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle"), time_varying = FALSE,
                        iter = 1000, chains = 3, control = list(adapt_delta = 0.9, max_treedepth = 20), ...) {
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- as.numeric(model.response(mf, "numeric"))

  if (time %in% colnames(data) == FALSE) {
    stop("The time variable needs to be specified as a named column in the data frame")
  }
  if (effort %in% colnames(data) == FALSE) {
    stop("The effort variable needs to be specified as a named column in the data frame")
  }
  if (family %in% c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle") == FALSE) {
    stop("The family must be specified as poisson, nbinom2, poisson-hurdle, or nbinom2-hurdle")
  }

  pars <- c("beta", "lambda", "log_lik")
  if (family != "poisson") pars <- c(pars, "nb2_phi")
  if (family %in% c("poisson-hurdle", "nbinom2-hurdle")) pars <- c(pars, "theta")

  family_id <- match(family, c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle")) # 1, 2, 3, 4
  datalist <- list(
    n_row = nrow(data),
    effort = data[, effort],
    yint = y,
    time = data[, time] - min(data[, time]) + 1,
    n_year = length(seq(min(data[, time]), max(data[, time]))),
    K = ncol(X),
    x = X,
    family = family_id,
    time_varying = as.numeric(time_varying)
  )

  mod <- rstan::sampling(
    object = stanmodels$bycatch,
    data = datalist,
    pars = pars,
    iter = iter,
    chains = chains,
    control = control, ...
  )

  return(list("data" = data, "effort" = effort, "events" = names(mf)[1], "time" = time, "fitted_model" = mod, "family" = family))
}
