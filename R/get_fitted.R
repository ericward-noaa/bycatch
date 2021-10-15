#' get_fitted returns df of observed bycatch estimates (lambda of Poisson), accounting for effort but not accounting for observer coverage
#'
#' @param fitted_model Data and fitted model returned from fit_bycatch(). If a hurdle model, then only
#' then the plot returns the total bycatch rate (including zero and non-zero components).
#' @param alpha The alpha level for the credible interval, defaults to 0.05
#' @return plot called from ggplot
#'
#' @export
#' @import ggplot2
#' @importFrom stats quantile
#' @examples
#' \donttest{
#' d <- data.frame(
#'   "Year" = 2002:2014,
#'   "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
#'   "expansionRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
#'   "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
#' )
#' fit <- fit_bycatch(Takes ~ 1,
#'   data = d, time = "Year", effort = "Sets",
#'   family = "poisson", time_varying = FALSE
#' )
#' get_fitted(fit)
#' }
get_fitted <- function(fitted_model, alpha = 0.05) {
  lambda <- rstan::extract(fitted_model$fitted_model, c("lambda"))$lambda

  if (fitted_model$family %in% c("poisson-hurdle", "nbinom2-hurdle", "gamma-hurdle", "lognormal-hurdle", "normal-hurdle")) {
    # adjust lambda estimates by including theta. theta = pr(0), (1-theta) = pr(>0)
    theta <- rstan::extract(fitted_model$fitted_model, c("theta"))$theta
    for (i in 1:nrow(lambda)) lambda[i, ] <- lambda[i, ] * (1 - theta[i, 1])
  }
  df <- data.frame(
    "time" = fitted_model$data[, fitted_model$time],
    "mean" = apply(lambda, 2, mean),
    "low" = apply(lambda, 2, quantile, alpha/2),
    "high" = apply(lambda, 2, quantile, 1-alpha/2),
    "obs" = fitted_model$data[, fitted_model$events]
  )

  return(df)
}
