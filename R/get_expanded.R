#' get_expanded is a helper function to return a matrix of posterior predictive values for unobserved bycatch
#'
#' @param fitted_model Data and fitted model returned from estimation

#' @return matrix (MCMC draws x time steps) of posterior predictive values for unobserved bycatch
#'
#' @export
#' @importFrom rstan extract
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
#'   effort = "Sets",
#'   family = "poisson",
#'   expansion_rate = "expansionRate",
#'   time_varying = FALSE
#' )
#' expanded <- get_expanded(fit)
#' }
get_expanded <- function(fitted_model) {
  expanded_estimates <- rstan::extract(fitted_model$fitted_model, "y_new")$y_new
  return(expanded_estimates)
}
