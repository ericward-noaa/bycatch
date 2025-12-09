#' get_total is a helper function to return a matrix of total estimated bycatch
#'
#' @param fitted_model Data and fitted model returned from estimation
#'
#' @return matrix (MCMC draws x time steps) of posterior predictive values for total bycatch (observed + unobserved)
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
#' total <- get_total(fit)
#'
#' # Calculate total bycatch summaries
#' total_mean <- colMeans(total)
#' total_quantiles <- apply(total, 2, quantile, probs = c(0.025, 0.975))
#' }
get_total <- function(fitted_model) {
  # Get expanded estimates for unobserved effort
  expanded_estimates <- get_expanded(fitted_model)

  # Get observed events
  observed_events <- fitted_model$data[[fitted_model$events]]

  # Vectorized addition: add observed events to each MCMC draw
  # sweep() applies the operation across columns (MARGIN = 2)
  # This is much faster than looping through columns
  total_estimates <- sweep(expanded_estimates, 2, observed_events, "+")

  return(total_estimates)
}
