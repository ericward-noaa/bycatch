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

  # Check if multi-stream mode
  if(fitted_model$multi_stream) {
    # Multi-stream: expanded_estimates is by year
    # Calculate observed takes by year
    time_values <- unique(fitted_model$data[[fitted_model$time]])
    time_values <- sort(time_values)
    n_year <- length(time_values)

    obs_by_year <- rep(0, n_year)

    # Add Observer takes
    for(i in 1:nrow(fitted_model$data)) {
      t_idx <- which(time_values == fitted_model$data[[fitted_model$time]][i])
      obs_by_year[t_idx] <- obs_by_year[t_idx] + fitted_model$data[[fitted_model$events]][i]
    }

    # Add EM takes if present
    if(!is.null(fitted_model$stream_info$takes_em)) {
      for(i in 1:nrow(fitted_model$data)) {
        t_idx <- which(time_values == fitted_model$data[[fitted_model$time]][i])
        obs_by_year[t_idx] <- obs_by_year[t_idx] + fitted_model$data[[fitted_model$stream_info$takes_em]][i]
      }
    }

    # Add Both takes if present
    if(!is.null(fitted_model$stream_info$takes_both)) {
      for(i in 1:nrow(fitted_model$data)) {
        t_idx <- which(time_values == fitted_model$data[[fitted_model$time]][i])
        obs_by_year[t_idx] <- obs_by_year[t_idx] + fitted_model$data[[fitted_model$stream_info$takes_both]][i]
      }
    }

    # Add observed to expanded (vectorized)
    total_estimates <- sweep(expanded_estimates, 2, obs_by_year, "+")

  } else {
    # Single-stream: expanded_estimates is by observation
    # Get observed events
    observed_events <- fitted_model$data[[fitted_model$events]]

    # Vectorized addition using sweep()
    total_estimates <- sweep(expanded_estimates, 2, observed_events, "+")
  }

  return(total_estimates)
}
