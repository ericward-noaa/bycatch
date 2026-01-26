#' get_fitted returns df of observed bycatch estimates (lambda of Poisson), accounting for effort but not accounting for observer coverage
#'
#' @param fitted_model Data and fitted model returned from fit_bycatch(). If a hurdle model, then the plot returns the total bycatch
#' rate (including zero and non-zero components).
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

  # Check if multi-stream mode
  if(fitted_model$multi_stream) {
    # Multi-stream: Extract lambda_base (per year)
    lambda_base <- rstan::extract(fitted_model$fitted_model, "lambda_base")$lambda_base

    # Adjust for hurdle models
    if (fitted_model$family %in% c("poisson-hurdle", "nbinom2-hurdle", "gamma-hurdle",
                                   "lognormal-hurdle", "normal-hurdle")) {
      theta <- rstan::extract(fitted_model$fitted_model, "theta")$theta
      # lambda_base already represents the rate, multiply by (1-theta) for hurdle
      for (i in 1:nrow(lambda_base)) {
        lambda_base[i, ] <- lambda_base[i, ] * (1 - theta[i, 1])
      }
    }

    # Get unique time values
    time_values <- unique(fitted_model$data[[fitted_model$time]])
    time_values <- sort(time_values)

    # Create data frame by year
    df <- data.frame(
      "time" = time_values,
      "mean" = apply(lambda_base, 2, mean),
      "low" = apply(lambda_base, 2, quantile, alpha/2),
      "high" = apply(lambda_base, 2, quantile, 1-alpha/2)
    )

    # Add observed takes (sum across streams for each year)
    obs_by_year <- rep(0, length(time_values))
    for(i in 1:nrow(fitted_model$data)) {
      t_idx <- which(time_values == fitted_model$data[[fitted_model$time]][i])
      obs_by_year[t_idx] <- obs_by_year[t_idx] + fitted_model$data[[fitted_model$events]][i]
    }

    # Add EM if present
    if(!is.null(fitted_model$stream_info$takes_em)) {
      for(i in 1:nrow(fitted_model$data)) {
        t_idx <- which(time_values == fitted_model$data[[fitted_model$time]][i])
        obs_by_year[t_idx] <- obs_by_year[t_idx] + fitted_model$data[[fitted_model$stream_info$takes_em]][i]
      }
    }

    # Add Both if present
    if(!is.null(fitted_model$stream_info$takes_both)) {
      for(i in 1:nrow(fitted_model$data)) {
        t_idx <- which(time_values == fitted_model$data[[fitted_model$time]][i])
        obs_by_year[t_idx] <- obs_by_year[t_idx] + fitted_model$data[[fitted_model$stream_info$takes_both]][i]
      }
    }

    df$obs <- obs_by_year

  } else {
    # Single-stream mode (original code)
    lambda <- rstan::extract(fitted_model$fitted_model, "lambda")$lambda

    if (fitted_model$family %in% c("poisson-hurdle", "nbinom2-hurdle", "gamma-hurdle",
                                   "lognormal-hurdle", "normal-hurdle")) {
      # Adjust lambda estimates by including theta
      theta <- rstan::extract(fitted_model$fitted_model, "theta")$theta
      lambda <- lambda * (1 - theta[, 1])
    }

    df <- data.frame(
      "time" = fitted_model$data[[fitted_model$time]],
      "mean" = apply(lambda, 2, mean),
      "low" = apply(lambda, 2, quantile, alpha/2),
      "high" = apply(lambda, 2, quantile, 1-alpha/2),
      "obs" = fitted_model$data[[fitted_model$events]]
    )
  }

  return(df)
}
