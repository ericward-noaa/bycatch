#' get_stream_summary provides a summary of monitoring coverage and bycatch by stream
#'
#' @param fitted_model Data and fitted model returned from fit_bycatch() with multi-stream monitoring
#' @param alpha The alpha level for the credible interval, defaults to 0.05
#'
#' @return data frame with summary statistics for each monitoring stream
#'
#' @export
#' @importFrom rstan extract
#' @importFrom stats quantile
#'
#' @examples
#' \donttest{
#' d <- data.frame(
#'   Year = 2002:2014,
#'   Takes_obs = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0),
#'   Sets_obs = c(200, 180, 150, 350, 250, 270, 180, 150, 400, 350, 280, 180, 250),
#'   Takes_em = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
#'   Sets_em = c(150, 120, 140, 250, 180, 190, 120, 100, 300, 270, 200, 130, 190),
#'   Sets_total = c(1000, 950, 900, 1800, 1200, 1300, 900, 750, 2000, 1800, 1400, 950, 1300)
#' )
#'
#' fit <- fit_bycatch(Takes_obs ~ 1,
#'   data = d, time = "Year",
#'   effort = "Sets_obs",
#'   takes_em = "Takes_em",
#'   effort_em = "Sets_em",
#'   effort_total = "Sets_total",
#'   family = "poisson"
#' )
#'
#' get_stream_summary(fit)
#' }
get_stream_summary <- function(fitted_model, alpha = 0.05) {

  if (!fitted_model$multi_stream) {
    warning("Model was fit in single-stream mode. Returning basic summary.")

    # Single stream summary
    total_effort <- sum(fitted_model$data[[fitted_model$effort]])
    total_takes <- sum(fitted_model$data[[fitted_model$events]])

    # Get expanded estimates
    expanded <- get_expanded(fitted_model)
    expanded_mean <- mean(rowSums(expanded))
    expanded_low <- quantile(rowSums(expanded), alpha / 2)
    expanded_high <- quantile(rowSums(expanded), 1 - alpha / 2)

    # Get total estimates
    total <- get_total(fitted_model)
    total_mean <- mean(rowSums(total))
    total_low <- quantile(rowSums(total), alpha / 2)
    total_high <- quantile(rowSums(total), 1 - alpha / 2)

    summary_df <- data.frame(
      stream = c("Observed", "Unobserved", "Total"),
      effort = c(total_effort, sum(fitted_model$new_effort),
                 total_effort + sum(fitted_model$new_effort)),
      observed_takes = c(total_takes, NA, total_takes),
      estimated_mean = c(total_takes, expanded_mean, total_mean),
      estimated_low = c(total_takes, expanded_low, total_low),
      estimated_high = c(total_takes, expanded_high, total_high),
      stringsAsFactors = FALSE
    )

    return(summary_df)
  }

  # Multi-stream mode
  stream_info <- fitted_model$stream_info

  # Extract effort and takes for each stream
  effort_obs <- sum(fitted_model$data[[fitted_model$effort]])
  takes_obs <- sum(fitted_model$data[[fitted_model$events]])

  # Initialize for optional streams
  effort_em <- 0
  takes_em <- 0
  effort_both <- 0
  takes_both <- 0

  streams <- c("Observer")
  efforts <- c(effort_obs)
  takes <- c(takes_obs)

  if (!is.null(stream_info$takes_em)) {
    effort_em <- sum(fitted_model$data[[stream_info$effort_em]])
    takes_em <- sum(fitted_model$data[[stream_info$takes_em]])
    streams <- c(streams, "EM")
    efforts <- c(efforts, effort_em)
    takes <- c(takes, takes_em)
  }

  if (!is.null(stream_info$takes_both)) {
    effort_both <- sum(fitted_model$data[[stream_info$effort_both]])
    takes_both <- sum(fitted_model$data[[stream_info$takes_both]])
    streams <- c(streams, "Both")
    efforts <- c(efforts, effort_both)
    takes <- c(takes, takes_both)
  }

  # Pooled observed
  effort_pooled <- effort_obs + effort_em + effort_both
  takes_pooled <- takes_obs + takes_em + takes_both
  streams <- c(streams, "Pooled Observed")
  efforts <- c(efforts, effort_pooled)
  takes <- c(takes, takes_pooled)

  # Get expanded estimates (unobserved)
  expanded <- get_expanded(fitted_model)
  expanded_mean <- mean(rowSums(expanded))
  expanded_low <- quantile(rowSums(expanded), alpha / 2)
  expanded_high <- quantile(rowSums(expanded), 1 - alpha / 2)

  # Get total effort (from data if available)
  if (!is.null(stream_info$effort_total)) {
    effort_total <- sum(fitted_model$data[[stream_info$effort_total]])
  } else {
    effort_total <- effort_pooled  # Assume 100% coverage
  }
  effort_unobs <- effort_total - effort_pooled

  streams <- c(streams, "Unobserved")
  efforts <- c(efforts, effort_unobs)
  takes <- c(takes, NA)

  # Get total estimates
  total <- get_total(fitted_model)
  total_mean <- mean(rowSums(total))
  total_low <- quantile(rowSums(total), alpha / 2)
  total_high <- quantile(rowSums(total), 1 - alpha / 2)

  streams <- c(streams, "Total Fishery")
  efforts <- c(efforts, effort_total)
  takes <- c(takes, NA)

  # Build summary data frame
  # Determine indices for observed streams and pooled observed
  n_streams <- length(streams)
  n_obs_streams <- n_streams - 3  # Observer (+ optional EM, Both), plus Pooled, Unobserved, Total
  obs_idx <- seq_len(n_obs_streams)
  pooled_idx <- n_obs_streams + 1

  summary_df <- data.frame(
    stream = streams,
    effort = efforts,
    observed_takes = c(takes[obs_idx], takes[pooled_idx], NA, NA),
    estimated_mean = c(takes[obs_idx], takes_pooled, expanded_mean, total_mean),
    estimated_low = c(rep(NA, n_streams - 2), expanded_low, total_low),
    estimated_high = c(rep(NA, n_streams - 2), expanded_high, total_high),
    stringsAsFactors = FALSE
  )

  # Add coverage percentage
  summary_df$coverage_pct <- round(100 * summary_df$effort / effort_total, 1)

  return(summary_df)
}
