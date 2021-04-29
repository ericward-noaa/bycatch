#' plot_expanded is makes plots of the expanded bycatch estimates, accounting for observer coverage and effort
#'
#' @param fitted_model Data and fitted model returned from estimation
#' @param expanded_estimates values returned from a call to expansion()
#' @param xlab X-axis label for plot
#' @param ylab Y-axis label for plot
#' @param include_points whether or not to include raw bycatch events on plots, defaults to FALSE
#'
#' @return plot called from ggplot
#'
#' @export
#' @import ggplot2
#' @importFrom stats quantile
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
#'   time_varying = FALSE
#' )
#' expanded <- expand(fit, coverage = d$expansionRate)
#' plot_expanded(
#'   fitted_model = fit,
#'   expanded_estimates = expanded,
#'   xlab = "Year",
#'   ylab = "Fleet-level bycatch",
#'   include_points = TRUE
#' )
#' }
plot_expanded <- function(fitted_model, expanded_estimates, xlab = "Time", ylab = "Events", include_points = FALSE) {
  df <- data.frame(
    "time" = fitted_model$data[, fitted_model$time],
    "mean" = apply(expanded_estimates, 2, mean),
    "low" = apply(expanded_estimates, 2, quantile, 0.025),
    "high" = apply(expanded_estimates, 2, quantile, 0.975),
    "obs" = fitted_model$data[, fitted_model$events]
  )

  g1 <- ggplot(df, aes(.data$time, mean)) +
    geom_ribbon(aes(ymin = .data$low, ymax = .data$high), fill = "blue", alpha = 0.3) +
    geom_line(color = "blue") +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw() +
    theme(
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
    )
  if (include_points == TRUE) g1 <- g1 + geom_point(aes(.data$time, .data$obs), size = 2, color = "blue")
  return(g1)
}
