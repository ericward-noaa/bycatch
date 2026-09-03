#' plot_pred_distribution plots the distribution of bycatch prediction for a specific year
#'
#' @param fitted_model Data and fitted model returned from fit_bycatch().
#' @param year Specified year for the prediction distribution
#' @param xlab X-axis label for plot
#' @param ylab Y-axis label for plot
#' @param breaks Breaks for bars in the plot, defaults to 1
#'
#' @return A ggplot object.
#'
#' @export
#' @import ggplot2
#' @importFrom rlang .data
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
#' plot_pred_distribution(fit, year = 2002)
#' }
#'
plot_pred_distribution <- function(fitted_model, year, xlab = "Predicted Bycatch", ylab = "Frequency", breaks = NULL) {
  exp <- bycatch::get_total(fitted_model)

  model_years <- fitted_model$data$Year
  year_index <- which(model_years == year)

  if (!(year %in% model_years)) {
    stop(paste("Year", year, "not found in data"))
  }

  predicted <- as.numeric(exp[, year_index])

  if (is.null(breaks)) {
    max_value <- ceiling(max(predicted))
    if (max_value < 1) {
      breaks <- c(0, 1)
    } else {
      breaks <- seq(0, max_value, by = 1)
    }
  }

  # Create a data frame for ggplot
  df <- data.frame(predicted = predicted)

  # Plot histogram using ggplot2
  ggplot(df, aes(x = predicted)) +
    geom_histogram(breaks = breaks, fill = "grey80", color = "black") +
    stat_bin(
      breaks = breaks,
      geom = "text",
      aes(label = after_stat(.data$count)),
      vjust = -0.5
    ) +
    scale_x_continuous(limits = c(0, max(predicted)), breaks = breaks) +
    labs(
      title = paste("Distribution of Predicted Bycatch -", year),
      x = xlab,
      y = ylab
    ) +
    theme_bw()
}
