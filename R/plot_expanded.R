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
plot_expanded = function(fitted_model, expanded_estimates, xlab = "Time", ylab="Events", include_points=FALSE) {

  df = data.frame("time" = fitted_model$data[,fitted_model$time],
    "mean" = apply(expanded_estimates, 2, mean),
    "low" = apply(expanded_estimates, 2, quantile, 0.025),
    "high" = apply(expanded_estimates, 2, quantile, 0.975),
    "obs" = fitted_model$data[,fitted_model$events])

  g1 = ggplot(df, aes(time, mean)) +
    geom_ribbon(aes(ymin = low, ymax=high), fill="blue", alpha=0.3) +
    geom_line(color = "blue") +
    xlab(xlab) + ylab(ylab) + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  if(include_points==TRUE) g1 = g1 + geom_point(aes(time, obs), size=2, color = "blue")
  return(g1)
}
