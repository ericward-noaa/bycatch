#' plot_expanded is makes plots of the expanded bycatch estimates, accounting for observer coverage and effort
#'
#' @param fitted_model Data and fitted model returned from bycatch_expansion()
#' @param xlab X-axis label for plot
#' @param ylab Y-axis label for plot
#'
#' @return plot called from ggplot
#'
#' @export
plot_expanded = function(fitted_model, xlab = "Time", ylab="Events") {

  df = data.frame("time" = fitted_model$data$time,
    "mean" = apply(fitted_model$expanded_estimates, 2, mean),
    "low" = apply(fitted_model$expanded_estimates, 2, quantile, 0.025),
    "high" = apply(fitted_model$expanded_estimates, 2, quantile, 0.975),
    "obs" = fitted_model$data$events)

  g1 = ggplot(df, aes(time, mean)) +
    geom_ribbon(aes(ymin = low, ymax=high), fill="grey60") +
    geom_line() +
    geom_point(aes(time, obs), size=2) +
    xlab(xlab) + ylab(ylab) + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  return(g1)
}
