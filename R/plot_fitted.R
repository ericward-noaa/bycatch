#' plot_fitted makes plots bycatch estimates (lambda of Poisson), accounting for effort but not accounting for observer coverage
#'
#' @param fitted_model Data and fitted model returned from bycatch_expansion()
#' @param xlab X-axis label for plot
#' @param ylab Y-axis label for plot
#' @param include_points whether or not to include raw bycatch events on plots, defaults to FALSE
#'
#' @return plot called from ggplot
#'
#' @export
#' @import ggplot2
#' @import viridis
plot_fitted = function(fitted_model, xlab = "Time", ylab="Events", include_points=FALSE) {

  lambda = rstan::extract(fitted_model$fitted_model, c("lambda"))$lambda

  df = data.frame("time" = fitted_model$data$time,
    "mean" = apply(lambda, 2, mean),
    "low" = apply(lambda, 2, quantile, 0.025),
    "high" = apply(lambda, 2, quantile, 0.975),
    "obs" = fitted_model$data$events)

  # generate intervals for new data
  #n_sim = 10000
  #rand = matrix(rpois(n = length(df$mean)*n_sim, lambda=rep(df$mean, n_sim)), ncol = n_sim)
  #df$low_obs = apply(rand, 1, quantile, 0.025)
  #df$high_obs = apply(rand, 1, quantile, 0.975)

  g1 = ggplot(df, aes(time, mean)) +
    geom_ribbon(aes(ymin = low, ymax=high), fill="blue", alpha=0.3) +
    geom_line(color="blue") +
    xlab(xlab) + ylab(ylab) + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  if(include_points==TRUE) g1 = g1 + geom_point(aes(time, obs), size=2, color="blue")
  return(g1)
}
