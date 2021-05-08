#' plot_fitted makes plots bycatch estimates (lambda of Poisson), accounting for effort but not accounting for observer coverage
#'
#' @param fitted_model Data and fitted model returned from fit_bycatch(). If a hurdle model, then only
#' then the plot returns the total bycatch rate (including zero and non-zero components).
#' @param xlab X-axis label for plot
#' @param ylab Y-axis label for plot
#' @param include_points whether or not to include raw bycatch events on plots, defaults to FALSE
#'
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
#' plot_fitted(fit,
#'   xlab = "Year", ylab = "Fleet-level bycatch",
#'   include_points = TRUE
#' )
#'
#' # fit a negative binomial model, with more chains and control arguments
#' fit_nb <- fit_bycatch(Takes ~ 1,
#'   data = d, time = "Year",
#'   effort = "Sets", family = "nbinom2",
#'   time_varying = FALSE, iter = 2000, chains = 4,
#'   control = list(adapt_delta = 0.99, max_treedepth = 20)
#' )
#'
#' # fit a time varying model
#' fit <- fit_bycatch(Takes ~ 1,
#'   data = d, time = "Year",
#'   effort = "Sets", family = "poisson", time_varying = TRUE
#' )
#'
#' # include data for expansion to unobserved sets
#' fit_nb <- fit_bycatch(Takes ~ 1,
#'   data = d, time = "Year",
#'   effort = "Sets", family = "nbinom2",
#'   expansion_rate = "expansionRate",
#'   time_varying = FALSE, iter = 2000, chains = 4,
#'   control = list(adapt_delta = 0.99, max_treedepth = 20)
#' )
#' }
plot_fitted <- function(fitted_model, xlab = "Time", ylab = "Events", include_points = FALSE) {
  lambda <- rstan::extract(fitted_model$fitted_model, c("lambda"))$lambda

  if (fitted_model$family %in% c("poisson-hurdle", "nbinom2-hurdle", "gamma-hurdle", "lognormal-hurdle", "normal-hurdle")) {
    # adjust lambda estimates by including theta. theta = pr(0), (1-theta) = pr(>0)
    theta <- rstan::extract(fitted_model$fitted_model, c("theta"))$theta
    for (i in 1:nrow(lambda)) lambda[i, ] <- lambda[i, ] * (1 - theta[i, 1])
  }
  df <- data.frame(
    "time" = fitted_model$data[, fitted_model$time],
    "mean" = apply(lambda, 2, mean),
    "low" = apply(lambda, 2, quantile, 0.025),
    "high" = apply(lambda, 2, quantile, 0.975),
    "obs" = fitted_model$data[, fitted_model$events]
  )

  # generate intervals for new data
  # n_sim = 10000
  # rand = matrix(rpois(n = length(df$mean)*n_sim, lambda=rep(df$mean, n_sim)), ncol = n_sim)
  # df$low_obs = apply(rand, 1, quantile, 0.025)
  # df$high_obs = apply(rand, 1, quantile, 0.975)

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
