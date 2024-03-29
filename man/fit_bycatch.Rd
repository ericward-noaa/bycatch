% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{fit_bycatch}
\alias{fit_bycatch}
\title{fit_bycatch is the primary function for fitting bycatch models to time series of takes and effort}
\usage{
fit_bycatch(
  formula,
  data,
  time = "year",
  effort = "effort",
  expansion_rate = NULL,
  family = c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle", "lognormal",
    "gamma", "lognormal-hurdle", "gamma-hurdle", "normal", "normal-hurdle"),
  time_varying = FALSE,
  iter = 1000,
  chains = 3,
  control = list(adapt_delta = 0.9, max_treedepth = 20),
  ...
)
}
\arguments{
\item{formula}{The model formula.}

\item{data}{A data frame.}

\item{time}{Named column of the 'data' data frame with the label for the time (e.g. year) variable}

\item{effort}{Named column of the 'effort' variable in the data frame with the label for the fishing effort to be used in estimation of mean bycatch rate. This
represents total observed effort}

\item{expansion_rate}{The expansion rate to be used in generating distributions for unobserved sets. If NULL, defaults to 100\% coverage (= 100)}

\item{family}{Family for response distribution can be discrete ("poisson",
"nbinom2", "poisson-hurdle","nbinom2-hurdle"), or continuous ("normal",
"gamma","lognormal", "normal-hurdle", "gamma-hurdle", "lognormal-hurdle"). The
default distribution is "poisson". The hurdle variants estimate the
probability of zeros (theta) separately from the other models and use
truncated distribution to model positive counts. All use a log
link function.}

\item{time_varying}{boolean TRUE/FALSE, whether to include time varying component (this is a random walk, analogous to making this a Dynamic linear model)}

\item{iter}{the number of mcmc iterations, defaults to 1000}

\item{chains}{the number of mcmc chains, defaults to 3}

\item{control}{List to pass to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}. For example,
increase \code{adapt_delta} if there are warnings about divergent
transitions: \code{control = list(adapt_delta = 0.99)}. By default,
\pkg{bycatch} sets \code{adapt_delta = 0.9}.}

\item{...}{Any other arguments to pass to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}.}
}
\value{
list of the data used to fit the model, the matrix of covariates, the expanded bycatch generated via the fit and simulations, and the fitted stan model
}
\description{
fit_bycatch is the primary function for fitting bycatch models to time series of takes and effort
}
\examples{
\donttest{
d <- data.frame(
  "Year" = 2002:2014,
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "expansionRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
)
fit <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year",
  effort = "Sets", family = "poisson", time_varying = FALSE
)
loo::loo(fit$fitted_model)$estimates

fit <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year", effort = "Sets",
  family = "poisson", time_varying = FALSE
)
plot_fitted(fit,
  xlab = "Year", ylab = "Fleet-level bycatch",
  include_points = TRUE
)

# fit a negative binomial model, with more chains and control arguments
fit_nb <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year",
  effort = "Sets", family = "nbinom2",
  time_varying = FALSE, iter = 2000, chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

# fit a time varying model
fit <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year",
  effort = "Sets", family = "poisson", time_varying = TRUE
)

# include data for expansion to unobserved sets
fit_nb <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year",
  effort = "Sets", family = "nbinom2",
  expansion_rate = "expansionRate",
  time_varying = FALSE, iter = 2000, chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

# fit a model with a lognormal distribution
d$Takes <- rnorm(nrow(d), 5, 0.1)
fit_ln <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year",
  effort = "Sets", family = "lognormal",
  expansion_rate = "expansionRate",
  time_varying = FALSE, iter = 2000, chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

# add zeros and fit a delta-gamma distribution
d$Takes <- rnorm(nrow(d), 5, 0.1)
d$Takes[c(1, 5, 10)] <- 0
fit_ln <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year",
  effort = "Sets", family = "gamma-hurdle",
  expansion_rate = "expansionRate",
  time_varying = FALSE, iter = 2000, chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)
}
}
