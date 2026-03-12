# fit_bycatch - Unified version with coverage rate parameters

fit_bycatch - Unified version with coverage rate parameters

## Usage

``` r
fit_bycatch(
  formula,
  data,
  time = "year",
  effort = "effort",
  covrate = NULL,
  expansion_rate = NULL,
  takes_em = NULL,
  effort_em = NULL,
  covrate_em = NULL,
  takes_both = NULL,
  effort_both = NULL,
  covrate_obs = NULL,
  covrate_both = NULL,
  effort_total = NULL,
  family = c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle", "lognormal",
    "gamma", "lognormal-hurdle", "gamma-hurdle", "normal", "normal-hurdle"),
  time_varying = FALSE,
  iter = 1000,
  chains = 3,
  control = list(adapt_delta = 0.9, max_treedepth = 20),
  ...
)
```

## Arguments

- formula:

  The model formula.

- data:

  A data frame.

- time:

  Named column of the 'data' data frame with the label for the time
  (e.g. year) variable

- effort:

  Named column of the 'effort' variable in the data frame with the label
  for the fishing effort to be used in estimation of mean bycatch rate.
  This represents total observed effort

- covrate:

  Optional: Column name for coverage rate (percentage) for single-stream
  models. Replaces expansion_rate.

- expansion_rate:

  DEPRECATED: Use 'covrate' or 'covrate_obs' instead. The expansion rate
  to be used in generating distributions for unobserved sets.

- takes_em:

  Optional: Column name for bycatch takes recorded by EM (electronic
  monitoring). Default NULL.

- effort_em:

  Optional: Column name for effort monitored by EM. Required if takes_em
  is provided.

- covrate_em:

  Optional: Column name for EM coverage rate (percentage). Used to
  calculate unobserved effort.

- takes_both:

  Optional: Column name for bycatch takes when both Observer and EM are
  present. Default NULL.

- effort_both:

  Optional: Column name for effort with both Observer and EM monitoring.
  Required if takes_both is provided.

- covrate_obs:

  Optional: Column name for Observer coverage rate (percentage). Used to
  calculate unobserved effort.

- covrate_both:

  Optional: Column name for "Both" coverage rate (percentage). Used to
  calculate unobserved effort.

- effort_total:

  Optional: Column name for total fishery effort. If provided, takes
  precedence over coverage rates. Must be in same units as effort
  columns.

- family:

  Family for response distribution can be discrete ("poisson",
  "nbinom2", "poisson-hurdle","nbinom2-hurdle"), or continuous
  ("normal", "gamma","lognormal", "normal-hurdle", "gamma-hurdle",
  "lognormal-hurdle"). The default distribution is "poisson". The hurdle
  variants estimate the probability of zeros (theta) separately from the
  other models and use truncated distribution to model positive counts.
  All use a log link function.

- time_varying:

  boolean TRUE/FALSE, whether to include time varying component (this is
  a random walk, analogous to making this a Dynamic linear model)

- iter:

  the number of mcmc iterations, defaults to 1000

- chains:

  the number of mcmc chains, defaults to 3

- control:

  List to pass to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).
  For example, increase `adapt_delta` if there are warnings about
  divergent transitions: `control = list(adapt_delta = 0.99)`. By
  default, bycatch sets `adapt_delta = 0.9`.

- ...:

  Any other arguments to pass to
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Value

list of the data used to fit the model, the matrix of covariates, the
expanded bycatch generated via the fit and simulations, and the fitted
stan model

## Details

**Coverage Rates vs effort_total:**

This function supports two approaches for calculating unobserved effort:

1.  **Coverage rates** (recommended): Use `covrate`, `covrate_obs`,
    `covrate_em`, `covrate_both` to specify what percentage of the
    fishery is monitored. The function calculates unobserved effort as:
    `observed_effort * ((100 - coverage) / coverage)`

2.  **Total effort**: Use `effort_total` to provide the total fishery
    effort directly. Unobserved effort is calculated as:
    `effort_total - observed_effort` NOTE: effort_total must be in the
    same units as your effort columns

**Priority order**: effort_total \> coverage rates \> 100% coverage
assumed

**Multi-stream monitoring:** When using multiple monitoring streams
(Observer, EM, Both), the data are kept separate in the statistical
model but all share the same underlying bycatch rate. This approach
assumes all monitoring types have perfect/equal detection (detection =
100%) but avoids distributional issues that arise from pooling.

## Examples

``` r
# \donttest{
# Single stream example with coverage rate
d <- data.frame(
  "Year" = 2002:2014,
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "CovRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
)
fit <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year",
  effort = "Sets",
  covrate = "CovRate",
  family = "poisson",
  time_varying = FALSE
)
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.004 seconds (Warm-up)
#> Chain 1:                0.005 seconds (Sampling)
#> Chain 1:                0.009 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.004 seconds (Warm-up)
#> Chain 2:                0.004 seconds (Sampling)
#> Chain 2:                0.008 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.004 seconds (Warm-up)
#> Chain 3:                0.004 seconds (Sampling)
#> Chain 3:                0.008 seconds (Total)
#> Chain 3: 

# Multi-stream example with coverage rates (Observer + EM)
d_multi <- data.frame(
  Year = 2002:2014,
  Takes_obs = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0),
  Sets_obs = c(200, 180, 150, 350, 250, 270, 180, 150, 400, 350, 280, 180, 250),
  CovRate_obs = c(20, 19, 17, 19, 21, 21, 20, 20, 20, 19, 20, 19, 19),
  Takes_em = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
  Sets_em = c(150, 120, 140, 250, 180, 190, 120, 100, 300, 270, 200, 130, 190),
  CovRate_em = c(15, 13, 16, 14, 15, 15, 13, 13, 15, 15, 14, 14, 15)
)

fit_multi <- fit_bycatch(Takes_obs ~ 1,
  data = d_multi,
  time = "Year",
  effort = "Sets_obs",
  covrate_obs = "CovRate_obs",
  takes_em = "Takes_em",
  effort_em = "Sets_em",
  covrate_em = "CovRate_em",
  family = "poisson"
)
#> Observer stream: 13 observations
#> Total takes: 3
#> Total effort: 3190
#> EM stream: 13 observations
#> Total takes: 1
#> Total effort: 2340
#> Using coverage rates for expansion
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.005 seconds (Warm-up)
#> Chain 1:                0.005 seconds (Sampling)
#> Chain 1:                0.01 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.004 seconds (Warm-up)
#> Chain 2:                0.005 seconds (Sampling)
#> Chain 2:                0.009 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.005 seconds (Warm-up)
#> Chain 3:                0.005 seconds (Sampling)
#> Chain 3:                0.01 seconds (Total)
#> Chain 3: 
# }
```
