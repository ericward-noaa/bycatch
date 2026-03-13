# Diagnosing errors

``` r
library(bycatch)
library(loo)
set.seed(123)
```

In some cases, models might produce warnings about high Pareto K
statistics (\> 0.7). These may reflect a model that is too flexible for
the data being used, in which case the model might need to be changed.
In some cases, ‘moment matching’ may be used to improve diagnostics; an
example of this is below.

### Load data

``` r
# replace this with your own data frame
d = data.frame("Year"= 2002:2014, 
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486))
```

## Simple model with constant bycatch, no covariates

This model has a constant bycatch rate and Poisson distribution,

``` r
fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = FALSE)
```

If the ‘fit’ object produced errors about the Pareto K statistics being
high, we could use moment matching in the `loo` package,

``` r
loo_stats <- loo(fit$fitted_model, moment_match = TRUE)

print(loo_stats)
```

    ## 
    ## Computed from 1500 by 13 log-likelihood matrix.
    ## 
    ##          Estimate  SE
    ## elpd_loo      0.0 0.0
    ## p_loo         0.0 0.0
    ## looic         0.0 0.0
    ## ------
    ## MCSE of elpd_loo is NA.
    ## MCSE and ESS estimates assume independent draws (r_eff=1).
    ## 
    ## All Pareto k estimates are good (k < 0.69).
    ## See help('pareto-k-diagnostic') for details.
