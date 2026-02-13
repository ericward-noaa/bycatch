# get_fitted returns df of observed bycatch estimates (lambda of Poisson), accounting for effort but not accounting for observer coverage

get_fitted returns df of observed bycatch estimates (lambda of Poisson),
accounting for effort but not accounting for observer coverage

## Usage

``` r
get_fitted(fitted_model, alpha = 0.05)
```

## Arguments

- fitted_model:

  Data and fitted model returned from fit_bycatch(). If a hurdle model,
  then the plot returns the total bycatch rate (including zero and
  non-zero components).

- alpha:

  The alpha level for the credible interval, defaults to 0.05

## Value

plot called from ggplot

## Examples

``` r
# \donttest{
d <- data.frame(
  "Year" = 2002:2014,
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "expansionRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
)
fit <- fit_bycatch(Takes ~ 1,
  data = d, time = "Year", effort = "Sets",
  family = "poisson", time_varying = FALSE
)
#> No expansion information provided - assuming 100% coverage
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 9e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
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
#> Chain 1:  Elapsed Time: 0.016 seconds (Warm-up)
#> Chain 1:                0.015 seconds (Sampling)
#> Chain 1:                0.031 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 7e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.07 seconds.
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
#> Chain 2:  Elapsed Time: 0.014 seconds (Warm-up)
#> Chain 2:                0.014 seconds (Sampling)
#> Chain 2:                0.028 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 7e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.07 seconds.
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
#> Chain 3:  Elapsed Time: 0.012 seconds (Warm-up)
#> Chain 3:                0.013 seconds (Sampling)
#> Chain 3:                0.025 seconds (Total)
#> Chain 3: 
get_fitted(fit)
#>    time      mean        low      high obs
#> 1  2002 0.2811414 0.08014784 0.6027706   0
#> 2  2003 0.2444708 0.06969377 0.5241483   0
#> 3  2004 0.2372804 0.06764395 0.5087322   0
#> 4  2005 0.4745609 0.13528791 1.0174644   0
#> 5  2006 0.3379449 0.09634139 0.7245580   0
#> 6  2007 0.3595158 0.10249084 0.7708064   0
#> 7  2008 0.2372804 0.06764395 0.5087322   0
#> 8  2009 0.2063621 0.05882974 0.4424429   0
#> 9  2010 0.5435879 0.15496615 1.1654592   1
#> 10 2011 0.4839083 0.13795267 1.0375054   3
#> 11 2012 0.3825248 0.10905025 0.8201380   0
#> 12 2013 0.2523801 0.07194857 0.5411061   0
#> 13 2014 0.3494494 0.09962110 0.7492238   0
# }
```
