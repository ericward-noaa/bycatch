# Expanding bycatch estimates

``` r
#library(devtools)
#devtools::install_github("ericward-noaa/bycatch")
library(bycatch)
set.seed(123)
```

### Load data

For our dummy dataaset, the covRate column contains the expansion rate.
In other words, the number of total sets is the sets divided by
expansion rate.

``` r
# replace this with your own data frame
d = data.frame("Year"= 2002:2014, 
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "covRate" = c(24, 22, 14, 32, 28, 25, 30,  7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486))
```

## Simple model with constant bycatch, no covariates

We’ll start by fitting a model with constant bycatch rate,

``` r
fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = FALSE)
```

## Expanding bycatch estimates

By default the model assumes coverage is 100%, and no expansion is done.
But to do the expansion, we can include the `covrate` argument,

``` r
d = data.frame("Year"= 2002:2014, 
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "covRate" = c(24, 22, 14, 32, 28, 25, 30,  7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486))

fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = FALSE, covrate = "covRate")
pars <- rstan::extract(fit$fitted_model)
# beta is the underlying bycatch rate before scaling for effort
rpois(1000, lambda = exp(pars$beta + log(500)))
```

And we can then plot these estimates. Like the previous function we can
specify whether to include the raw points or not. First, we’ll show
estimates of total bycatch (observed + unobserved)

``` r
plot_expanded(fitted_model=fit, xlab="Year", ylab = "Fleet-level bycatch", include_points = TRUE)
```

Alternatively, we may be interested in just summarizing the unobserved
estimates. In this case, we can set the parameter `show_total` = FALSE.

``` r
plot_expanded(fitted_model=fit, xlab="Year", ylab = "Unobserved bycatch", include_points = TRUE,show_total = FALSE)
```

### Make table of expanded bycatch estimates

We can also do things like summarize the expanded estimates in table
form. This takes use of the `get_expanded` helper function,

``` r
expanded = get_expanded(fit)

df = data.frame("time" = d[,"Year"], 
  "mean" = apply(expanded, 2, mean),
  "median" = apply(expanded, 2, quantile, 0.5),
  "lower95" = apply(expanded, 2, quantile, 0.025),
  "upper95" = apply(expanded, 2, quantile, 0.975))
```

We could do all kinds of things with this kind of table, like output it
to a .csv,

``` r
write.table(df, "estimated_bycatch.csv", row.names=F, col.names=T, sep=",")
```

### Make table of total bycatch estimates

Like the `get_expanded` helper function, we have a `get_total` helper
function to summarize the total observed + expanded unobserved bycatch.
This function also takes a fitted model object,

``` r
total = get_total(fit)

df = data.frame("time" = d[,"Year"], 
  "mean" = apply(total, 2, mean),
  "median" = apply(total, 2, quantile, 0.5),
  "lower95" = apply(total, 2, quantile, 0.025),
  "upper95" = apply(total, 2, quantile, 0.975))
```

### Derived quantities

Or calculate the probability of getting values \> 0, by year

``` r
pr_zero = function(x) {
  return(length(which(x==0))/length(x))
}

pr_pos = 1 - apply(expanded,2, pr_zero)
print(pr_pos)
```

    ##  [1] 0.5566667 0.5493333 0.7080000 0.6140000 0.5466667 0.6266667 0.4093333
    ##  [8] 0.8933333 0.7233333 0.7800000 0.6820000 0.5193333 0.5866667

Because values in `total` represent a probability distribution, we could
alternatively ask what’s the probability of getting a value above some
threshold. For this example, we’ll arbitrarily use a threshold of 4, but
you can calculate

``` r
apply(total, 2, function(x) sum(x > 4) / length(x))
```

    ##  [1] 0.0040000000 0.0040000000 0.0346666667 0.0106666667 0.0060000000
    ##  [6] 0.0133333333 0.0006666667 0.1693333333 0.0860000000 0.5033333333
    ## [11] 0.0280000000 0.0026666667 0.0080000000

``` r
# or equivalently
#colMeans(total > 4)
```

### Simulating future values

One potential use of this output is to also generate distributions of
future bycatch. The workflow to do this is (1) fit the model, (2)
extract the parameters (focusing on `beta`, the bycatch rate unadjusted
for effort), and (3) simulating the distribution of future bycatch for
known effort.

``` r
d = data.frame("Year"= 2002:2014,
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "covRate" = c(24, 22, 14, 32, 28, 25, 30,  7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486))

# fit the model
fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = FALSE, covrate = "covRate")
```

    ## 
    ## SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 5e-06 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 0.004 seconds (Warm-up)
    ## Chain 1:                0.003 seconds (Sampling)
    ## Chain 1:                0.007 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 1e-06 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 0.004 seconds (Warm-up)
    ## Chain 2:                0.004 seconds (Sampling)
    ## Chain 2:                0.008 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL 'bycatch' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 1e-06 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
    ## Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
    ## Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
    ## Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
    ## Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
    ## Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
    ## Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
    ## Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
    ## Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
    ## Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
    ## Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
    ## Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 0.004 seconds (Warm-up)
    ## Chain 3:                0.004 seconds (Sampling)
    ## Chain 3:                0.008 seconds (Total)
    ## Chain 3:

``` r
# extract parameters
pars <- rstan::extract(fit$fitted_model)

# beta is the underlying bycatch rate before scaling for effort - 500 is the future sets
future_bycatch <- rpois(1000, lambda = exp(pars$beta + log(500)))
```
