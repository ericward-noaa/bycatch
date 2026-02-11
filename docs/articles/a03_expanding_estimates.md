# Expanding bycatch estimates

``` r
#library(devtools)
#devtools::install_github("ericward-noaa/bycatch")
library(bycatch)
set.seed(123)
```

### Load data

For our dummy dataaset, the expansionRate column contains the expansion
rate. In other words, the number of total sets is the sets divided by
expansion rate.

``` r
# replace this with your own data frame
d = data.frame("Year"= 2002:2014, 
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "expansionRate" = c(24, 22, 14, 32, 28, 25, 30,  7, 26, 21, 22, 23, 27),
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
But to do the expansion, we can include the `expansion_rate` argument,

``` r
fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = FALSE, expansion_rate = "expansionRate")
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
