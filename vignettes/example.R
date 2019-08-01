## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- results="hide", message=FALSE, warning=FALSE-----------------------
library(devtools)
devtools::install_github("ericward-noaa/bycatch")
library(bycatch)
set.seed(123)


## ----data----------------------------------------------------------------
# replace this with your own data frame
d = data.frame("Year"= 2002:2014, 
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "expansionRate" = c(24, 22, 14, 32, 28, 25, 30,  7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486))


## ---- results="hide", message=FALSE, warning=FALSE-----------------------
fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = FALSE)


## ----eval=FALSE----------------------------------------------------------
## fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
##   time_varying = FALSE, control=list(adapt_delta=0.99,max_treedepth=20))


## ----eval=FALSE----------------------------------------------------------
## fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
##   time_varying = FALSE, iter=3000, chains=4)


## ---- fig.pos="placeHere", fig.cap = "Estimated bycatch for observed vessels (not expanded by observer coverage), but including observed takes and effort."----
plot_fitted(fit, xlab="Year", ylab = "Bycatch (observed vessels)")


## ---- fig.pos="placeHere", fig.cap = "Observed bycatch (not expanded by observer coverage), incorporating data on observed takes and effort. Dots represent observed bycatch events."----
plot_fitted(fit, xlab="Year", ylab = "Estimated bycatch for observed vessels", include_points = TRUE)


## ------------------------------------------------------------------------
loo::loo(fit$fitted_model)$estimates


## ------------------------------------------------------------------------
expanded = expand(fit, coverage = d$expansionRate)


## ---- fig.pos="placeHere", fig.cap = "Estimated fleet-level expanded bycatch, incorporating data on takes, effort, and observer coverage. Dots represent observed bycatch events."----
plot_expanded(fitted_model=fit, expanded_estimates = expanded, xlab="Year", ylab = "Fleet-level bycatch", include_points = TRUE)


## ---- eval=FALSE---------------------------------------------------------
## df = data.frame("time" = d[,"Year"],
##   "mean" = apply(expanded, 2, mean),
##   "median" = apply(expanded, 2, quantile, 0.5),
##   "lower95" = apply(expanded, 2, quantile, 0.025),
##   "upper95" = apply(expanded, 2, quantile, 0.975))
## 
## write.table(df, "estimated_bycatch.csv", row.names=F, col.names=T, sep=",")


## ---- results="hide", message=FALSE, warning=FALSE-----------------------
fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="nbinom2",
  time_varying = FALSE)


## ------------------------------------------------------------------------
phi = rstan::extract(fit$fitted_model)$nb2_phi


## ------------------------------------------------------------------------
d$Day = sample(seq(220,280),size=nrow(d),replace=T)
d$Reg = ifelse(d$Year < 2008, 0, 1)


## ----message=FALSE, warning=FALSE,results="hide"-------------------------
fit = fit_bycatch(Takes ~ Day + Reg, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = FALSE)


## ------------------------------------------------------------------------
betas = rstan::extract(fit$fitted_model)$beta


## ---- results="hide", message=FALSE, warning=FALSE-----------------------
fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
  time_varying = TRUE)


## ---- fig.pos="placeHere", fig.cap = "Estimated fleet-level expanded bycatch from the model with time-varying effects, incorporating data on takes, effort, and observer coverage. Dots represent observed bycatch events."----
expanded = expand(fit, coverage=d$expansionRate)
plot_expanded(fitted_model=fit, expanded_estimates = expanded, xlab="Year", ylab = "Fleet-level bycatch", include_points = TRUE)


## ----data2---------------------------------------------------------------
# replace this with your own data frame
d = data.frame("Year"= 2002:2014, 
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  "expansionRate" = c(24, 22, 14, 32, 28, 25, 30,  7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486))


## ---- eval = FALSE-------------------------------------------------------
## fit = fit_bycatch(Takes ~ 1, data=d, time="Year", effort="Sets", family="poisson",
##   time_varying = FALSE)


## ----fig.pos="placeHere", fig.cap = "Estimated fleet-level expanded bycatch from the dataset with no events, incorporating effort, and observer coverage. Dots represent observed bycatch events."----
expanded = expand(fit, coverage=d$expansionRate)
plot_expanded(fitted_model=fit, expanded_estimates = expanded, xlab="Year", ylab = "Fleet-level bycatch", include_points = TRUE)

