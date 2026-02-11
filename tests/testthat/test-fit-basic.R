# Test file for bycatch model fitting

library(testthat)
library(bycatch)

# Setup test data
d_pois <- data.frame(
  "Year" = 2002:2014,
  "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
  "CovRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
)

d_nb <- data.frame(
  "Year" = 2002:2014,
  "Takes" = c(0, 0, 0, 3, 0, 2, 0, 5, 1, 3, 0, 0, 2),
  "CovRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
)

d_pos <- data.frame(
  "Year" = 2002:2014,
  "Takes" = c(0.5, 0.2, 0.8, 1.3, 0.4, 1.2, 0.7, 1.5, 1.1, 2.3, 0.6, 0.9, 1.4),
  "CovRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
  "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
)

# Multi-stream test data
d_multi <- data.frame(
  Year = 2002:2014,
  Takes_obs = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0),
  Sets_obs = c(200, 180, 150, 350, 250, 270, 180, 150, 400, 350, 280, 180, 250),
  CovRate_obs = c(20, 19, 17, 19, 21, 21, 20, 20, 20, 19, 20, 19, 19),
  Takes_em = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
  Sets_em = c(150, 120, 140, 250, 180, 190, 120, 100, 300, 270, 200, 130, 190),
  CovRate_em = c(15, 13, 16, 14, 15, 15, 13, 13, 15, 15, 14, 14, 15)
)

test_that("fitting function works for poisson model with covrate", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pois, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "poisson",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "poisson")
  expect_false(fit$multi_stream)
})

test_that("fitting function works for negative binomial2 model with covrate", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_nb, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "nbinom2",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "nbinom2")
})

test_that("fitting function works for hurdle poisson model with covrate", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pois, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "poisson-hurdle",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "poisson-hurdle")
})

test_that("backward compatibility: expansion_rate still works (with warning)", {
  expect_warning(
    fit <- fit_bycatch(Takes ~ 1,
                       data = d_pois, time = "Year",
                       effort = "Sets",
                       expansion_rate = "CovRate",  # Old parameter name
                       family = "poisson",
                       time_varying = FALSE
    ),
    "expansion_rate.*deprecated"
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "poisson")
})

test_that("fitting function works for time varying poisson model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pois, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "poisson",
                     time_varying = TRUE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "poisson")

  # Check that time-varying parameters exist
  params <- names(rstan::extract(fit$fitted_model))
  expect_true("sigma_rw" %in% params)
  expect_true("est_time_dev" %in% params)
})

test_that("fitting function works without expansion (100% coverage)", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pois, time = "Year",
                     effort = "Sets",
                     # No covrate or expansion_rate - should assume 100% coverage
                     family = "poisson",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")

  # Check that new_effort is all zeros (no expansion)
  expanded <- get_expanded(fit)
  # With 100% coverage, expanded estimates should be zero or very small
})

test_that("fitting function works for time varying poisson hurdle model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pois, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "poisson-hurdle",
                     time_varying = TRUE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "poisson-hurdle")
})

test_that("fitting function works for gamma model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pos, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "gamma",
                     time_varying = FALSE,
                     chains = 1, iter = 1000
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "gamma")
})

test_that("fitting function works for normal model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pos, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "normal",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "normal")
})

test_that("fitting function works for lognormal model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pos, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "lognormal",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "lognormal")
})

test_that("fitting function works for normal-hurdle model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pos, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "normal-hurdle",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "normal-hurdle")
})

test_that("fitting function works for lognormal-hurdle model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pos, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "lognormal-hurdle",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "lognormal-hurdle")
})

test_that("fitting function works for gamma-hurdle model", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pos, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "gamma-hurdle",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_equal(fit$family, "gamma-hurdle")
})

test_that("multi-stream model works with coverage rates", {
  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi,
                     time = "Year",
                     effort = "Sets_obs",
                     covrate_obs = "CovRate_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     covrate_em = "CovRate_em",
                     family = "poisson",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_true(fit$multi_stream)
  expect_equal(fit$stream_info$n_obs, 13)
  expect_equal(fit$stream_info$n_em, 13)
})

test_that("multi-stream model works without coverage rates (100% coverage)", {
  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi,
                     time = "Year",
                     effort = "Sets_obs",
                     # No coverage rates provided
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     family = "poisson",
                     time_varying = FALSE
  )

  expect_equal(class(fit$fitted_model)[[1]], "stanfit")
  expect_true(fit$multi_stream)
})

test_that("error handling: missing required columns", {
  expect_error(
    fit_bycatch(Takes ~ 1,
                data = d_pois, time = "BadColumn",
                effort = "Sets",
                family = "poisson"
    ),
    "time variable"
  )

  expect_error(
    fit_bycatch(Takes ~ 1,
                data = d_pois, time = "Year",
                effort = "BadColumn",
                family = "poisson"
    ),
    "effort variable"
  )
})

test_that("error handling: invalid family", {
  expect_error(
    fit_bycatch(Takes ~ 1,
                data = d_pois, time = "Year",
                effort = "Sets",
                family = "not_a_family"
    ),
    "family must be specified"
  )
})

test_that("error handling: multi-stream missing effort", {
  expect_error(
    fit_bycatch(Takes_obs ~ 1,
                data = d_multi,
                time = "Year",
                effort = "Sets_obs",
                takes_em = "Takes_em",
                # Missing effort_em
                family = "poisson"
    ),
    "effort_em must also be provided"
  )
})

test_that("covrate and covrate_obs are equivalent for single-stream", {
  # Using covrate
  fit1 <- fit_bycatch(Takes ~ 1,
                      data = d_pois, time = "Year",
                      effort = "Sets",
                      covrate = "CovRate",
                      family = "poisson",
                      time_varying = FALSE
  )

  # Using covrate_obs
  fit2 <- fit_bycatch(Takes ~ 1,
                      data = d_pois, time = "Year",
                      effort = "Sets",
                      covrate_obs = "CovRate",
                      family = "poisson",
                      time_varying = FALSE
  )

  # Both should work and give similar results
  expect_equal(class(fit1$fitted_model)[[1]], "stanfit")
  expect_equal(class(fit2$fitted_model)[[1]], "stanfit")
})

test_that("extract_log_lik_for_loo works", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pois, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "poisson",
                     time_varying = FALSE
  )

  log_lik <- extract_log_lik_for_loo(fit)

  expect_true(is.matrix(log_lik))
  expect_equal(ncol(log_lik), nrow(d_pois))  # One column per observation
  expect_true(all(!is.na(log_lik)))  # No NAs after filtering
})

test_that("get_expanded works", {
  fit <- fit_bycatch(Takes ~ 1,
                     data = d_pois, time = "Year",
                     effort = "Sets",
                     covrate = "CovRate",
                     family = "poisson",
                     time_varying = FALSE
  )

  expanded <- get_expanded(fit)

  expect_true(is.matrix(expanded))
  # Should have one column per year
  expect_equal(ncol(expanded), length(unique(d_pois$Year)))
})
