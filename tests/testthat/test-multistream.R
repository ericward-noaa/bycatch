context("multi-stream monitoring tests")

# Set up test data
sample_size <- 20

# Two-stream data (Observer + EM)
set.seed(456)
d_multi_pois <- data.frame(
  Year = 1:sample_size,
  Takes_obs = rpois(sample_size, 1.0),
  Sets_obs = rep(100, sample_size),
  Takes_em = rpois(sample_size, 0.8),
  Sets_em = rep(80, sample_size),
  Sets_total = rep(500, sample_size)
)

# Three-stream data (Observer + EM + Both)
set.seed(456)
d_multi_three <- data.frame(
  Year = 1:sample_size,
  Takes_obs = rpois(sample_size, 1.0),
  Sets_obs = rep(100, sample_size),
  Takes_em = rpois(sample_size, 0.8),
  Sets_em = rep(80, sample_size),
  Takes_both = rpois(sample_size, 0.5),
  Sets_both = rep(50, sample_size),
  Sets_total = rep(500, sample_size)
)

# Continuous data for Gamma/Lognormal
set.seed(456)
d_multi_cont <- data.frame(
  Year = 1:sample_size,
  Takes_obs = abs(rnorm(sample_size, 5, 0.5)),
  Sets_obs = rep(100, sample_size),
  Takes_em = abs(rnorm(sample_size, 4.8, 0.5)),
  Sets_em = rep(80, sample_size),
  Sets_total = rep(500, sample_size)
)

# ===== BASIC FUNCTIONALITY TESTS =====

test_that("multi-stream mode activates correctly for two streams", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
  expect_equal(fit$stream_info$n_obs, sample_size)
  expect_equal(fit$stream_info$n_em, sample_size)
  expect_equal(fit$stream_info$n_both, 0)
})

test_that("multi-stream mode activates correctly for three streams", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_three,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     takes_both = "Takes_both",
                     effort_both = "Sets_both",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
  expect_equal(fit$stream_info$n_obs, sample_size)
  expect_equal(fit$stream_info$n_em, sample_size)
  expect_equal(fit$stream_info$n_both, sample_size)
})

test_that("single-stream mode still works (backwards compatible)", {
  set.seed(123)

  d_single <- data.frame(
    Year = 1:sample_size,
    Takes = rpois(sample_size, 1.0),
    Sets = rep(100, sample_size)
  )

  fit <- fit_bycatch(Takes ~ 1,
                     data = d_single,
                     time = "Year",
                     effort = "Sets",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_false(fit$multi_stream)
  expect_null(fit$stream_info)
})

# ===== DISTRIBUTIONAL FAMILY TESTS =====

test_that("multi-stream works with Poisson", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)

  # Check lambda_base was estimated
  pars <- rstan::extract(fit$fitted_model, "lambda_base")
  expect_type(pars$lambda_base, "double")
})

test_that("multi-stream works with Negative Binomial", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "nbinom2",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
})

test_that("multi-stream works with Gamma (this is the key test!)", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_cont,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "gamma",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)

  # This would have FAILED before Phase 1.5!
  pars <- rstan::extract(fit$fitted_model, "lambda_base")
  expect_type(pars$lambda_base, "double")
})

test_that("multi-stream works with Lognormal", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_cont,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "lognormal",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
})

test_that("multi-stream works with Normal", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_cont,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "normal",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
})

test_that("multi-stream works with Poisson-Hurdle", {
  set.seed(123)

  # Add some zeros
  d_hurdle <- d_multi_pois
  d_hurdle$Takes_obs[c(1, 5, 10)] <- 0
  d_hurdle$Takes_em[c(2, 7, 12)] <- 0

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_hurdle,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson-hurdle",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
})

test_that("multi-stream works with Gamma-Hurdle", {
  set.seed(123)

  # Add some zeros
  d_hurdle <- d_multi_cont
  d_hurdle$Takes_obs[c(1, 5, 10)] <- 0
  d_hurdle$Takes_em[c(2, 7, 12)] <- 0

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_hurdle,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "gamma-hurdle",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
})

# ===== TIME-VARYING TESTS =====

test_that("multi-stream works with time-varying effects", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = TRUE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)

  # Check time deviations were estimated
  pars <- rstan::extract(fit$fitted_model, "est_time_dev")
  expect_type(pars$est_time_dev, "double")
})

# ===== EXTRACTION FUNCTION TESTS =====

test_that("get_expanded works with multi-stream", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expanded <- get_expanded(fit)

  expect_equal(ncol(expanded), sample_size)
  expect_true(nrow(expanded) > 0)
})

test_that("get_fitted works with multi-stream", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  fitted <- get_fitted(fit)

  expect_s3_class(fitted, "data.frame")
  expect_true("time" %in% names(fitted))
  expect_true("mean" %in% names(fitted))
  expect_true("low" %in% names(fitted))
  expect_true("high" %in% names(fitted))
  expect_true("obs" %in% names(fitted))
  expect_equal(nrow(fitted), sample_size)
})

test_that("get_total works with multi-stream", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  total <- get_total(fit)
  expanded <- get_expanded(fit)

  expect_equal(dim(total), dim(expanded))

  # Total should be greater than expanded (since it includes observed)
  expect_true(mean(rowSums(total)) > mean(rowSums(expanded)))
})

test_that("get_stream_summary works with multi-stream", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  summary <- get_stream_summary(fit)

  expect_s3_class(summary, "data.frame")
  expect_true("stream" %in% names(summary))
  expect_true("effort" %in% names(summary))
  expect_true("coverage_pct" %in% names(summary))

  # Should have at least: Observer, EM, Pooled, Unobserved, Total
  expect_gte(nrow(summary), 5)

  # Check that efforts sum correctly
  total_row <- summary[summary$stream == "Total Fishery", ]
  expect_equal(total_row$effort, sum(d_multi_pois$Sets_total))
})

test_that("get_stream_summary works with three streams", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_three,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     takes_both = "Takes_both",
                     effort_both = "Sets_both",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  summary <- get_stream_summary(fit)

  # Should include "Both" stream
  expect_true("Both" %in% summary$stream)

  # Should have: Observer, EM, Both, Pooled, Unobserved, Total
  expect_gte(nrow(summary), 6)
})

# ===== PLOTTING FUNCTION TESTS =====

test_that("plot_fitted works with multi-stream", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  p <- plot_fitted(fit)

  expect_s3_class(p, "ggplot")
  expect_false(is.null(p))
})

test_that("plot_expanded works with multi-stream", {
  set.seed(123)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_multi_pois,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  p <- plot_expanded(fit)

  expect_s3_class(p, "ggplot")
  expect_false(is.null(p))
})

# ===== INPUT VALIDATION TESTS =====

test_that("error when takes_em provided without effort_em", {
  expect_error(
    fit_bycatch(Takes_obs ~ 1,
                data = d_multi_pois,
                time = "Year",
                effort = "Sets_obs",
                takes_em = "Takes_em",
                # Missing effort_em
                family = "poisson"
    ),
    "effort_em must also be provided"
  )
})

test_that("error when effort_em provided without takes_em", {
  expect_error(
    fit_bycatch(Takes_obs ~ 1,
                data = d_multi_pois,
                time = "Year",
                effort = "Sets_obs",
                effort_em = "Sets_em",
                # Missing takes_em
                family = "poisson"
    ),
    "takes_em must also be provided"
  )
})

test_that("error when takes_both provided without effort_both", {
  expect_error(
    fit_bycatch(Takes_obs ~ 1,
                data = d_multi_three,
                time = "Year",
                effort = "Sets_obs",
                takes_em = "Takes_em",
                effort_em = "Sets_em",
                takes_both = "Takes_both",
                # Missing effort_both
                family = "poisson"
    ),
    "effort_both must also be provided"
  )
})

test_that("error when column not found", {
  expect_error(
    fit_bycatch(Takes_obs ~ 1,
                data = d_multi_pois,
                time = "Year",
                effort = "Sets_obs",
                takes_em = "NonexistentColumn",
                effort_em = "Sets_em",
                family = "poisson"
    ),
    "NonexistentColumn not found"
  )
})

test_that("warning when observed exceeds total effort", {
  d_bad <- d_multi_pois
  d_bad$Sets_total <- rep(50, sample_size)  # Less than Sets_obs + Sets_em

  expect_warning(
    fit_bycatch(Takes_obs ~ 1,
                data = d_bad,
                time = "Year",
                effort = "Sets_obs",
                takes_em = "Takes_em",
                effort_em = "Sets_em",
                effort_total = "Sets_total",
                family = "poisson",
                iter = 100,
                chains = 1
    ),
    "observed effort exceeding total effort"
  )
})

# ===== COMPARISON TESTS =====

test_that("multi-stream estimates similar rate to single-stream pooled data", {
  set.seed(123)

  # Multi-stream fit
  fit_multi <- fit_bycatch(Takes_obs ~ 1,
                           data = d_multi_pois,
                           time = "Year",
                           effort = "Sets_obs",
                           takes_em = "Takes_em",
                           effort_em = "Sets_em",
                           effort_total = "Sets_total",
                           family = "poisson",
                           time_varying = FALSE,
                           iter = 500,
                           chains = 1
  )

  # Manual pooling for comparison
  d_pooled <- d_multi_pois
  d_pooled$Takes_pooled <- d_pooled$Takes_obs + d_pooled$Takes_em
  d_pooled$Sets_pooled <- d_pooled$Sets_obs + d_pooled$Sets_em

  fit_pooled <- fit_bycatch(Takes_pooled ~ 1,
                            data = d_pooled,
                            time = "Year",
                            effort = "Sets_pooled",
                            family = "poisson",
                            time_varying = FALSE,
                            iter = 500,
                            chains = 1
  )

  # Extract rates
  lambda_multi <- rstan::extract(fit_multi$fitted_model, "lambda_base")$lambda_base
  lambda_pooled <- rstan::extract(fit_pooled$fitted_model, "lambda")$lambda

  # Multi-stream: lambda_base is already a rate per unit effort
  rate_multi <- mean(lambda_multi)

  # Single-stream: lambda is expected count, divide by effort to get rate
  # Use first observation's effort as representative
  rate_pooled <- mean(lambda_pooled[, 1]) / d_pooled$Sets_pooled[1]

  # Should be similar (within 20% given sampling variation)
  expect_equal(rate_multi, rate_pooled, tolerance = 0.2)
})

# ===== STRESS TESTS =====

test_that("multi-stream works with only EM (no observer)", {
  set.seed(123)

  # Set observer to zero
  d_em_only <- d_multi_pois
  d_em_only$Takes_obs <- rep(0, sample_size)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_em_only,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
})

test_that("multi-stream works with sparse data", {
  set.seed(123)

  # Very sparse data
  d_sparse <- d_multi_pois
  d_sparse$Takes_obs <- rep(0, sample_size)
  d_sparse$Takes_obs[c(5, 10, 15)] <- c(1, 2, 1)
  d_sparse$Takes_em <- rep(0, sample_size)
  d_sparse$Takes_em[c(7, 12)] <- c(1, 1)

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_sparse,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")
  expect_true(fit$multi_stream)
})

test_that("multi-stream works with 100% coverage (no expansion)", {
  set.seed(123)

  # Total effort equals observed effort (no unobserved)
  d_full <- d_multi_pois
  d_full$Sets_total <- d_full$Sets_obs + d_full$Sets_em

  fit <- fit_bycatch(Takes_obs ~ 1,
                     data = d_full,
                     time = "Year",
                     effort = "Sets_obs",
                     takes_em = "Takes_em",
                     effort_em = "Sets_em",
                     effort_total = "Sets_total",
                     family = "poisson",
                     time_varying = FALSE,
                     iter = 100,
                     chains = 1
  )

  expect_type(fit, "list")

  # Unobserved should be zero or very small
  expanded <- get_expanded(fit)
  expect_true(mean(rowSums(expanded)) < 1)
})
