context("model fitting tests")

# this controls how large to make the data below
sample_size <- 20

set.seed(123)
# simulate data
d_nb <- data.frame(
  "Year" = 1:sample_size,
  "Takes" = rnbinom(sample_size, size = 10, prob = 0.9),
  "expansionRate" = 100,
  "Sets" = 100
)

# simulate data
set.seed(123)
d_pois <- data.frame(
  "Year" = 1:sample_size,
  "Takes" = rnbinom(sample_size, size = 10, prob = 0.9),
  "expansionRate" = 100,
  "Sets" = 100
)

test_that("fitting function works for poisson model", {
  set.seed(123)
  # simulate data
  fit <- fit_bycatch(Takes ~ 1,
    data = d_pois, time = "Year",
    effort = "Sets", family = "poisson", time_varying = FALSE
  )

  pars <- rstan::extract(fit$fitted_model, "lambda")
  expect_equal(apply(pars$lambda, 2, mean), rep(1.367753, sample_size), tol = 0.1)
})


test_that("fitting function works for negative binomial2 model", {
  set.seed(123)
  fit <- fit_bycatch(Takes ~ 1,
    data = d_nb, time = "Year",
    effort = "Sets", family = "nbinom2", time_varying = FALSE
  )

  pars <- rstan::extract(fit$fitted_model, "lambda")
  expect_equal(apply(pars$lambda, 2, mean), rep(1.401189, sample_size), tol = 0.1)
})

test_that("fitting function works for hurdle poisson model", {
  set.seed(123)

  fit <- fit_bycatch(Takes ~ 1,
    data = d_pois, time = "Year",
    effort = "Sets", family = "poisson-hurdle", time_varying = FALSE
  )

  pars <- rstan::extract(fit$fitted_model, "lambda")
  expect_equal(apply(pars$lambda, 2, mean), rep(1.353893, sample_size), tol = 0.1)
})

# Commented out for speed
# test_that("fitting function works for hurdle negative binomial2 model", {
#   set.seed(123)
#   fit = fit_bycatch(Takes ~ 1, data=d_nb,iter=100,chains=1, time="Year",
#                     effort="Sets", family="nbinom2-hurdle", time_varying = FALSE)
#
#   pars = rstan::extract(fit$fitted_model, "lambda")
#   expect_equal(apply(pars$lambda,2,mean), rep(3.485921,sample_size), tol = 0.001)
# })


test_that("fitting function works for time varying poisson model", {
  set.seed(123)

  fit <- fit_bycatch(Takes ~ 1,
    data = d_pois, time = "Year",
    effort = "Sets", family = "poisson", time_varying = TRUE,
    iter = 100, chains = 1
  )

  # pars <- rstan::extract(fit$fitted_model, "lambda")
  # expect_equal(apply(pars$lambda, 2, mean)[1:5], c(2.591904, 1.691243, 1.730186, 1.729729, 1.580735), tol = 0.001)
  expect_type(fit, "list")
})

test_that("fitting function works for time varying nbinom2 model", {
  set.seed(123)
  # simulate data
  fit <- fit_bycatch(Takes ~ 1,
    data = d_nb, time = "Year",
    effort = "Sets", family = "nbinom2", time_varying = TRUE,
    iter = 100, chains = 1
  )

  # pars <- rstan::extract(fit$fitted_model, "lambda")
  # expect_equal(apply(pars$lambda, 2, mean)[1:5], c(3.165161, 2.258439, 2.151508, 2.364193, 2.560882), tol = 0.2)
  expect_type(fit, "list")
})

test_that("fitting function works for time varying poisson hurdle model", {
  set.seed(123)

  fit <- fit_bycatch(Takes ~ 1,
    data = d_pois, time = "Year",
    effort = "Sets", family = "poisson-hurdle", time_varying = TRUE,
    iter = 100, chains = 1
  )
  expect_type(fit, "list")
  # pars <- rstan::extract(fit$fitted_model, "lambda")
  # expect_equal(apply(pars$lambda, 2, mean)[1:5], c(2.101424, 1.401246, 1.305310, 1.370120, 1.645385), tol = 0.001)
})

test_that("fitting function works for time varying nbinom2 hurdle model", {
  set.seed(123)

  fit <- fit_bycatch(Takes ~ 1,
    data = d_pois, time = "Year",
    effort = "Sets", family = "nbinom2-hurdle", time_varying = TRUE,
    iter = 100, chains = 1
  )
  expect_type(fit, "list")
  # pars <- rstan::extract(fit$fitted_model, "lambda")
  # expect_equal(apply(pars$lambda, 2, mean)[1:5], c(0.5550384, 2.1295139, 2.2619993, 1.9444669, 1.7304549), tol = 0.001)
})
