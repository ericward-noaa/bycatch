#' fit_bycatch is the primary function for fitting bycatch models to time series of takes and effort
#'
#' @param formula The model formula.
#' @param data A data frame.
#' @param time Named column of the 'data' data frame with the label for the time (e.g. year) variable
#' @param effort Named column of the 'effort' variable in the data frame with the label for the fishing effort to be used in estimation of mean bycatch rate. This
#' represents total observed effort
#' @param expansion_rate The expansion rate to be used in generating distributions for unobserved sets. If NULL, defaults to 100% coverage (= 100). Should be the same as the coverage rate. Deprecated if using effort_total.
#' @param takes_em Optional: Column name for bycatch takes recorded by EM (electronic monitoring). Default NULL.
#' @param effort_em Optional: Column name for effort monitored by EM. Required if takes_em is provided.
#' @param takes_both Optional: Column name for bycatch takes when both Observer and EM are present. Default NULL.
#' @param effort_both Optional: Column name for effort with both Observer and EM monitoring. Required if takes_both is provided.
#' @param effort_total Optional: Column name for total fishery effort. If provided, replaces expansion_rate for calculating unobserved effort.
#' @param family Family for response distribution can be discrete ("poisson",
#' "nbinom2", "poisson-hurdle","nbinom2-hurdle"), or continuous ("normal",
#' "gamma","lognormal", "normal-hurdle", "gamma-hurdle", "lognormal-hurdle"). The
#' default distribution is "poisson". The hurdle variants estimate the
#' probability of zeros (theta) separately from the other models and use
#' truncated distribution to model positive counts. All use a log
#' link function.
#' @param time_varying boolean TRUE/FALSE, whether to include time varying component (this is a random walk, analogous to making this a Dynamic linear model)
#' @param iter the number of mcmc iterations, defaults to 1000
#' @param chains the number of mcmc chains, defaults to 3
#' @param control List to pass to [rstan::sampling()]. For example,
#'   increase \code{adapt_delta} if there are warnings about divergent
#'   transitions: \code{control = list(adapt_delta = 0.99)}. By default,
#'   \pkg{bycatch} sets \code{adapt_delta = 0.9}.
#' @param ... Any other arguments to pass to [rstan::sampling()].
#'
#' @return list of the data used to fit the model, the matrix of covariates, the expanded bycatch generated via the fit and simulations, and the fitted stan model
#'
#' @export
#'
#' @importFrom rstan sampling vb
#' @importFrom cli cli_inform cli_warn
#' @import Rcpp
#' @importFrom stats model.frame model.matrix model.response
#'
#' @details
#' **Multi-stream monitoring:** When using multiple monitoring streams (Observer, EM, Both),
#' the data are kept separate in the statistical model but all share the same underlying bycatch rate.
#' This approach assumes all monitoring types have perfect/equal detection (detection = 100\%) but
#' avoids distributional issues that arise from pooling.
#'
#' @examples
#' \donttest{
#' # Single stream example (backwards compatible)
#' d <- data.frame(
#'   "Year" = 2002:2014,
#'   "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
#'   "expansionRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
#'   "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
#' )
#' fit <- fit_bycatch(Takes ~ 1,
#'   data = d, time = "Year",
#'   effort = "Sets", family = "poisson", time_varying = FALSE
#' )
#'
#' # Multi-stream example (Observer + EM)
#' d_multi <- data.frame(
#'   Year = 2002:2014,
#'   Takes_obs = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0),
#'   Sets_obs = c(200, 180, 150, 350, 250, 270, 180, 150, 400, 350, 280, 180, 250),
#'   Takes_em = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
#'   Sets_em = c(150, 120, 140, 250, 180, 190, 120, 100, 300, 270, 200, 130, 190),
#'   Sets_total = c(1000, 950, 900, 1800, 1200, 1300, 900, 750, 2000, 1800, 1400, 950, 1300)
#' )
#'
#' fit_multi <- fit_bycatch(Takes_obs ~ 1,
#'   data = d_multi,
#'   time = "Year",
#'   effort = "Sets_obs",
#'   takes_em = "Takes_em",
#'   effort_em = "Sets_em",
#'   effort_total = "Sets_total",
#'   family = "poisson"
#' )
#' }
fit_bycatch <- function(formula, data, time = "year", effort = "effort", expansion_rate = NULL,
                        takes_em = NULL,
                        effort_em = NULL,
                        takes_both = NULL,
                        effort_both = NULL,
                        effort_total = NULL,
                        family = c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle",
                                   "lognormal", "gamma", "lognormal-hurdle", "gamma-hurdle",
                                   "normal", "normal-hurdle"),
                        time_varying = FALSE,
                        iter = 1000,
                        chains = 3,
                        control = list(adapt_delta = 0.9, max_treedepth = 20), ...) {

  # Extract formula components
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)

  # Validate basic inputs
  if (time %in% colnames(data) == FALSE) {
    stop("The time variable needs to be specified as a named column in the data frame")
  }
  if (effort %in% colnames(data) == FALSE) {
    stop("The effort variable needs to be specified as a named column in the data frame")
  }
  if (family %in% c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle",
                    "lognormal", "gamma", "lognormal-hurdle", "gamma-hurdle",
                    "normal", "normal-hurdle") == FALSE) {
    stop("The family must be specified as a distribution in the recognized list")
  }

  # Validate multi-stream inputs (ALWAYS check, regardless of mode)
  if (!is.null(takes_em) && is.null(effort_em)) {
    stop("If takes_em is provided, effort_em must also be provided")
  }
  if (!is.null(effort_em) && is.null(takes_em)) {
    stop("takes_em must also be provided")
  }
  if (!is.null(takes_both) && is.null(effort_both)) {
    stop("If takes_both is provided, effort_both must also be provided")
  }
  if (!is.null(effort_both) && is.null(takes_both)) {
    stop("If effort_both is provided, takes_both must also be provided")
  }

  # Check if multi-stream mode
  multi_stream <- !is.null(takes_em) || !is.null(takes_both)

  # Additional multi-stream validation
  if (multi_stream) {
    if (!is.null(effort_total) && !is.null(expansion_rate)) {
      cli_warn("Both effort_total and expansion_rate provided. Using effort_total and ignoring expansion_rate.")
    }

    # Check that multi-stream columns exist
    if (!is.null(takes_em) && !(takes_em %in% colnames(data))) {
      stop(paste("Column", takes_em, "not found in data"))
    }
    if (!is.null(effort_em) && !(effort_em %in% colnames(data))) {
      stop(paste("Column", effort_em, "not found in data"))
    }
    if (!is.null(takes_both) && !(takes_both %in% colnames(data))) {
      stop(paste("Column", takes_both, "not found in data"))
    }
    if (!is.null(effort_both) && !(effort_both %in% colnames(data))) {
      stop(paste("Column", effort_both, "not found in data"))
    }
    if (!is.null(effort_total) && !(effort_total %in% colnames(data))) {
      stop(paste("Column", effort_total, "not found in data"))
    }
  }

  # Get response variable
  y <- as.numeric(model.response(mf, "numeric"))

  # Determine family ID
  family_id <- match(family, c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle",
                               "lognormal", "gamma", "lognormal-hurdle", "gamma-hurdle",
                               "normal", "normal-hurdle"))

  # Prepare data list for Stan
  if (multi_stream) {

    # Prepare time index
    time_values <- data[[time]]
    time_min <- min(time_values)
    time_idx <- time_values - time_min + 1
    n_year <- length(seq(min(time_values), max(time_values)))

    # Get observer data (from formula)
    takes_obs_vec <- y
    effort_obs_vec <- data[[effort]]
    n_obs <- length(takes_obs_vec)
    time_idx_obs <- time_idx

    cli_inform("Observer stream: {n_obs} observations")
    cli_inform("  Total takes: {sum(takes_obs_vec)}")
    cli_inform("  Total effort: {sum(effort_obs_vec)}")

    # Get EM data (if provided)
    if (!is.null(takes_em)) {
      takes_em_vec <- data[[takes_em]]
      effort_em_vec <- data[[effort_em]]
      n_em <- length(takes_em_vec)
      time_idx_em <- time_idx

      cli_inform("EM stream: {n_em} observations")
      cli_inform("  Total takes: {sum(takes_em_vec)}")
      cli_inform("  Total effort: {sum(effort_em_vec)}")
    } else {
      takes_em_vec <- numeric(0)
      effort_em_vec <- numeric(0)
      n_em <- 0
      time_idx_em <- integer(0)
    }

    # Get Both data (if provided)
    if (!is.null(takes_both)) {
      takes_both_vec <- data[[takes_both]]
      effort_both_vec <- data[[effort_both]]
      n_both <- length(takes_both_vec)
      time_idx_both <- time_idx

      cli_inform("Both stream: {n_both} observations")
      cli_inform("  Total takes: {sum(takes_both_vec)}")
      cli_inform("  Total effort: {sum(effort_both_vec)}")
    } else {
      takes_both_vec <- numeric(0)
      effort_both_vec <- numeric(0)
      n_both <- 0
      time_idx_both <- integer(0)
    }

    # Calculate total observed effort by year
    effort_by_year <- rep(0, n_year)
    for(i in 1:n_obs) {
      t <- time_idx_obs[i]
      effort_by_year[t] <- effort_by_year[t] + effort_obs_vec[i]
    }
    if(n_em > 0) {
      for(i in 1:n_em) {
        t <- time_idx_em[i]
        effort_by_year[t] <- effort_by_year[t] + effort_em_vec[i]
      }
    }
    if(n_both > 0) {
      for(i in 1:n_both) {
        t <- time_idx_both[i]
        effort_by_year[t] <- effort_by_year[t] + effort_both_vec[i]
      }
    }

    # Calculate unobserved effort by year
    if (!is.null(effort_total)) {
      # Use total effort directly
      total_effort_by_year <- rep(0, n_year)
      for(i in 1:nrow(data)) {
        t <- time_idx[i]
        total_effort_by_year[t] <- data[[effort_total]][i]
      }
      new_effort_by_year <- total_effort_by_year - effort_by_year

      # Check for negative unobserved effort
      if (any(new_effort_by_year < 0, na.rm=T)) {
        cli_warn("Some years have observed effort exceeding total effort. Setting unobserved effort to 0 for these years.")
        new_effort_by_year[new_effort_by_year < 0] <- 0
      }

      if (any(is.na(new_effort_by_year))) {
        new_effort_by_year[which(is.na(new_effort_by_year))] <- 0
      }

      cli_inform("Total fishery effort: {sum(total_effort_by_year)}")
      cli_inform("Observed effort: {sum(effort_by_year)}")
      cli_inform("Unobserved effort: {sum(new_effort_by_year)}")

    } else if (!is.null(expansion_rate)) {
      # Use expansion rate (backwards compatible)
      total_effort_by_year <- rep(0, n_year)
      for(i in 1:nrow(data)) {
        t <- time_idx[i]
        p_expansion <- data[[expansion_rate]][i] / 100
        total_effort_i <- ceiling(effort_by_year[t] / p_expansion)
        total_effort_by_year[t] <- max(total_effort_by_year[t], total_effort_i)
      }
      new_effort_by_year <- total_effort_by_year - effort_by_year

    } else {
      # No expansion (100% coverage)
      new_effort_by_year <- rep(0, n_year)
      cli_inform("No expansion (100% coverage assumed)")
    }

    # Determine if discrete or continuous
    if (family %in% c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle")) {
      # Discrete family
      yint_obs <- as.integer(takes_obs_vec)
      yint_em <- if(n_em > 0) as.integer(takes_em_vec) else integer(0)
      yint_both <- if(n_both > 0) as.integer(takes_both_vec) else integer(0)
      yreal_obs <- takes_obs_vec  # Not used but need for Stan
      yreal_em <- if(n_em > 0) takes_em_vec else numeric(0)
      yreal_both <- if(n_both > 0) takes_both_vec else numeric(0)
    } else {
      # Continuous family
      yint_obs <- ifelse(takes_obs_vec > 0, 1L, 0L)
      yint_em <- if(n_em > 0) ifelse(takes_em_vec > 0, 1L, 0L) else integer(0)
      yint_both <- if(n_both > 0) ifelse(takes_both_vec > 0, 1L, 0L) else integer(0)
      yreal_obs <- takes_obs_vec
      yreal_em <- if(n_em > 0) takes_em_vec else numeric(0)
      yreal_both <- if(n_both > 0) takes_both_vec else numeric(0)
    }

    # For multi-stream mode, we need year-level covariates (n_year x K)
    # For intercept-only model, this is just a column of 1s
    # For models with covariates, aggregate by year (use first value for each year)
    X_year <- matrix(nrow = n_year, ncol = ncol(X))
    for(k in 1:ncol(X)) {
      for(t in 1:n_year) {
        # Get first observation for this year
        idx <- which(time_idx == t)[1]
        X_year[t, k] <- X[idx, k]
      }
    }

    # Prepare data list for Stan
    datalist <- list(
      # Multi-stream indicator
      multi_stream = 1,

      # Stream counts - use actual counts (can be 0)
      n_obs = n_obs,
      n_em = n_em,
      n_both = n_both,

      # Stream data - discrete
      yint_obs = yint_obs,
      yint_em = yint_em,
      yint_both = yint_both,

      # Stream data - continuous
      yreal_obs = yreal_obs,
      yreal_em = yreal_em,
      yreal_both = yreal_both,

      # Stream effort
      effort_obs = effort_obs_vec,
      effort_em = if(n_em > 0) effort_em_vec else numeric(0),
      effort_both = if(n_both > 0) effort_both_vec else numeric(0),

      # Time indices
      time_idx_obs = time_idx_obs,
      time_idx_em = if(n_em > 0) time_idx_em else integer(0),
      time_idx_both = if(n_both > 0) time_idx_both else integer(0),

      # Unobserved effort by year
      new_effort_by_year = new_effort_by_year,

      # Backwards compatibility (not used in multi-stream mode)
      n_row = nrow(data),
      effort = data[[effort]],
      new_effort = rep(0, nrow(data)),  # Not used
      yint = as.integer(y),
      yreal = y,
      time = time_idx,

      # Common data - use year level covariates
      n_year = n_year,
      K = ncol(X_year),
      x = X_year,  # Changed from X to X_year
      family = family_id,
      time_varying = as.numeric(time_varying)
    )

    # Store stream info for later use
    stream_info <- list(
      "takes_em" = takes_em,
      "effort_em" = effort_em,
      "takes_both" = takes_both,
      "effort_both" = effort_both,
      "effort_total" = effort_total,
      "n_obs" = n_obs,
      "n_em" = n_em,
      "n_both" = n_both
    )

  } else {
    # Single stream mode (backwards compatible)

    observed_effort <- data[[effort]]

    # Calculate unobserved effort (original logic)
    if (is.null(expansion_rate)) {
      p_expansion <- rep(1, nrow(data))
    } else {
      p_expansion <- data[[expansion_rate]] / 100
    }

    total_effort <- ceiling(observed_effort / p_expansion)
    new_effort <- ceiling(total_effort * (1 - p_expansion))

    # Prepare response variable
    if (family %in% c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle") == FALSE) {
      yint <- ifelse(y > 0, 1, 0)
    } else {
      yint <- y
    }

    # Calculate new_effort by year
    time_values <- data[[time]]
    time_min <- min(time_values)
    time_idx <- time_values - time_min + 1
    n_year <- length(seq(min(time_values), max(time_values)))

    new_effort_by_year <- rep(0, n_year)
    for(i in 1:nrow(data)) {
      t <- time_idx[i]
      new_effort_by_year[t] <- new_effort_by_year[t] + new_effort[i]
    }

    # Create year-level covariate matrix (same as multi-stream)
    X_year <- matrix(nrow = n_year, ncol = ncol(X))
    for(k in 1:ncol(X)) {
      for(t in 1:n_year) {
        # Get first observation for this year
        idx <- which(time_idx == t)[1]
        X_year[t, k] <- X[idx, k]
      }
    }

    datalist <- list(
      # Single stream mode
      multi_stream = 0,

      # Zero-length arrays for multi-stream placeholders
      n_obs = 0,
      n_em = 0,
      n_both = 0,
      yint_obs = integer(0),
      yint_em = integer(0),
      yint_both = integer(0),
      yreal_obs = numeric(0),
      yreal_em = numeric(0),
      yreal_both = numeric(0),
      effort_obs = numeric(0),
      effort_em = numeric(0),
      effort_both = numeric(0),
      time_idx_obs = integer(0),
      time_idx_em = integer(0),
      time_idx_both = integer(0),
      new_effort_by_year = new_effort_by_year,

      # Single stream data
      n_row = nrow(data),
      effort = observed_effort,
      new_effort = new_effort,
      yint = yint,
      yreal = y,
      time = time_idx,
      n_year = n_year,
      K = ncol(X_year),
      x = X_year,
      family = family_id,
      time_varying = as.numeric(time_varying)
    )

    stream_info <- NULL
  }

  # Define parameters to extract
  if(multi_stream) {
    pars <- c("beta", "lambda_base", "log_lambda_base", "log_lik", "y_new",
              "est_time_dev", "sigma_rw", "sigma_logn", "cv_gamma",
              "nb2_phi", "theta", "y_new_real")
  } else {
    pars <- c("beta", "lambda", "log_lambda", "log_lik", "y_new",
              "est_time_dev", "sigma_rw", "sigma_logn", "cv_gamma",
              "nb2_phi", "theta", "y_new_real")
  }

  # Fit Stan model
  mod <- rstan::sampling(
    object = stanmodels$bycatch,
    data = datalist,
    pars = pars,
    iter = iter,
    chains = chains,
    control = control, ...
  )

  # Prepare return object
  result <- list(
    "data" = data,
    "effort" = effort,
    "events" = names(mf)[1],
    "time" = time,
    "fitted_model" = mod,
    "family" = family,
    "multi_stream" = multi_stream,
    "stream_info" = stream_info
  )

  return(result)
}
