#' Aggregate self-reporting probability sub-model inputs to the year level
#'
#' Internal helper. `vessels_observed` gives, per row, the number of observed
#' (obs/em/both-covered) vessels; `covrate_selfreport` gives, per row, the
#' percentage of those vessels that also self-report (0-100). Together these
#' reconstruct a binomial numerator/denominator, which the Stan model uses to
#' estimate a single pooled p_report with a posterior (rather than trusting a
#' fixed rate the way covrate_obs/covrate_em are used).
#'
#' @param data The original data frame passed to fit_bycatch.
#' @param time_idx Integer vector mapping each row of data to a 1:n_year index.
#' @param n_year Number of unique years.
#' @param vessels_observed Column name (or NULL) for observed-vessel counts (denominator).
#' @param covrate_selfreport Column name (or NULL) for self-report coverage rate,
#'   0-100 (numerator, as a percentage of vessels_observed).
#'
#' @return list of Stan data fields: estimate_self_report,
#'   n_vessels_observed, n_vessels_self_report
#' @noRd
.prepare_self_report_data <- function(data, time_idx, n_year,
                                       vessels_observed, covrate_selfreport) {

  estimate_self_report <- !is.null(vessels_observed) && !is.null(covrate_selfreport)

  n_vessels_observed <- rep(0L, n_year)
  n_vessels_self_report <- rep(0L, n_year)

  if (estimate_self_report) {
    row_n_observed <- data[[vessels_observed]]
    row_n_self_report <- round(data[[covrate_selfreport]] / 100 * row_n_observed)

    for (i in seq_len(nrow(data))) {
      t <- time_idx[i]
      n_vessels_observed[t] <- n_vessels_observed[t] + row_n_observed[i]
      n_vessels_self_report[t] <- n_vessels_self_report[t] + row_n_self_report[i]
    }
    cli_inform("Self-report sub-model: {sum(n_vessels_observed)} observed vessels, {sum(n_vessels_self_report)} implied self-reporters")
  }

  list(
    estimate_self_report = as.integer(estimate_self_report),
    n_vessels_observed = as.integer(n_vessels_observed),
    n_vessels_self_report = as.integer(n_vessels_self_report)
  )
}

#' fit_bycatch - Unified version with coverage rate parameters
#'
#' @param formula The model formula.
#' @param data A data frame.
#' @param time Named column of the 'data' data frame with the label for the time (e.g. year) variable
#' @param effort Named column of the 'effort' variable in the data frame with the label for the fishing effort to be used in estimation of mean bycatch rate. This represents total observed effort
#' @param covrate Optional: Column name for coverage rate (percentage) for single-stream models. Replaces expansion_rate.
#' @param expansion_rate DEPRECATED: Use 'covrate' or 'covrate_obs' instead. The expansion rate to be used in generating distributions for unobserved sets.
#' @param takes_em Optional: Column name for bycatch takes recorded by EM (electronic monitoring). Default NULL.
#' @param effort_em Optional: Column name for effort monitored by EM. Required if takes_em is provided.
#' @param covrate_em Optional: Column name for EM coverage rate (percentage). Used to calculate unobserved effort.
#' @param takes_both Optional: Column name for bycatch takes when both Observer and EM are present. Default NULL.
#' @param effort_both Optional: Column name for effort with both Observer and EM monitoring. Required if takes_both is provided.
#' @param covrate_obs Optional: Column name for Observer coverage rate (percentage). Used to calculate unobserved effort.
#' @param covrate_both Optional: Column name for "Both" coverage rate (percentage). Used to calculate unobserved effort.
#' @param effort_total Optional: Column name for total fishery effort. If provided, takes precedence over coverage rates. Must be in same units as effort columns.
#' @param takes_selfreport Optional: Column name for bycatch takes recorded via self-report. Default NULL. Treated as a 4th monitoring stream, exactly like 'takes_em'/'takes_both' - contributes directly to the shared bycatch rate.
#' @param effort_selfreport Optional: Column name for effort covered by self-report. Required if takes_selfreport is provided.
#' @param covrate_selfreport Optional: Column name for the percentage (0-100) of observed (obs/em/both-covered) vessels that also self-report, per row. Used with 'vessels_observed' to estimate the probability of self-reporting, p_report, with a full posterior (unlike covrate_obs/covrate_em/covrate_both, which are treated as fixed known rates).
#' @param vessels_observed Optional: Column name for the count of observed (obs/em/both-covered) vessels per row - the denominator for the p_report binomial sub-model. Required together with 'covrate_selfreport'.
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
#' **Coverage Rates vs effort_total:**
#'
#' This function supports two approaches for calculating unobserved effort:
#'
#' 1. **Coverage rates** (recommended): Use `covrate`, `covrate_obs`, `covrate_em`, `covrate_both`
#'    to specify what percentage of the fishery is monitored. The function calculates
#'    unobserved effort as: `observed_effort * ((100 - coverage) / coverage)`
#'
#' 2. **Total effort**: Use `effort_total` to provide the total fishery effort directly.
#'    Unobserved effort is calculated as: `effort_total - observed_effort`
#'    NOTE: effort_total must be in the same units as your effort columns
#'
#' **Priority order**: effort_total > coverage rates > 100% coverage assumed
#'
#' **Multi-stream monitoring:** When using multiple monitoring streams (Observer, EM, Both),
#' the data are kept separate in the statistical model but all share the same underlying bycatch rate.
#' This approach assumes all monitoring types have perfect/equal detection (detection = 100%) but
#' avoids distributional issues that arise from pooling.
#'
#' **Self-reporting as a monitoring stream:** Providing `takes_selfreport`/`effort_selfreport`
#' adds self-report as a 4th stream, structurally identical to `em`/`both` - its catches feed
#' directly into the shared bycatch rate.
#'
#' **Self-reporting probability sub-model:** A vessel can be covered by EM/Observer *and*
#' self-report, so summing streams naively can double-count effort. Providing
#' `covrate_selfreport` (percentage of observed vessels that also self-report, per row) and
#' `vessels_observed` (the denominator) fits a binomial sub-model for the probability of
#' self-reporting, p_report (pooled across all years). Within each posterior draw, p_report
#' is used to estimate how much of the self-report stream's effort overlaps with vessels
#' already covered by obs/em/both; only the genuinely incremental portion is netted out of
#' the unobserved effort before expansion — so the posterior uncertainty in the self-report
#' rate shows up directly in the width of the expanded bycatch estimates, rather than being
#' applied as a fixed correction factor after fitting.
#'
#' @examples
#' \donttest{
#' # Single stream example with coverage rate
#' d <- data.frame(
#'   "Year" = 2002:2014,
#'   "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
#'   "CovRate" = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
#'   "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
#' )
#' fit <- fit_bycatch(Takes ~ 1,
#'   data = d, time = "Year",
#'   effort = "Sets",
#'   covrate = "CovRate",
#'   family = "poisson",
#'   time_varying = FALSE
#' )
#'
#' # Multi-stream example with coverage rates (Observer + EM)
#' d_multi <- data.frame(
#'   Year = 2002:2014,
#'   Takes_obs = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0),
#'   Sets_obs = c(200, 180, 150, 350, 250, 270, 180, 150, 400, 350, 280, 180, 250),
#'   CovRate_obs = c(20, 19, 17, 19, 21, 21, 20, 20, 20, 19, 20, 19, 19),
#'   Takes_em = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
#'   Sets_em = c(150, 120, 140, 250, 180, 190, 120, 100, 300, 270, 200, 130, 190),
#'   CovRate_em = c(15, 13, 16, 14, 15, 15, 13, 13, 15, 15, 14, 14, 15)
#' )
#'
#' fit_multi <- fit_bycatch(Takes_obs ~ 1,
#'   data = d_multi,
#'   time = "Year",
#'   effort = "Sets_obs",
#'   covrate_obs = "CovRate_obs",
#'   takes_em = "Takes_em",
#'   effort_em = "Sets_em",
#'   covrate_em = "CovRate_em",
#'   family = "poisson"
#' )
#'
#' # Multi-stream example with self-report as a 4th stream plus a reporting-probability sub-model
#' d_sr <- data.frame(
#'   Year = 2002:2014,
#'   Takes_obs = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0),
#'   Sets_obs = c(200, 180, 150, 350, 250, 270, 180, 150, 400, 350, 280, 180, 250),
#'   CovRate_obs = c(20, 19, 17, 19, 21, 21, 20, 20, 20, 19, 20, 19, 19),
#'   # Self-report stream: its own catches and effort, like takes_em/effort_em
#'   Takes_selfreport = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
#'   Sets_selfreport = c(60, 55, 45, 80, 65, 70, 48, 40, 90, 78, 62, 50, 58),
#'   # Of the observed (obs/em/both) vessels each year, how many also self-report?
#'   VesselsObserved = c(18, 17, 15, 22, 20, 21, 17, 15, 24, 22, 20, 17, 19),
#'   CovRate_selfreport = c(67, 65, 60, 68, 70, 71, 65, 67, 71, 73, 70, 65, 68) # percent
#' )
#'
#' fit_sr <- fit_bycatch(Takes_obs ~ 1,
#'   data = d_sr,
#'   time = "Year",
#'   effort = "Sets_obs",
#'   covrate_obs = "CovRate_obs",
#'   takes_selfreport = "Takes_selfreport",
#'   effort_selfreport = "Sets_selfreport",
#'   covrate_selfreport = "CovRate_selfreport",
#'   vessels_observed = "VesselsObserved",
#'   family = "poisson"
#' )
#'
#' # Posterior of the self-report probability
#' p_report_draws <- rstan::extract(fit_sr$fitted_model, "p_report")$p_report
#' }
fit_bycatch <- function(formula, data, time = "year", effort = "effort",
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
                        takes_selfreport = NULL,
                        effort_selfreport = NULL,
                        covrate_selfreport = NULL,
                        vessels_observed = NULL,
                        family = c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle",
                                   "lognormal", "gamma", "lognormal-hurdle", "gamma-hurdle",
                                   "normal", "normal-hurdle"),
                        time_varying = FALSE,
                        iter = 1000,
                        chains = 3,
                        control = list(adapt_delta = 0.9, max_treedepth = 20), ...) {

  # Handle backward compatibility for expansion_rate
  if (!is.null(expansion_rate)) {
    cli_warn("'expansion_rate' parameter is deprecated. Please use 'covrate' for single-stream models or 'covrate_obs' for multi-stream models.")

    # Map expansion_rate to appropriate coverage parameter
    if (is.null(covrate) && is.null(covrate_obs)) {
      # Determine if this is multi-stream
      if (!is.null(takes_em) || !is.null(takes_both)) {
        covrate_obs <- expansion_rate
      } else {
        covrate <- expansion_rate
      }
    }
  }

  # For single-stream models, covrate and covrate_obs are equivalent
  if (!is.null(covrate) && is.null(covrate_obs) && is.null(takes_em) && is.null(takes_both)) {
    covrate_obs <- covrate
  }

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

  # Validate multi-stream inputs
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

  # Validate self-report stream inputs (a 4th stream, like em/both)
  if (!is.null(takes_selfreport) && is.null(effort_selfreport)) {
    stop("If takes_selfreport is provided, effort_selfreport must also be provided")
  }
  if (!is.null(effort_selfreport) && is.null(takes_selfreport)) {
    stop("If effort_selfreport is provided, takes_selfreport must also be provided")
  }

  # Validate self-report probability sub-model inputs
  if (!is.null(covrate_selfreport) && is.null(vessels_observed)) {
    stop("If covrate_selfreport is provided, vessels_observed must also be provided")
  }
  if (!is.null(vessels_observed) && is.null(covrate_selfreport)) {
    stop("If vessels_observed is provided, covrate_selfreport must also be provided")
  }
  if (!is.null(covrate_selfreport) && !(covrate_selfreport %in% colnames(data))) {
    stop(paste("Column", covrate_selfreport, "not found in data"))
  }
  if (!is.null(vessels_observed) && !(vessels_observed %in% colnames(data))) {
    stop(paste("Column", vessels_observed, "not found in data"))
  }
  if (!is.null(covrate_selfreport)) {
    if (any(data[[covrate_selfreport]] < 0 | data[[covrate_selfreport]] > 100, na.rm = TRUE)) {
      stop("covrate_selfreport must be a percentage between 0 and 100")
    }
  }

  # The self-report probability sub-model piggybacks on the self-report stream
  # (it estimates overlap between that stream and obs/em/both), so require both together
  if (!is.null(covrate_selfreport) && is.null(takes_selfreport)) {
    stop("covrate_selfreport/vessels_observed require takes_selfreport/effort_selfreport to also be provided")
  }

  # Check if multi-stream mode
  multi_stream <- !is.null(takes_em) || !is.null(takes_both) || !is.null(takes_selfreport)

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

    # Self-report probability sub-model inputs (aggregated to year level)
    self_report_data <- .prepare_self_report_data(
      data = data, time_idx = time_idx, n_year = n_year,
      vessels_observed = vessels_observed, covrate_selfreport = covrate_selfreport
    )

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

    # Get Self-report data (if provided) - a 4th stream, same shape as EM/Both
    if (!is.null(takes_selfreport)) {
      takes_selfreport_vec <- data[[takes_selfreport]]
      effort_selfreport_vec <- data[[effort_selfreport]]
      n_selfreport <- length(takes_selfreport_vec)
      time_idx_selfreport <- time_idx

      cli_inform("Self-report stream: {n_selfreport} observations")
      cli_inform("  Total takes: {sum(takes_selfreport_vec)}")
      cli_inform("  Total effort: {sum(effort_selfreport_vec)}")
    } else {
      takes_selfreport_vec <- numeric(0)
      effort_selfreport_vec <- numeric(0)
      n_selfreport <- 0
      time_idx_selfreport <- integer(0)
    }

    # Calculate total observed effort by year (obs+em+both only - self-report
    # is netted against this inside Stan via p_report to avoid double-counting
    # effort from vessels that are both monitored AND self-report)
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
    # Priority: effort_total > coverage rates > 100% coverage

    if (!is.null(effort_total)) {
      # Original approach: Use external effort_total column
      cli_inform("Using effort_total column for expansion")

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

    } else if (!is.null(covrate_obs) || !is.null(covrate_em) || !is.null(covrate_both)) {
      # Coverage rate approach: Calculate expansion from coverage percentages
      cli_inform("Using coverage rates for expansion")

      # Calculate total coverage rate by year
      total_covrate_by_year <- rep(0, n_year)

      # Add observer coverage
      if (!is.null(covrate_obs)) {
        if (!(covrate_obs %in% colnames(data))) {
          stop(paste("Column", covrate_obs, "not found in data"))
        }
        for(i in 1:n_obs) {
          t <- time_idx_obs[i]
          total_covrate_by_year[t] <- total_covrate_by_year[t] + data[[covrate_obs]][i]
        }
      }

      # Add EM coverage
      if (!is.null(covrate_em)) {
        if (!(covrate_em %in% colnames(data))) {
          stop(paste("Column", covrate_em, "not found in data"))
        }
        if (n_em > 0) {
          for(i in 1:n_em) {
            t <- time_idx_em[i]
            total_covrate_by_year[t] <- total_covrate_by_year[t] + data[[covrate_em]][i]
          }
        }
      }

      # Add Both coverage
      if (!is.null(covrate_both)) {
        if (!(covrate_both %in% colnames(data))) {
          stop(paste("Column", covrate_both, "not found in data"))
        }
        if (n_both > 0) {
          for(i in 1:n_both) {
            t <- time_idx_both[i]
            total_covrate_by_year[t] <- total_covrate_by_year[t] + data[[covrate_both]][i]
          }
        }
      }

      # Calculate new_effort using coverage rates
      # Formula: new_effort = observed_effort * ((100 - coverage_rate) / coverage_rate)
      new_effort_by_year <- rep(0, n_year)
      for(t in 1:n_year) {
        if (total_covrate_by_year[t] > 0 && total_covrate_by_year[t] < 100) {
          # Calculate unobserved proportion
          unobserved_pct <- 100 - total_covrate_by_year[t]
          new_effort_by_year[t] <- effort_by_year[t] * (unobserved_pct / total_covrate_by_year[t])
        } else if (total_covrate_by_year[t] >= 100) {
          new_effort_by_year[t] <- 0
        } else {
          # Coverage is 0 - no information to expand
          new_effort_by_year[t] <- 0
          if (effort_by_year[t] > 0) {
            cli_warn("Year {time_min + t - 1}: Coverage is 0% but effort > 0. Cannot calculate unobserved effort.")
          }
        }
      }

    } else {
      # No expansion information provided
      cli_inform("No expansion information provided - assuming 100% coverage")
      new_effort_by_year <- rep(0, n_year)
    }

    # Prepare response arrays for Stan
    if (family %in% c("poisson", "nbinom2", "poisson-hurdle", "nbinom2-hurdle")) {
      # Discrete family
      yint_obs <- as.integer(takes_obs_vec)
      yint_em <- if(n_em > 0) as.integer(takes_em_vec) else integer(0)
      yint_both <- if(n_both > 0) as.integer(takes_both_vec) else integer(0)
      yint_selfreport <- if(n_selfreport > 0) as.integer(takes_selfreport_vec) else integer(0)
      yreal_obs <- takes_obs_vec
      yreal_em <- if(n_em > 0) takes_em_vec else numeric(0)
      yreal_both <- if(n_both > 0) takes_both_vec else numeric(0)
      yreal_selfreport <- if(n_selfreport > 0) takes_selfreport_vec else numeric(0)
    } else {
      # Continuous family
      yint_obs <- ifelse(takes_obs_vec > 0, 1L, 0L)
      yint_em <- if(n_em > 0) ifelse(takes_em_vec > 0, 1L, 0L) else integer(0)
      yint_both <- if(n_both > 0) ifelse(takes_both_vec > 0, 1L, 0L) else integer(0)
      yint_selfreport <- if(n_selfreport > 0) ifelse(takes_selfreport_vec > 0, 1L, 0L) else integer(0)
      yreal_obs <- takes_obs_vec
      yreal_em <- if(n_em > 0) takes_em_vec else numeric(0)
      yreal_both <- if(n_both > 0) takes_both_vec else numeric(0)
      yreal_selfreport <- if(n_selfreport > 0) takes_selfreport_vec else numeric(0)
    }

    # For multi-stream mode, we need year-level covariates (n_year x K)
    X_year <- matrix(nrow = n_year, ncol = ncol(X))
    for(k in 1:ncol(X)) {
      for(t in 1:n_year) {
        idx <- which(time_idx == t)[1]
        X_year[t, k] <- X[idx, k]
      }
    }

    # Prepare data list for Stan
    datalist <- list(
      multi_stream = 1,
      n_obs = n_obs,
      n_em = n_em,
      n_both = n_both,
      n_selfreport = n_selfreport,
      yint_obs = yint_obs,
      yint_em = yint_em,
      yint_both = yint_both,
      yint_selfreport = yint_selfreport,
      yreal_obs = yreal_obs,
      yreal_em = yreal_em,
      yreal_both = yreal_both,
      yreal_selfreport = yreal_selfreport,
      effort_obs = effort_obs_vec,
      effort_em = if(n_em > 0) effort_em_vec else numeric(0),
      effort_both = if(n_both > 0) effort_both_vec else numeric(0),
      effort_selfreport = if(n_selfreport > 0) effort_selfreport_vec else numeric(0),
      time_idx_obs = time_idx_obs,
      time_idx_em = if(n_em > 0) time_idx_em else integer(0),
      time_idx_both = if(n_both > 0) time_idx_both else integer(0),
      time_idx_selfreport = if(n_selfreport > 0) time_idx_selfreport else integer(0),
      new_effort_by_year = new_effort_by_year,
      n_row = nrow(data),
      effort = data[[effort]],
      new_effort = rep(0, nrow(data)),
      yint = as.integer(y),
      yreal = y,
      time = time_idx,
      n_year = n_year,
      K = ncol(X_year),
      x = X_year,
      family = family_id,
      time_varying = as.numeric(time_varying),
      estimate_self_report = self_report_data$estimate_self_report,
      n_vessels_observed = self_report_data$n_vessels_observed,
      n_vessels_self_report = self_report_data$n_vessels_self_report
    )

    # Store stream info
    stream_info <- list(
      "takes_em" = takes_em,
      "effort_em" = effort_em,
      "covrate_em" = covrate_em,
      "takes_both" = takes_both,
      "effort_both" = effort_both,
      "covrate_obs" = covrate_obs,
      "covrate_both" = covrate_both,
      "effort_total" = effort_total,
      "n_obs" = n_obs,
      "n_em" = n_em,
      "n_both" = n_both,
      "takes_selfreport" = takes_selfreport,
      "effort_selfreport" = effort_selfreport,
      "n_selfreport" = n_selfreport,
      "covrate_selfreport" = covrate_selfreport,
      "vessels_observed" = vessels_observed,
      "estimate_self_report" = as.logical(self_report_data$estimate_self_report)
    )

  } else {
    # Single stream mode

    observed_effort <- data[[effort]]

    # Calculate unobserved effort using coverage rate
    if (!is.null(covrate_obs)) {
      # Use coverage rate
      if (!(covrate_obs %in% colnames(data))) {
        stop(paste("Column", covrate_obs, "not found in data"))
      }

      p_coverage <- data[[covrate_obs]] / 100

      # Calculate total and unobserved effort from coverage
      total_effort <- observed_effort / p_coverage
      new_effort <- total_effort - observed_effort

    } else if (!is.null(effort_total)) {
      # Use effort_total directly
      total_effort <- data[[effort_total]]
      new_effort <- total_effort - observed_effort

      # Check for negatives
      if (any(new_effort < 0, na.rm = TRUE)) {
        cli_warn("Some observations have effort exceeding effort_total. Setting unobserved effort to 0.")
        new_effort[new_effort < 0] <- 0
      }
    } else {
      # No expansion (100% coverage)
      total_effort <- observed_effort
      new_effort <- rep(0, nrow(data))
      cli_inform("No expansion information provided - assuming 100% coverage")
    }

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

    # Create year-level covariate matrix
    X_year <- matrix(nrow = n_year, ncol = ncol(X))
    for(k in 1:ncol(X)) {
      for(t in 1:n_year) {
        idx <- which(time_idx == t)[1]
        X_year[t, k] <- X[idx, k]
      }
    }

    datalist <- list(
      multi_stream = 0,
      n_obs = 0,
      n_em = 0,
      n_both = 0,
      n_selfreport = 0,
      yint_obs = integer(0),
      yint_em = integer(0),
      yint_both = integer(0),
      yint_selfreport = integer(0),
      yreal_obs = numeric(0),
      yreal_em = numeric(0),
      yreal_both = numeric(0),
      yreal_selfreport = numeric(0),
      effort_obs = numeric(0),
      effort_em = numeric(0),
      effort_both = numeric(0),
      effort_selfreport = numeric(0),
      time_idx_obs = integer(0),
      time_idx_em = integer(0),
      time_idx_both = integer(0),
      time_idx_selfreport = integer(0),
      new_effort_by_year = new_effort_by_year,
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
      time_varying = as.numeric(time_varying),
      # Self-reporting is only supported alongside the multi-stream self-report
      # stream (see validation above), so this is always off in single-stream mode.
      estimate_self_report = 0L,
      n_vessels_observed = rep(0L, n_year),
      n_vessels_self_report = rep(0L, n_year)
    )

    stream_info <- NULL
  }

  # Define parameters to extract
  if(multi_stream) {
    pars <- c("beta", "lambda_base", "log_lambda_base", "log_lik", "y_new",
              "est_time_dev", "sigma_rw", "sigma_logn", "cv_gamma",
              "nb2_phi", "theta", "y_new_real", "logit_p_report",
              "p_report", "new_effort_adjusted", "effort_selfreport_new_coverage")
  } else {
    pars <- c("beta", "lambda", "log_lambda", "log_lik", "y_new",
              "est_time_dev", "sigma_rw", "sigma_logn", "cv_gamma",
              "nb2_phi", "theta", "y_new_real", "logit_p_report",
              "p_report", "new_effort_adjusted", "effort_selfreport_new_coverage")
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
    "stream_info" = stream_info,
    "estimate_self_report" = as.logical(datalist$estimate_self_report)
  )

  return(result)
}


#' Extract log-likelihood matrix for loo, removing NA observations
#'
#' This function extracts the log-likelihood array from a fitted bycatch model
#' and removes any columns (observations) that have NA values in any MCMC iteration.
#' This is necessary because loo::loo() cannot handle NA values.
#'
#' @param fit A bycatch model fit object (output from fit_bycatch or fit_bycatch_modified)
#' @return A matrix suitable for loo::loo() with dimensions (n_iterations x n_valid_observations)
#' @export
#'
#' @importFrom rstan extract
#'
#' @examples
#' \donttest{
#' # Create example data
#' d <- data.frame(
#'   Year = 2002:2014,
#'   Takes = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
#'   CovRate = c(24, 22, 14, 32, 28, 25, 30, 7, 26, 21, 22, 23, 27),
#'   Sets = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486)
#' )
#'
#' # Fit a model
#' fit <- fit_bycatch(Takes ~ 1,
#'   data = d,
#'   time = "Year",
#'   effort = "Sets",
#'   covrate = "CovRate",
#'   family = "poisson",
#'   time_varying = FALSE
#' )
#'
#' # Extract log-likelihood for LOO-CV
#' log_lik <- extract_log_lik_for_loo(fit)
#'
#' # Run LOO-CV
#' library(loo)
#' loo_result <- loo(log_lik)
#' }
extract_log_lik_for_loo <- function(fit) {
  # Extract log-likelihood array from Stan fit
  # Dimensions: (n_iterations x n_observations)
  log_lik_array <- rstan::extract(fit$fitted_model, pars = "log_lik")$log_lik

  # Check for NAs in any iteration for each observation
  na_cols <- apply(log_lik_array, 2, function(x) any(is.na(x)))

  # Remove observations with NA values and warn user
  if (any(na_cols)) {
    warning(paste("Removing", sum(na_cols), "observations with NA log-likelihood values"))
    log_lik_array <- log_lik_array[, !na_cols, drop = FALSE]
  }

  return(log_lik_array)
}
