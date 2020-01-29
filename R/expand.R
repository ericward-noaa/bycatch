#' expand does the expansion from a fitted object
#' @param fit the estimated bycatch rate from a call to fit
#' @param coverage the coverage rate (0-100). This can be a scalar if all values are the same, or a vector corresponding to each row of the data frame used in fitting
#' @param control list of 3 elements used for control: sigfig_multiplier (used for adjusting precision of estimates, defaults to 1000), mcmc_samples (number of mcmc samples randomly selected from the posterior for expansion), maxX (upper bound for drawing from latent bycatch events, also related to sigfig_multiplier)
#'
#' @return expanded_estimates the estimates of expanded estimates across rows of the dataframe and mcmc samples
#'
#' @export
#' @importFrom rstan extract
#' @importFrom stats dbinom rpois
#'
#' @examples
#' \donttest{
#' d = data.frame("Year"= 2002:2014,
#' "Takes" = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0),
#' "expansionRate" = c(24, 22, 14, 32, 28, 25, 30,  7, 26, 21, 22, 23, 27),
#' "Sets" = c(391, 340, 330, 660, 470, 500, 330, 287, 756, 673, 532, 351, 486))
#' fit = fit_bycatch(Takes ~ 1, data=d,
#' time="Year", effort="Sets", family="poisson", time_varying = FALSE)
#' expanded = expand(fit, coverage = d$expansionRate)
#' }
expand <- function(fit, coverage = NULL,
  control = list(sigfig_multiplier = 100, mcmc_samples = 1000, maxX = 20000)) {

  if(is.null(coverage)) {
    stop("Please remember to enter a vector of numbers (0-100) that is a scalar representing constant coverage or a vector representing coverage for each row of the data frame that was used to fit the model")
  }

  if(max(coverage) < 1) coverage = coverage * 100
  fit$data$coverage = coverage

  # extract parameters
  pars = rstan::extract(fit$fitted_model)
  fit$data$lambda = apply(pars$lambda, 2, mean)

  # do binomial expansion of estimated takes
  # estimated takes can be thought of as the observed draw from binomial distribution where p (but not N) is known
  binom_p = fit$data$coverage/100
  # point estimates may be unstable - see Raftery (1988) - so we can simulate distribution

  # sample ~ 1000 random draws from posterior
  sigfig_multiplier = control$sigfig_multiplier # default to 100
  mcmc_samples = control$mcmc_samples # default to 1000
  maxX = control$maxX # default to 20000

  samples = sample(1:dim(pars$lambda)[1], size=mcmc_samples, replace=F)
  expanded_estimates = matrix(0, length(samples), nrow(fit$data))

  for(y in 1:nrow(fit$data)) {
    for(mcmc in 1:length(samples)) {
      # We observe takes, df$events[y]
      # pmf of N | p
      # In other words, we need to sample from the density of Binomial N
      # given Binomial p and X

      # X needs to be an integer, so expand to a series of large numbers as approximation
      X = round(pars$lambda[samples[mcmc],y] * sigfig_multiplier)
      # calculate the probabilites of these Ns given the mean observed takes and observer coverage
      prob_N = dbinom(x=X, size = seq(X, maxX), prob = binom_p[y])
      if(X/binom_p[y] < maxX) {
        warning(
          paste("Warning: maxX in the control list needs to be increased",
            " to at least ", ceiling(X/binom_p[y]),". Or try smaller values of ",
            "the sigfig_multiplier parameter. Please try running the expand function again.")
        )
      }
      # sample from distribution of N, representing total mean takes
      mean_takes = try(sample(seq(X, maxX),size=1, prob=prob_N),
        silent=TRUE)
      if(class(mean_takes)=="try-error") {
        warning(paste("Warning: maxX in the control list needs to be increased",
          " to at least ", ceiling(X/binom_p[y]),". Or try smaller values of ",
          "the sigfig_multiplier parameter. Please try running the expand function again.")
        ))
      }
      # sample from poisson to convert mean -> observed data with observation model
      if(fit$family == "poisson") {
        N = rpois(n = 1, lambda = mean_takes)
      } else {
        N = MASS::rnegbin(n = 1, mu = mean_takes, theta = pars$nb2_phi[samples[mcmc]])
      }
      # use the expansion for unobserved sets, use observed as known perfectly for observed
      N = (1-binom_p[y]) * N / sigfig_multiplier + fit$data[y,fit$events]
      expanded_estimates[mcmc,y] = N
    }
  }

  return(expanded_estimates)
}
