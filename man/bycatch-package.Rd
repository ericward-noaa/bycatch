\name{bycatch-package}
\alias{bycatch-package}
\alias{bycatch}
\docType{package}
\title{
\packageTitle{bycatch}
}
\description{
\packageDescription{bycatch}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bycatch}
\packageIndices{bycatch}
Following Martin et al. (2015), our models use a rate parameter \eqn{\theta} representing the per set (or hook, or other metric observed) take rate. The take rate is combined with measures of effort \eqn{E_{t}} (observed sets) to calculate the expected number of observed takes, \deqn{{ \lambda  }_{ i }=\theta E_{ t }}. The observed takes are related to expected by assuming a Poisson distribution, \eqn{Y_{t} ~ Poisson(\lambda_{t})}. Covariates can be included in this framework in several ways. First, we can include a covariate matrix \eqn{X} and estimate the effects of each covariate on \eqn{\theta} using a log-link, \deqn{log(\theta_{t}) = B_{0} + B_{1}*X_{t,1} + B_{2}*X_{t,2}}. A second way to include covariates is as a random walk in time (analogous to a DLM, using the `time_varying=TRUE` argument). In this setting, each year has an estimated temporal deviation modeled as a random walk,
\deqn{{ \delta  }_{ t }~ Normal\left( { \delta  }_{ t-1 },\sigma  \right)} and these are incorporated into the bycatch rate as \deqn{log(\theta_{t}) = B_{0} + B_{1}*X_{t,1} + \delta_{t}...}.

This framework multiplies the estimated per set bycatch rate \eqn{\theta_{t}} by the known effort \eqn{E_{t}} to generate the mean bycatch in each year. This expectation represents the mean, as if this were an experiment, and multiple samples were generated repeatedly. In scenarios where observer coverage is 100 percent, no expansion is necessary, and these estimates will be used as the total expected bycatch. Because intervals around these estimates represent intervals around a mean, we don't expect the observed data to be contained 95 percent of the time (or any other threshold). *If we included an observation model on top of these estimates, we would generate intervals with these properties*

In situations where observer coverage is \eqn{0 < p < 1}, we need to extrapolate these estimates to the unobserved sets. One approach for doing this is to view each estimated number of takes (corrected for effort) \eqn{\lambda_{t}} as the \eqn{x} of a Binomial distribution with known \eqn{p} (observer coverage is known) but unknown N. This algorithm works by iterating over MCMC draws (\eqn{i = 1, 2, ...,n}). For each draw and in each year, we can generate \eqn{{ N }_{ i,t }|p,{ x }_{ i,t }={ \lambda  }_{ i,t }}. because the relative probabilities of each \eqn{{ N }_{ i,t }} can be computed using `dbinom(x = x, size = N, prob = p)`.

In terms of interpretation, the only constraints on \eqn{{ N }_{ i,t }} are that they have to be greater than \eqn{{ x }_{ i,t }}. Like the \eqn{{ x }_{ i,t }}, the \eqn{{ N }_{ i,t }} represents the expanded mean (and thus, the 95 percent CIs are not expected to contain most of the observations).
}
\author{
\packageAuthor{bycatch}

Maintainer: \packageMaintainer{bycatch}
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
Optional links to other man pages, e.g.
}
\examples{
#simple examples of the most important functions
}
