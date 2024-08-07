% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/si_uniroot_conf_int.R
\name{univariate_si_conf_ints}
\alias{univariate_si_conf_ints}
\title{Calculate uniroot CIs on selection intensities}
\usage{
univariate_si_conf_ints(fit, lik_fn, min_si, max_si, conf)
}
\arguments{
\item{fit}{From bbmle}

\item{lik_fn}{likelihood function}

\item{min_si}{lower limit on SI/CI}

\item{max_si}{upper limit on SI/CI}

\item{conf}{e.g., .95 -> 95\% CIs}
}
\description{
Given a model fit, calculate univariate confidence intervals for each parameter.
Returns a list of low/high bounds.
}
\keyword{internal}
