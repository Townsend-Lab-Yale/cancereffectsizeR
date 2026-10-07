#' Linear continuous-covariate SSWM likelihood
#'
#' Constructs a negative log-likelihood function for site-level selection under
#' the strong-selection, weak-mutation (SSWM) model, with selection intensity
#' modeled as a linear function of a continuous covariate:
#'
#' \deqn{\gamma(x) = \beta_0 + \beta_1 x,}
#'
#' where \eqn{x} is the sample-level continuous covariate, \eqn{\beta_0} is
#' the intercept, and \eqn{\beta_1} is the covariate-dependent slope.
#'
#' This function is a likelihood-function factory intended for use with
#' \code{ces_variant_linear()}. For each variant, the mutation rates and sample
#' information required by the likelihood are supplied automatically by the
#' cancer effect size inference workflow.
#'
#' The returned function accepts a parameter vector
#' \code{c(beta0, beta1)} and returns the negative log-likelihood.
#'
#' @param rates_tumors_with vector of site-specific mutation rates for all
#' tumors with variant
#'
#' @param rates_tumors_without vector of site-specific mutation rates for all
#' eligible tumors without variant
#'
#' @param covariate_data A \code{data.table} containing
#'   \code{Unique_Patient_Identifier} and \code{covariate_value}. The
#'   \code{covariate_value} column contains the sample-level continuous
#'   covariate and must be coercible to numeric.
#'
#' @return A function of a numeric parameter vector
#'   \code{c(beta0, beta1)} that returns the negative log-likelihood for the
#'   variant under the linear continuous-covariate SSWM model.
#'
#' @seealso \code{\link{ces_variant_linear}},
#'   \code{\link{sswm_age_lik_logistic}},
#'   \code{\link{sswm_age_lik_sigmoid}}
#'
#' @family continuous selection likelihoods
#'
#' @export
sswm_age_lik <- function(rates_tumors_with, rates_tumors_without, covariate_data, selection_results=NULL) {
  fn <- function(params) {
    beta0 <- params[1]
    beta1 <- params[2]

    # Retrieve covariate values for the tumors with and without the variant
    covariate_tumors_with <- covariate_data[names(rates_tumors_with), covariate_value, on = "Unique_Patient_Identifier"]
    covariate_tumors_without <- covariate_data[names(rates_tumors_without), covariate_value, on = "Unique_Patient_Identifier"]
    covariate_tumors_with <- as.numeric(covariate_tumors_with)
    covariate_tumors_without <- as.numeric(covariate_tumors_without)

    gamma_tumors_with <- beta0 + beta1 * covariate_tumors_with
    gamma_tumors_without <- beta0 + beta1 * covariate_tumors_without

    # ISRES may evaluate parameter values outside the nonlinear constraint
    # region. Reject invalid selection intensities before evaluating the
    # likelihood.
    if (any(!is.finite(gamma_tumors_with)) ||
        any(!is.finite(gamma_tumors_without)) ||
        any(gamma_tumors_with <= 0) ||
        any(gamma_tumors_without <= 0)) {
      return(1e100)
    }

    z_tumors_with <- gamma_tumors_with * rates_tumors_with
    z_tumors_without <- gamma_tumors_without * rates_tumors_without

    if (any(!is.finite(z_tumors_with)) ||
        any(!is.finite(z_tumors_without)) ||
        any(z_tumors_with <= 0) ||
        any(z_tumors_without < 0)) {
      return(1e100)
    }

    sum_log_lik <- 0

    # Calculate likelihood for tumors without the variant
    if (length(rates_tumors_without) > 0) {
      sum_log_lik <- sum_log_lik - sum(z_tumors_without)
    }

    # Calculate likelihood for tumors with the variant
    if (length(rates_tumors_with) > 0) {
      sum_log_lik <- sum_log_lik + sum(log(-expm1(-z_tumors_with)))
    }

    # Return negative log-likelihood
    return(-1 * sum_log_lik)
  }
  # Set default values for beta0 and beta1
  formals(fn)[["params"]] <- c(beta0 = 1, beta1 = 0)
  bbmle::parnames(fn) <- c("beta0", "beta1")
  return(fn)
}

#' Logistic continuous-covariate SSWM likelihood
#'
#' Constructs a negative log-likelihood function for site-level selection under
#' the strong-selection, weak-mutation (SSWM) model, with selection intensity
#' modeled as a logistic function of a continuous covariate:
#'
#' \deqn{
#' \gamma(x) = \frac{\exp(L)}{1 + \exp[-k(x - m)]}.
#' }
#'
#' Here, \eqn{\exp(L)} determines the upper asymptote, \eqn{k} determines the
#' steepness or growth/decay rate of the transition, and \eqn{m} is the
#' midpoint of the transition. A positive \eqn{k} gives an increasing curve and
#' a negative \eqn{k} gives a decreasing curve.
#'
#' This function is a likelihood-function factory intended for use with
#' \code{ces_variant_logistic()}. For each variant, baseline mutation rates and
#' sample information are supplied automatically by the cancer effect size
#' inference workflow.
#'
#'
#' @param rates_tumors_with vector of site-specific mutation rates for all
#' tumors with variant
#'
#' @param rates_tumors_without vector of site-specific mutation rates for all
#' eligible tumors without variant
#'
#' @param covariate_data A \code{data.table} containing
#'   \code{Unique_Patient_Identifier} and \code{covariate_value}. The
#'   \code{covariate_value} column contains the sample-level continuous
#'   covariate and must be coercible to numeric.
#'
#' @return A function of a three-element numeric parameter vector corresponding
#'   to \code{L}, \code{k}, and \code{m}. The returned function evaluates
#'   the negative log-likelihood for the variant under the logistic
#'   continuous-covariate SSWM model.
#'
#' @seealso \code{\link{ces_variant_logistic}},
#'   \code{\link{sswm_age_lik}},
#'   \code{\link{sswm_age_lik_sigmoid}}
#'
#' @family continuous selection likelihoods
#'
#' @export
sswm_age_lik_logistic <- function(rates_tumors_with, rates_tumors_without, covariate_data, selection_results=NULL) {
  fn <- function(params) {
    L <- params[1]
    k <- params[2]
    m <- params[3]

    # Retrieve covariate values for the tumors with and without the variant
    covariate_tumors_with <- covariate_data[names(rates_tumors_with), covariate_value, on = "Unique_Patient_Identifier"]
    covariate_tumors_without <- covariate_data[names(rates_tumors_without), covariate_value, on = "Unique_Patient_Identifier"]
    covariate_tumors_with <- as.numeric(covariate_tumors_with)
    covariate_tumors_without <- as.numeric(covariate_tumors_without)

    sum_log_lik <- 0

    # Calculate likelihood for tumors without the variant
    if (length(rates_tumors_without) > 0) {
      sum_log_lik <- sum_log_lik - sum(mapply(function(rate, covariate_value) {
        gamma_covariate <- exp(L) / (1 + exp(-k * (covariate_value - m)))
        return(gamma_covariate * rate)
      }, rates_tumors_without, covariate_tumors_without))
    }

    # Calculate likelihood for tumors with the variant
    if (length(rates_tumors_with) > 0) {
      sum_log_lik <- sum_log_lik + sum(mapply(function(rate, covariate_value) {
        gamma_covariate <- exp(L) / (1 + exp(-k * (covariate_value - m)))
        return(log(1 - exp(-1 * gamma_covariate * rate)))
      }, rates_tumors_with, covariate_tumors_with))
    }

    # Return negative log-likelihood
    return(-1 * sum_log_lik)
  }
  # Set default values for L, k, m
  formals(fn)[["params"]] <- c(L = 1, k = 0, m = median(as.numeric(covariate_data$covariate_value)))
  bbmle::parnames(fn) <- c("L", "k", "m")
  return(fn)
}



#' Generalized sigmoid continuous-covariate SSWM likelihood
#'
#' Constructs a negative log-likelihood function for site-level selection under
#' the strong-selection, weak-mutation (SSWM) model, with selection intensity
#' modeled as a four-parameter generalized sigmoid function of a continuous covariate:
#'
#' \deqn{
#' \gamma(x) =
#' \exp(C) + (\exp(L)-\exp(C))\frac{x^s}{x^s + m^s}.
#' }
#'
#' The parameters \eqn{\exp(C)} and \eqn{\exp(L)} define the lower and upper
#' asymptotic selection intensities, \eqn{m} specifies the covariate value at
#' which selection intensity is halfway between the two asymptotes, and \eqn{s}
#' is the shape exponent controlling how sharply the curve transitions around
#' \eqn{m}.
#'
#' This function is a likelihood-function factory intended for use with
#' \code{ces_variant_sigmoid()}. For each variant, baseline mutation rates and
#' sample information are supplied automatically by the cancer effect size
#' inference workflow.
#'
#' @param rates_tumors_with vector of site-specific mutation rates for all
#' tumors with variant
#'
#' @param rates_tumors_without vector of site-specific mutation rates for all
#' eligible tumors without variant
#'
#' @param covariate_data A \code{data.table} containing
#'   \code{Unique_Patient_Identifier} and \code{covariate_value}. The
#'   \code{covariate_value} column contains the sample-level continuous
#'   covariate and must be coercible to numeric. For the generalized sigmoid
#'   model, covariate values must be strictly positive.
#'
#' @return A function of a four-element numeric parameter vector corresponding
#'   to \code{C}, \code{L}, \code{s}, and \code{m}. The returned function
#'   evaluates the negative log-likelihood for the variant under the generalized sigmoid
#'   continuous-covariate SSWM model.
#'
#' @seealso \code{\link{ces_variant_sigmoid}},
#'   \code{\link{sswm_age_lik}},
#'   \code{\link{sswm_age_lik_logistic}}
#'
#' @family continuous selection likelihoods
#'
#' @export
sswm_age_lik_sigmoid <- function(rates_tumors_with, rates_tumors_without, covariate_data, selection_results=NULL) {
  fn <- function(params) {
    C <- params[1]
    L <- params[2]
    s <- params[3]
    m <- params[4]

    # Retrieve covariate values for the tumors with and without the variant
    covariate_tumors_with <- covariate_data[names(rates_tumors_with), covariate_value, on = "Unique_Patient_Identifier"]
    covariate_tumors_without <- covariate_data[names(rates_tumors_without), covariate_value, on = "Unique_Patient_Identifier"]
    covariate_tumors_with <- as.numeric(covariate_tumors_with)
    covariate_tumors_without <- as.numeric(covariate_tumors_without)

    # The generalized sigmoid model requires positive covariate values and a
    # positive midpoint m.
    if (m <= 0 ||
        any(!is.finite(covariate_tumors_with)) ||
        any(!is.finite(covariate_tumors_without)) ||
        any(covariate_tumors_with <= 0) ||
        any(covariate_tumors_without <= 0)) {
      return(1e100)
    }

    # Numerically stable form of x^s / (x^s + m^s).
    sigmoid_weight_with <- plogis(s * (log(covariate_tumors_with) - log(m)))
    sigmoid_weight_without <- plogis(s * (log(covariate_tumors_without) - log(m)))

    gamma_tumors_with <- exp(C) + (exp(L) - exp(C)) * sigmoid_weight_with
    gamma_tumors_without <- exp(C) + (exp(L) - exp(C)) * sigmoid_weight_without

    if (any(!is.finite(gamma_tumors_with)) ||
        any(!is.finite(gamma_tumors_without)) ||
        any(gamma_tumors_with <= 0) ||
        any(gamma_tumors_without <= 0)) {
      return(1e100)
    }

    z_tumors_with <- gamma_tumors_with * rates_tumors_with
    z_tumors_without <- gamma_tumors_without * rates_tumors_without

    if (any(!is.finite(z_tumors_with)) ||
        any(!is.finite(z_tumors_without)) ||
        any(z_tumors_with <= 0) ||
        any(z_tumors_without < 0)) {
      return(1e100)
    }

    sum_log_lik <- 0

    # Calculate likelihood for tumors without the variant
    if (length(rates_tumors_without) > 0) {
      sum_log_lik <- sum_log_lik - sum(z_tumors_without)
    }

    # Calculate likelihood for tumors with the variant
    if (length(rates_tumors_with) > 0) {
      sum_log_lik <- sum_log_lik + sum(log(-expm1(-z_tumors_with)))
    }

    # Return negative log-likelihood
    return(-1 * sum_log_lik)
  }
  # Set default values for C, L, s, m
  formals(fn)[["params"]] <- c(C = 1, L = 2, s = 1, m = median(as.numeric(covariate_data$covariate_value)))
  bbmle::parnames(fn) <- c("C", "L", "s", "m")
  return(fn)
}
