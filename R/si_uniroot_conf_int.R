#' Calculate uniroot CIs on selection intensities
#' 
#' Given a model fit, calculate univariate confidence intervals for each parameter.
#' Returns a list of low/high bounds.
#' 
#' @param fit From bbmle
#' @param lik_fn likelihood function
#' @param min_si lower limit on SI/CI
#' @param max_si upper limit on SI/CI
#' @param conf e.g., .95 -> 95\% CIs
#' @keywords internal
univariate_si_conf_ints = function(fit, lik_fn, min_si, max_si, conf) {
  max_ll = -1 * as.numeric(bbmle::logLik(fit))
  offset = stats::qchisq(conf, 1)/2
  selection_intensity = bbmle::coef(fit)
  num_pars = length(selection_intensity)
  conf_ints = list(num_pars * 2)
  for (i in 1:num_pars) {
    if(is.na(selection_intensity[i])) {
      lower = NA_real_
      upper = NA_real_
    } else {
      # univariate likelihood function freezes all but one SI at MLE
      # offset makes output at MLE negative; function should be positive at lower/upper boundaries,
      # and uniroot will find the zeroes, which should represent the lower/uppper CIs
      ulik = function(x) { 
        pars = selection_intensity
        pars[i] = x
        return(lik_fn(pars) - max_ll - offset)
      }
      # if ulik() of the floor SI is negative, no root on [floor MLE], so setting an NA lower bound
      if(ulik(min_si) < 0) {
        lower = NA_real_
      } else {
        lower = max(stats::uniroot(ulik, lower = min_si, upper = selection_intensity[i])$root, min_si)
      }
      if(ulik(max_si) < 0){
        # this really shouldn't happen
        upper = NA_real_
      } else {
        upper = stats::uniroot(ulik, lower = selection_intensity[i], upper = max_si)$root
      }
    }
    conf_ints[[i * 2 - 1]] = lower
    conf_ints[[i * 2]] = upper
  }
  si_names = bbmle::parnames(lik_fn)
  ci_high_colname = paste0("ci_high_", conf * 100)
  ci_low_colname = paste0("ci_low_", conf * 100)
  if (num_pars == 1) {
    ci_colnames = c(ci_low_colname, ci_high_colname)
  } else {
    low_colnames = paste(ci_low_colname, si_names, sep = "_")
    high_colnames = paste(ci_high_colname, si_names, sep = "_")
    ci_colnames = unlist(S4Vectors::zipup(low_colnames, high_colnames))
  }
  names(conf_ints) = ci_colnames
  return(conf_ints)
}