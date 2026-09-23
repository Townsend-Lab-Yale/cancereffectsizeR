#' step_selection_lik
#'
#' Generates a log-likelihood function for step-specific selection: a variant's scaled selection
#' coefficient is allowed to differ across an ordered sequence of tumor progression stages (for
#' example, normal tissue and primary tumor), rather than assuming one constant coefficient across
#' all samples (compare \code{sswm_lik()}). All arguments to this likelihood function factory are
#' automatically supplied by \code{ces_variant_step()}.
#'
#' @param rates_tumors_with vector of site-specific mutation rates for all tumors with variant
#' @param rates_tumors_without vector of site-specific mutation rates for all eligible tumors
#'   without variant
#' @param sample_index data.table with columns Unique_Patient_Identifier, group_index, and
#'   group_name, associating each sample with its (1-indexed) stage and the stage's display name.
#'   See \code{assign_stage_index()}.
#' @param stage_mut_prop Numeric vector, in stage order, giving the proportion of the variant's
#'   gene-level mutation rate estimated to accumulate during each stage. Should sum to 1. See
#'   \code{stage_mutation_proportions()}.
#' @export
step_selection_lik = function(rates_tumors_with, rates_tumors_without, sample_index, stage_mut_prop) {
  stage_mut_prop = as.numeric(stage_mut_prop)
  num_pars = sample_index[, uniqueN(group_index)]
  if (length(stage_mut_prop) != num_pars) {
    stop("stage_mut_prop has ", length(stage_mut_prop), " value(s) but sample_index has ",
         num_pars, " stage(s); these must match.")
  }

  stages_tumors_with = sample_index[names(rates_tumors_with), group_index, on = "Unique_Patient_Identifier"]
  stages_tumors_without = sample_index[names(rates_tumors_without), group_index, on = "Unique_Patient_Identifier"]

  fn = function(gamma) {
    gamma = unname(gamma) # math faster on unnamed vectors

    # samples without the variant: expected mutational flux across every stage they've reached
    sum_log_lik = -1 * sum(mapply(
      function(rate, stage) {
        flux = gamma[1:stage] * rate * stage_mut_prop[1:stage]
        return(sum(flux))
      }, rates_tumors_without, stages_tumors_without))

    # samples with the variant: likelihood of first acquiring it during the stage it was observed,
    # having survived without it through every earlier stage
    if (length(rates_tumors_with) > 0) {
      sum_log_lik = sum_log_lik + sum(mapply(
        function(rate, stage) {
          lik_no_mutation = exp(-1 * gamma * rate * stage_mut_prop)
          lik_mutation = 1 - lik_no_mutation
          cum_lik_no_mut = c(1, cumprod(lik_no_mutation))
          return(log(sum(cum_lik_no_mut[1:stage] * lik_mutation[1:stage])))
        }, rates_tumors_with, stages_tumors_with))
    }

    # guards against the optimizer trying all-boundary values at once
    if (! is.finite(sum_log_lik)) {
      return(-1e200)
    }
    return(-1 * sum_log_lik) # negative log-likelihood, as bbmle::mle2 expects
  }

  # Set default starting values for all stage SIs, which ces_variant_step() uses to initialize the optimizer
  formals(fn)[["gamma"]] = rep.int(1000, num_pars)

  # bbmle::mle2 requires named parameters when vecpar = TRUE
  group_names = unique(sample_index[, .(group_index, group_name)], by = "group_index")[order(group_index), group_name]
  bbmle::parnames(fn) = paste0("si_", group_names)
  return(fn)
}
