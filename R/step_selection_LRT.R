#' Likelihood ratio test for step-specific selection
#'
#' Tests, for each variant (or compound variant), whether a step-specific model of selection (see
#' \code{ces_variant_step()}) fits significantly better than the constrained model in which one
#' selection intensity is shared across all stages.
#'
#' By default, the constrained model is the one fit alongside the step-specific model by
#' \code{ces_variant_step()} (\code{loglikelihood_constant} column). It uses the same likelihood
#' and gene mutation rates, so the two models are nested by construction. Alternatively, supply
#' \code{simple_run_name} to compare against a separate \code{ces_variant()} run (with
#' \code{return_fit = TRUE}) over the same variants. Only do this if that run's model is nested in
#' the step-specific model: each sample's gene mutation rate in the \code{ces_variant()} run must
#' equal the cumulative rate for that sample's own stage.
#'
#' The step-specific model has one additional free parameter per extra stage beyond the first
#' (that is, \code{num_stages - 1} more parameters than the constrained model's single selection
#' intensity), which sets the test's degrees of freedom.
#'
#' @param cesa CESAnalysis object containing the selection_results run(s).
#' @param step_run_name run_name of a previous \code{ces_variant_step()} call.
#' @param simple_run_name Optional run_name of a previous \code{ces_variant()} call over the same
#'   variants, to use instead of the built-in constrained fit (see details).
#' @param max_loglik Log-likelihoods above this value (default 1e5) are treated as implausible
#'   (usually a sign of optimizer non-convergence) and set to NA, along with the LRT statistic and
#'   p-value for that variant.
#' @return A data.table with columns variant_name, loglik_step, loglik_simple, df, LRT_stat, and
#'   p_value.
#' @export
step_selection_LRT = function(cesa = NULL, step_run_name = NULL, simple_run_name = NULL, max_loglik = 1e5) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (! is.character(step_run_name) || length(step_run_name) != 1) {
    stop("step_run_name should be 1-length character.")
  }
  if (! is.null(simple_run_name) && (! is.character(simple_run_name) || length(simple_run_name) != 1)) {
    stop("simple_run_name should be 1-length character, or NULL to use the constrained fit from ces_variant_step().")
  }
  if (! is.numeric(max_loglik) || length(max_loglik) != 1) {
    stop("max_loglik should be 1-length numeric.")
  }
  missing_runs = setdiff(c(step_run_name, simple_run_name), names(cesa@selection_results))
  if (length(missing_runs) > 0) {
    stop("Run name(s) not found in cesa's selection results: ", paste(missing_runs, collapse = ", "), ".")
  }

  step_results = cesa@selection_results[[step_run_name]]
  si_cols = grep("^si_", names(step_results), value = TRUE)
  if (length(si_cols) < 2 || ! "loglikelihood" %in% names(step_results)) {
    stop("step_run_name doesn't look like ces_variant_step() output (no si_<stage> columns found).")
  }
  df = length(si_cols) - 1

  if (is.null(simple_run_name)) {
    if (! "loglikelihood_constant" %in% names(step_results)) {
      stop("step_run_name has no loglikelihood_constant column; supply simple_run_name instead.")
    }
    lrt_table = step_results[, .(variant_name, loglik_step = loglikelihood, loglik_simple = loglikelihood_constant)]
  } else {
    simple_results = cesa@selection_results[[simple_run_name]]
    if (! "loglikelihood" %in% names(simple_results)) {
      stop("simple_run_name has no loglikelihood column (make sure it was run with return_fit = TRUE).")
    }
    shared_variants = intersect(step_results$variant_name, simple_results$variant_name)
    if (length(shared_variants) == 0) {
      stop("No shared variant_name values between step_run_name and simple_run_name results.")
    }
    if (length(shared_variants) < step_results[, .N]) {
      pretty_message(paste0("Note: ", step_results[, .N] - length(shared_variants),
                            " variant(s) in step_run_name have no match in simple_run_name and will be excluded."))
    }
    step_ll = step_results[shared_variants, .(variant_name, loglik_step = loglikelihood), on = "variant_name"]
    simple_ll = simple_results[shared_variants, .(variant_name, loglik_simple = loglikelihood), on = "variant_name"]
    lrt_table = step_ll[simple_ll, on = "variant_name"]
  }

  # Guard against implausible loglikelihoods from failed optimizer convergence
  lrt_table[loglik_step > max_loglik, loglik_step := NA_real_]
  lrt_table[loglik_simple > max_loglik, loglik_simple := NA_real_]

  lrt_table[, df := df]
  # the step model nests the simple model, so a (tiny) negative statistic is optimizer tolerance
  lrt_table[, LRT_stat := pmax(-2 * (loglik_simple - loglik_step), 0)]
  lrt_table[, p_value := stats::pchisq(LRT_stat, df = df, lower.tail = FALSE)]

  return(lrt_table[])
}
