#' Estimate selection under pairwise epistasis model
#' 
#' @param cesa CESAnalysis
#' @param samples Validated samples data.table (as from select_samples())
#' @param variant_pair 2-length character of variant IDs, or 2-length numeric giving
#'   indices of VariantSetList for the current two compound variants
#' @param compound_variants If testing a pair of compound variants, the VariantSetList defining them
#' @param model Passed from ces_epistasis or ces_gene_epistasis. Set to "default" to use built-in
#'   model of epistatic selection, or supply a custom function factory (see details).
#' @param lik_args Extra arguments, given as a list, passed from ces_epistasis or ces_gene_epistasis
#' @param optimizer_args List of arguments to pass to the optimizer (bbmle::mle2).
#' @param pval_calc_fn For use with custom models; optional. A function that takes an epistasis
#'   model fit as input and returns p-values and other descriptives.
#' 
#' @keywords internal
pairwise_variant_epistasis = function(cesa, variant_pair, samples, conf, compound_variants = NULL, model = "default", 
                                      lik_args = list(), optimizer_args = list(), pval_calc_fn = NULL) {
  
  running_compound = FALSE
  if (is(compound_variants, "VariantSetList")) {
    running_compound = TRUE
    compound_variants = compound_variants[variant_pair]
    joint_coverage = c("genome", intersect(compound_variants@compounds$shared_cov[[1]], compound_variants@compounds$shared_cov[[2]]))
    v1 = compound_variants@compounds$set_id[[1]]
    v2 = compound_variants@compounds$set_id[[2]]
    v1_ids = compound_variants@sbs[set_id == v1, sbs_id]
    v2_ids = compound_variants@sbs[set_id == v2, sbs_id]
    variant_ids = c(v1_ids, v2_ids)
  } else {
    v1 = variant_pair[1]
    v2 = variant_pair[2]
    
    v1_coverage = cesa@mutations$variants_to_cov[[v1]]
    v2_coverage = cesa@mutations$variants_to_cov[[v2]]
    
    # Samples have to have v1 and v2 coverage (and samples with covered_regions == "genome" always have coverage)
    joint_coverage = c("genome", intersect(v1_coverage, v2_coverage))
    variant_ids = c(v1, v2)
  }
  
  eligible_tumors = samples[covered_regions %in% joint_coverage, patient_id]
  if (length(eligible_tumors) == 0) {
    warning(sprintf("No samples have coverage at both %s and %s, so this variant pair had to be skipped.", v1, v2), immediate. = T, call. = F)
  }
  all_rates = baseline_mutation_rates(cesa = cesa, variant_ids = variant_ids, samples = eligible_tumors)
  
  if (running_compound) {
    # sample_calls could include tumors that don't have coverage across both in pair
    tumors_with_v1 = intersect(compound_variants@sample_calls[[v1]], eligible_tumors)
    tumors_with_v2 = intersect(compound_variants@sample_calls[[v2]], eligible_tumors)
    tumors_with_both = intersect(tumors_with_v1, tumors_with_v2)
    tumors_just_v1 = setdiff(tumors_with_v1, tumors_with_both)
    tumors_just_v2 = setdiff(tumors_with_v2, tumors_with_both)
    tumors_with_neither = setdiff(eligible_tumors, c(tumors_with_v1, tumors_with_v2))
    
    v1_rates = all_rates[, ..v1_ids]
    v1_rates = rowSums(v1_rates)
    names(v1_rates) = all_rates$patient_id
    v2_rates = all_rates[, ..v2_ids]
    v2_rates = rowSums(v2_rates)
    names(v2_rates) = all_rates$patient_id
    
    # 2-item lists: first item is baseline rates for v1; second for v2
    with_just_1 = NULL; with_just_2 = NULL; with_both = NULL; with_neither = NULL
    if(length(tumors_just_v1) > 0){
      with_just_1 = list(v1_rates[tumors_just_v1], v2_rates[tumors_just_v1])
    }
    if(length(tumors_just_v2) > 0) {
      with_just_2 = list(v1_rates[tumors_just_v2], v2_rates[tumors_just_v2])
    }
    if (length(tumors_with_both) > 0) {
      with_both = list(v1_rates[tumors_with_both], v2_rates[tumors_with_both])
    }
    if(length(tumors_with_neither) > 0) {
      with_neither = list(v1_rates[tumors_with_neither], v2_rates[tumors_with_neither])
    }
  } else {
    setcolorder(all_rates, c("patient_id", v1, v2))
    setkey(all_rates, "patient_id")
    covered_maf = cesa@maf[eligible_tumors, on = "patient_id", nomatch = NULL]
    tumors_with_v1 = intersect(samples_with(cesa, v1), eligible_tumors)
    tumors_with_v2 = intersect(samples_with(cesa, v2), eligible_tumors)
    tumors_with_both = intersect(tumors_with_v1, tumors_with_v2)
    tumors_just_v1 = setdiff(tumors_with_v1, tumors_with_both)
    tumors_just_v2 = setdiff(tumors_with_v2, tumors_with_both)
    tumors_with_neither = setdiff(eligible_tumors, c(tumors_with_v1, tumors_with_v2))
    
    # 2-item lists: first item is baseline rates for v1; second for v2
    with_just_1 = lapply(as.list(all_rates[tumors_just_v1])[2:3], setNames, tumors_just_v1)
    with_just_2 = lapply(as.list(all_rates[tumors_just_v2])[2:3], setNames, tumors_just_v2)
    with_both = lapply(as.list(all_rates[tumors_with_both])[2:3], setNames, tumors_with_both)
    with_neither = lapply(as.list(all_rates[tumors_with_neither])[2:3], setNames, tumors_with_neither)
  }
  
  validate_optimizer_args(optimizer_args)
  running_default = FALSE
  if(is(model, "character")) {
    if(length(model) != 1 || ! model %in% c("default")) {
      stop("model should specify a built-in selection model (i.e., \"default\") or a custom function factory.")
    } else {
      if (model == 'default') {
        if(! is.null(pval_calc_fn)) {
          warning('pval_calc_fn should usually be NULL when using the default epistasis model.')
          if(! is(pval_calc_fn, 'function') || ! identical(names(formals(pval_calc_fn)), 'fit')) {
            stop('Invalid pval_calc_fn: Must be a function with a single named argument (fit).')
          }
        } else {
          pval_calc_fn = default_epistasis_pvalue_calc
        }
        running_default = TRUE
        lik_factory = pairwise_epistasis_lik
      } else {
        stop("Unrecognized model")
      }
    }
  } else if (! is(model, "function")) {
    stop("model should specify a built-in selection model (\"default\") or a custom function factory.")
  } else {
    lik_factory = model
    if(! is.null(pval_calc_fn) && 
       (! is(pval_calc_fn, 'function') || ! identical(names(formals(pval_calc_fn)), 'fit'))) {
      stop('pval_calc_fn should be NULL or a function that takes a single named argument (fit, the epistatic model fit object).')
    }
  }
  
  if(! is(lik_args, "list")) {
    stop("lik_args should be named list") 
  }
  
  if(length(lik_args) > 0 && is.character(model)){
    if(model %in% c('default')) {
      stop("lik_args aren't used in the chosen model.")
    }
  }
  if(length(lik_args) != uniqueN(names(lik_args))) {
    stop('lik_args should be a named list without repeated names.')
  }
  
  lik_args = c(list(with_just_1 = with_just_1, with_just_2 = with_just_2, with_both = with_both, with_neither = with_neither), lik_args)
  
  # call factory function to get variant-specific likelihood function
  epistasis_lik_fn = do.call(lik_factory, lik_args)
  par_init = formals(epistasis_lik_fn)[[1]]
  names(par_init) = bbmle::parnames(epistasis_lik_fn)
  
  
  # No point of testing epistasis if either variant doesn't appear
  n_total = length(tumors_just_v1) + length(tumors_just_v2) + length(tumors_with_neither) + length(tumors_with_both)
  if (length(tumors_with_v1) == 0 || length(tumors_with_v2) == 0) {
    early_output = list(variant_A = v1, variant_B = v2,
                        nA0 = length(tumors_just_v1),
                        nB0 = length(tumors_just_v2),
                        nAB = length(tumors_with_both),
                        n00 = length(tumors_with_neither),
                        n_total = n_total)
    return(list(summary = early_output, fit = NULL)) # fit list will have a NULL entry
  }
  
  # find optimized selection intensities
  # the selection intensity when some group has 0 variants will be on the lower boundary; will muffle the associated warning
  final_optimizer_args = c(list(minuslogl = epistasis_lik_fn, start = par_init, vecpar = T), optimizer_args)
  if(is.null(final_optimizer_args$control$ndeps)) {
    final_optimizer_args$control$ndeps = rep(.1, length(par_init))
  }
  if(is.null(final_optimizer_args$method)) {
    final_optimizer_args$method = 'L-BFGS-B'
  }
  if(is.null(final_optimizer_args$lower)) {
    final_optimizer_args$lower = 1e-3
  }
  if(is.null(final_optimizer_args$upper)) {
    final_optimizer_args$upper = 1e9
  }
  
  withCallingHandlers(
    {
      fit = do.call(bbmle::mle2, final_optimizer_args)
    },
    warning = function(w) {
      if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
        invokeRestart("muffleWarning")
      }
    }
  )

  # Calculate p-values and other descriptives
  p_values_and_more = list()
  if(! is.null(pval_calc_fn)) {
    # Convenient for fit to know the names of its variants.
    attr(fit, 'v1') = v1
    attr(fit, 'v2') = v2
    p_values_and_more = pval_calc_fn(fit)
  }
  
  params = bbmle::coef(fit)
  if(running_default) {
    si_output = list(ces_A0 = params[1], ces_B0 = params[2],
                     ces_A_on_B = params[3], ces_B_on_A = params[4])
  } else {
    si_output = as.list(bbmle::coef(fit))
  }
  
  variant_ep_results = c(list(variant_A = v1, variant_B = v2),
                         si_output,
                         p_values_and_more,
                         list(nA0 = length(tumors_just_v1), nB0 = length(tumors_just_v2), nAB = length(tumors_with_both),
                              n00 = length(tumors_with_neither), n_total = n_total))
  if(! is.null(conf)) {
    if(running_default) {
      bbmle::parnames(epistasis_lik_fn) = c('ces_A0', 'ces_B0', 'ces_A_on_B', 'ces_B_on_A')
      ci = univariate_si_conf_ints(fit, epistasis_lik_fn, final_optimizer_args$lower, final_optimizer_args$upper, conf)
    } else {
      withCallingHandlers(
        {
          ci = bbmle::confint(fit, level = conf, quietly = TRUE)
        },
        warning = function(w) {
          if (conditionMessage(w) %like% 'reverting from spline to linear approximation') {
            invokeRestart("muffleWarning")
          }
        })
      
      ci_dt = as.data.table(ci, keep.rownames = 'param')
      ci_dt[, low_name := paste('ci_low', conf * 100, param, sep = '_')]
      ci_dt[, high_name := paste('ci_high', conf * 100, param, sep = '_')]
      ci = unlist(S4Vectors::zipup(ci_dt[[2]], ci_dt[[3]]))
      names(ci) = unlist(S4Vectors::zipup(ci_dt$low_name, ci_dt$high_name))
    }
    variant_ep_results = c(variant_ep_results, ci)
  }
  
  # To reduce size, remove copy of CESAnalysis from fit environment.
  rm('cesa', envir = environment(fit))
  
  return(list(summary = variant_ep_results, fit = fit))
}