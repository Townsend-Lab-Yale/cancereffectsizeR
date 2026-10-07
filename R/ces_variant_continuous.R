#' Estimate variant effects with a linear continuous-covariate model
#'
#' Estimates cancer effects for individual or compound variants when the
#' selection intensity is modeled as a linear function of a continuous
#' sample-level covariate:
#'
#' \deqn{\gamma(x) = \beta_0 + \beta_1 x.}
#'
#' For each variant, baseline mutation rates are obtained from the
#' \code{CESAnalysis} object and the parameters of the supplied likelihood
#' model are estimated by maximum likelihood. The fitted selection intensity
#' is constrained to remain positive over the observed covariate range.
#'
#' The likelihood function factory supplied through \code{model} is expected
#' to return a likelihood function whose first two parameters are
#' \code{beta0} and \code{beta1}, in that order. The sample-level continuous
#' covariate is supplied through \code{lik_args$covariate_data}; its
#' \code{covariate_value} column must be coercible to numeric.
#'
#' @inheritParams ces_variant
#'
#' @param model A likelihood-function factory defining the selection model.
#'   The package provides likelihood functions for continuous-covariate
#'   selection that can be used directly; for the linear model, use
#'   \code{sswm_age_lik()}. A custom likelihood-function factory may also be
#'   supplied. See \code{ces_variant()} for the general custom-model interface.
#'
#' @param lik_args Named list of additional arguments passed to \code{model}.
#'   For continuous-covariate inference, this must include
#'   \code{covariate_data}, a table associating
#'   \code{Unique_Patient_Identifier} values with a numeric continuous
#'   covariate in \code{covariate_value}.
#'
#' @param optimizer Character scalar specifying the optimization algorithm.
#'   Currently supported values are \code{"COBYLA"} and \code{"ISRES"},
#'   corresponding to the NLOPT algorithms \code{NLOPT_LN_COBYLA} and
#'   \code{NLOPT_GN_ISRES}, respectively.
#'
#' @param constraint Positive numeric value giving the minimum permitted
#'   selection intensity over the observed covariate range. The linear model
#'   is constrained such that
#'   \code{beta0 + beta1 * x >= constraint} at both the minimum and maximum
#'   observed covariate values.
#'
#' @param conf Nominal confidence level. Set to \code{NULL} to skip confidence
#'   interval calculation. For inference on continuous-model parameters,
#'   confidence intervals can also be calculated separately using the
#'   continuous-model confidence-interval routines.
#'
#' @return A \code{CESAnalysis} object with a new entry appended to
#'   \code{selection_results}. For each analyzed variant, the result contains
#'   the fitted linear-model parameters and maximized log-likelihood, together
#'   with variant annotations and applicable sample-count information.
#'
#' @seealso \code{\link{ces_variant}}
#'
#' @family continuous selection models
#'
#' @export
ces_variant_linear <- function(cesa = NULL,
                               variants = select_variants(cesa, min_freq = 2),
                               samples = character(),
                               model = "default",
                               run_name = "auto",
                               lik_args = list(),
                               optimizer_args = if(identical(model, 'default')) list(method = 'L-BFGS-B', lower = 1e-3, upper = 1e9) else list(),
                               return_fit = FALSE,
                               hold_out_same_gene_samples = "auto",
                               cores = 1,
                               conf = .95,
                               optimizer = c("COBYLA", "ISRES"),
                               constraint = 1e-3)
{
  if(! is.numeric(cores) || length(cores) != 1 || cores - as.integer(cores) != 0 || cores < 1) {
    stop('cores should be 1-length positive integer')
  }

  if(! rlang::is_bool(return_fit)) {
    stop('return_fit should be TRUE/FALSE.')
  }

  validate_optimizer_args(optimizer_args)

  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (length(hold_out_same_gene_samples) == 1) {
    if (! is.logical(hold_out_same_gene_samples)) {
      if (identical(hold_out_same_gene_samples, "auto")) {
        hold_out_same_gene_samples = is(variants, "data.table")
      } else {
        stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
      }
    }
  } else {
    stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
  }

  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if(run_name %in% names(cesa@selection_results)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if (run_name == "auto") {
    # sequentially name results, handling nefarious run naming
    run_number = length(cesa@selection_results) + 1
    run_name = paste0('selection.', run_number)
    while(run_name %in% names(cesa@selection_results)) {
      run_number = run_number + 1
      run_name = paste0('variant_effects_', run_number)
    }
  }

  if(is(model, "character")) {
    # old names were basic, sswm-sequential (no one could remember if hyphen or underscore)
    model = tolower(model)
    model[model %in% c('sswm', 'default')] = 'basic'
    model[model %like% 'sswm[-_]sequential'] = 'sequential'
    if(length(model) != 1 || ! model %in% c("basic", "sequential")) {
      stop("model should specify a built-in selection model (i.e., \"default\") or a custom function factory.")
    } else {
      if (model == 'basic') {
        lik_factory = sswm_lik
      } else if(model == 'sequential') {
        lik_factory = sswm_sequential_lik
      } else {
        stop("Unrecognized model")
      }
    }
  } else if (! is(model, "function")) {
    stop("model should specify a built-in selection model (\"default\") or a custom function factory.")
  } else {
    lik_factory = model
  }

  if(! is(lik_args, "list")) {
    stop("lik args should be named list")
  }

  if(length(lik_args) != uniqueN(names(lik_args))) {
    stop('lik_args should be a named list without repeated names.')
  }

  optimizer <- match.arg(optimizer)

  samples = select_samples(cesa, samples)
  if(samples[, .N] < cesa@samples[, .N]) {
    num_excluded = cesa@samples[, .N] - samples[, .N]
    pretty_message(paste0("Note that ", format(num_excluded, big.mark = ','), " samples are being excluded from selection inference."))
  }

  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())

  # Set keys in case they've been lost
  mutations = cesa@mutations
  if(! is.null(conf)) {
    if(is(model, 'function')) {
      if(! rlang::is_scalar_double(conf) || conf != .95) {
        warning('conf is ignored when running a custom model.')
      }
      conf = NULL
    } else {
      if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
        stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
      }
    }
  }

  running_compound = FALSE

  # If an input variant table came directly from select_variants() and the variants are non-overlapping,
  # just accept the table. Otherwise, re-select the variants with the variant_id field.
  if (is(variants, "data.table")) {
    if(! "variant_id" %in% names(variants)) {
      stop("variants table is missing a variant_id column. Typically, variants is generated using select_variants().")
    }
    nonoverlapping = attr(variants, "nonoverlapping")
    if(is.null(nonoverlapping)) {
      if ('variant_id' %in% names(variants)) {
        pretty_message('Taking variants from variant_id column of input table....')
      }
    } else if(! identical(nonoverlapping, TRUE)) {
      stop("Input variants table may contain overlapping variants; re-run select_variants() to get a non-overlapping table.")
    }

    # re-select variants for maximum safety
    variants = select_variants(cesa, variant_ids = variants[, variant_id])

  } else if (is(variants, "CompoundVariantSet")) {
    running_compound = TRUE
    if (cesa@advanced$uid != variants@cesa_uid) {
      stop("Input CompoundVariantSet does not appear to derive from the input CESAnalysis.")
    }
    if (cesa@samples[, .N] != variants@cesa_num_samples) {
      stop("The number of samples in the CESAnalysis has changed since the CompoundVariantSet was created. ",
           "Please re-generate it.")
    }
    compound_variants = variants
    variants = select_variants(cesa, variant_ids = compound_variants@snvs$snv_id, include_subvariants = TRUE)

    # copy in compound variant names and overwrite covered_in with value of shared_cov
    variants = variants[compound_variants@snvs, compound_name := compound_name, on = c(variant_id = "snv_id")]
    if(variants[, .N] != compound_variants@snvs[, .N]) {
      stop("Internal error: select_variants() didn't return variant info 1-to-1 with compound input.")
    }
    variants[compound_variants@compounds, covered_in := shared_cov, on = "compound_name"]
  } else {
    stop("variants expected to be a variant table (from select_variants(), usually) or a CompoundVariantSet")
  }
  if (variants[, .N] == 0) {
    stop("There are no variants in the input!")
  }
  aac_ids = variants[variant_type == "aac", variant_id]

  # By noncoding, we just mean that SIs are calculated at the SNV site rather than at the AAC level,
  # regardless of whether there's a CDS annotation.
  noncoding_snv_ids = variants[variant_type == "snv", variant_id]

  if(length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", call. = F)
  }

  # identify mutations by nearest gene(s)
  maf = cesa@maf[samples$Unique_Patient_Identifier, on = "Unique_Patient_Identifier", nomatch = NULL]
  tmp = unique(maf[, .(gene = unlist(genes)), by = "Unique_Patient_Identifier"])[, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  tumors_with_variants_by_gene = list2env(tumors_with_variants_by_gene)

  # identify mutations by sample
  snv_aac_of_interest = cesa@mutations$aac_snv_key[aac_ids, on = 'aac_id']
  tmp = maf[snv_aac_of_interest, .(variant_id), on = c(variant_id = 'snv_id'),
            by = "Unique_Patient_Identifier", nomatch = NULL]
  tmp[snv_aac_of_interest, aac_id := aac_id, on = c(variant_id = 'snv_id')]
  tmp = tmp[, .(samples = list(unique(Unique_Patient_Identifier))), by = "aac_id"]
  samples_by_aac = setNames(tmp$samples, tmp$aac_id)

  setkey(maf, "variant_id")
  # need nomatch because some noncoding SNVs may not be present in the samples
  tmp = maf[noncoding_snv_ids, variant_id, by = "Unique_Patient_Identifier", nomatch = NULL][, .(samples = list(Unique_Patient_Identifier)), by = "variant_id"]
  samples_by_snv = tmp$samples
  names(samples_by_snv) = tmp$variant_id
  samples_by_variant = list2env(c(samples_by_aac, samples_by_snv))


  setkey(samples, "covered_regions")
  # These are WGS samples with purportedly whole-genome coverage.
  # That is, for better or worse, assuming that any variant can be found in these samples.
  # (Trimmed-interval WGS samples will have coverage = "genome" and covered_regions != "genome.")
  genome_wide_cov_samples = samples["genome", Unique_Patient_Identifier, nomatch = NULL]


  # Will process variants by coverage group (i.e., groups of variants that have the same tumors covering them)
  selection_results = NULL

  all_coverage = rbind(cesa@mutations$snv[, .(variant_id = snv_id, covered_in)],
                       cesa@mutations$amino_acid_change[, .(variant_id = aac_id, covered_in)])
  variants[all_coverage, covered_in := covered_in, on = 'variant_id']
  coverage_groups = unique(variants$covered_in)
  num_coverage_groups = length(coverage_groups)

  setkey(maf, "variant_id")
  setkey(variants, "variant_id")

  selection_results = lapply(1:length(coverage_groups), function(i) {
    coverage_group = coverage_groups[[i]]
    # one coverage group may be NA, for variants that are not covered by any specific covered_regions
    if(length(coverage_group) == 1 && is.na(coverage_group)) {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, NA_character_)))]
    } else {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, coverage_group)))]
    }
    message(sprintf("Preparing to calculate cancer effects (batch %i of %i)...", i, num_coverage_groups))
    if(is.null(coverage_group)) {
      covered_samples = genome_wide_cov_samples
    } else {
      covered_samples = c(samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL], genome_wide_cov_samples)
    }
    variants[curr_variants$variant_id, num_covered_and_in_samples := length(covered_samples), on = 'variant_id']

    # When not all samples are used, it's possible that no sampels will have coverage at the input variants.
    if(length(covered_samples) == 0) {
      # message("Skipped batch ", i, " because no samples had coverage at the variant sites in the batch.")
      return(list(data.table(), NULL))
    }
    # rough size of baseline rates data.table in bytes, if all included in one table
    work_size = length(covered_samples) * curr_variants[,.N] * 8

    # we divide into subgroups to cap baseline rates table at around 1 GB
    num_proc_groups = ceiling(work_size / 1e9)
    curr_variants[, subgroup := ceiling(num_proc_groups * 1:.N / .N)]

    if (running_compound) {
      # can't have subvariants of the same compound variant ending up in different subgroups
      curr_variants[, subgroup := rep.int(subgroup[1], .N), by = "compound_name"]
      num_proc_groups = max(curr_variants$subgroup) # rarely, last subgroup dropped by above
    }


    curr_results = lapply(1:num_proc_groups, function(j) {
      if (num_proc_groups > 1) {
        message(sprintf("Working on sub-batch %i of %i...", j, num_proc_groups))
      }
      curr_subgroup = curr_variants[subgroup == j]
      aac_ids = curr_subgroup[variant_type == "aac", variant_id]
      snv_ids = curr_subgroup[variant_type == "snv", variant_id]

      baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, snv_ids = snv_ids, samples = covered_samples)
      # put gene(s) by variant into env for quick access
      gene_lookup = curr_subgroup[, all_genes]
      names(gene_lookup) = curr_subgroup[, variant_id]
      gene_lookup = list2env(gene_lookup)

      # function to run MLE on given variant_id (vector of IDs for compound variants)
      process_variant = function(variant_id) {
        if(running_compound) {
          compound_id = variant_id
          # Important: sample_calls includes samples that are not in shared coverage
          tumors_with_variant = intersect(compound_variants@sample_calls[[compound_id]], covered_samples)
          current_snvs = compound_variants@snvs[compound_name == compound_id]
          all_genes = current_snvs[, unique(unlist(genes))]
          variant_id = current_snvs$snv_id
          rates = baseline_rates[, ..variant_id]
          # Sum Poisson rates across variants
          rates = rowSums(rates)
        } else {
          tumors_with_variant = samples_by_variant[[variant_id]]
          all_genes = gene_lookup[[variant_id]]
          rates = baseline_rates[, ..variant_id][[1]]
        }
        names(rates) = baseline_rates[, Unique_Patient_Identifier]

        # usually but not always just 1 gene when not compound (when compound, anything possible)
        if (hold_out_same_gene_samples) {
          if (length(all_genes) == 1) {
            if(is.na(all_genes)) {
              tumors_with_gene_mutated = tumors_with_variant
            } else {
              tumors_with_gene_mutated = tumors_with_variants_by_gene[[all_genes]]
            }
          } else {
            tumors_with_gene_mutated = unique(unlist(sapply(all_genes, function(x) tumors_with_variants_by_gene[[x]])))
          }
          tumors_without = setdiff(covered_samples, tumors_with_gene_mutated)
        } else {
          tumors_without = setdiff(covered_samples, tumors_with_variant)
        }
        rates_tumors_with = rates[tumors_with_variant]
        rates_tumors_without = rates[tumors_without]


        lik_args = c(list(rates_tumors_with = rates_tumors_with, rates_tumors_without = rates_tumors_without),
                     lik_args)
        fn = do.call(lik_factory, lik_args)

        # Linear model, selection intensity > 0 over the observed covariate range
        min_covariate <- min(as.numeric(lik_args[["covariate_data"]]$covariate_value))
        max_covariate <- max(as.numeric(lik_args[["covariate_data"]]$covariate_value))

        eval_g0 <- function(params){
          beta0 <- params[1]
          beta1 <- params[2]
          return(rbind(-beta0 - min_covariate * beta1 + constraint, -beta0 - max_covariate * beta1 + constraint))
        }

        # Initial parameters
        par_init = formals(fn)[[1]]
        names(par_init) = bbmle::parnames(fn)

        if (optimizer == "COBYLA") {
          print("Method used: nloptr COBYLA")
          fit <- nloptr(x0=par_init,
                        eval_f=fn,
                        lb = c(-1e6,-1e7),
                        ub = c(1e9,1e7),
                        eval_g_ineq = eval_g0,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1.0e-8,
                                    "maxeval"=10000))

          selection_intensity <- fit$solution
          loglikelihood <- -fit$objective
          print(paste0("loglik: ", loglikelihood))

        } else if (optimizer == "ISRES"){
          print("Method used: nloptr ISRES")
          fit <- nloptr(x0=par_init,
                        eval_f=fn,
                        lb = c(-1e6,-1e7),
                        ub = c(1e9,1e7),
                        eval_g_ineq = eval_g0,
                        opts = list(algorithm    = "NLOPT_GN_ISRES",
                                    "xtol_rel"=1.0e-8,
                                    "maxeval"=10000))

          selection_intensity <- fit$solution
          loglikelihood <- -fit$objective
          print(paste0("loglik: ", loglikelihood))
        }

        names(selection_intensity) <- names(par_init)

        if (running_compound) {
          variant_id = compound_id
        }
        variant_output = c(list(variant_id = variant_id),
                           as.list(selection_intensity),
                           list(loglikelihood = loglikelihood))

        if(is.character(model) || is.null(lik_args$covariate_data)){

          if(is.character(model)){
            if (model == 'basic') {
              # Record counts of total samples included in inference and included samples with the variant.
              # This may vary from the naive output of variant_counts() due to issues of sample coverage and
              # (by default) the use of hold_out_same_gene_samples = TRUE.

              num_samples_with = length(tumors_with_variant)
              num_samples_total = num_samples_with + length(tumors_without)
              variant_output = c(variant_output, list(included_with_variant = num_samples_with,
                                                      included_total = num_samples_total))
            }
          }

          if(is(model, "function") && is.null(lik_args$covariate_data)){
            num_samples_with = length(tumors_with_variant)
            num_samples_total = num_samples_with + length(tumors_without)
            variant_output = c(variant_output, list(included_with_variant = num_samples_with,
                                                    included_total = num_samples_total))
          }
        }
        if(! is.null(conf)) {
          min_value = ifelse(is.null(final_optimizer_args$lower), -Inf, final_optimizer_args$lower)
          max_value = ifelse(is.null(final_optimizer_args$upper), Inf, final_optimizer_args$upper)
          variant_output = c(variant_output,
                             univariate_si_conf_ints(fit, fn, min_value, max_value, conf))
        }
        return(list(summary = variant_output, fit = if (return_fit) fit else NULL))
      }

      if (running_compound) {
        variants_to_run = curr_subgroup[, unique(compound_name)]
      } else {
        variants_to_run = curr_subgroup$variant_id
      }
      message("Calculating cancer effects...")
      subgroup_results = pbapply::pblapply(variants_to_run, process_variant, cl = cores)
      subgroup_selection = rbindlist(lapply(subgroup_results, '[[', 1))
      subgroup_fit = lapply(subgroup_results, '[[', 2)
      return(list(subgroup_selection, subgroup_fit))
    })
    group_selection = rbindlist(lapply(curr_results, '[[', 1))
    group_fit = lapply(curr_results, '[[', 2)
    return(list(group_selection, group_fit))
  })

  fits = unlist(lapply(selection_results, '[[', 2))
  selection_results = rbindlist(lapply(selection_results, '[[', 1))

  if(selection_results[, .N] == 0) {
    msg = paste0("No selection inference was performed, so returning CESAnalysis unaltered without selection output. ",
                 "Perhaps none of the variants had coverage in the specified samples?")
    pretty_message(msg)
    return(cesa)
  }

  if (running_compound) {
    selection_results[, variant_type := "compound"]
    setnames(selection_results, "variant_id", "variant_name")
    selection_results = selection_results[compound_variants@compounds, on = c(variant_name = "compound_name")]
    setattr(selection_results, "is_compound", TRUE)
    setcolorder(selection_results, c("variant_name", "variant_type"))
  } else {
    # Fill in top-priority gene for variants that are in coding regions or essential splice (or within 1 bp)
    # Other variants will get NA gene.
    selection_results[variants, c("variant_type", "variant_name", "gene", "intergenic") :=
                        list(variant_type, variant_name, gene, intergenic), on = "variant_id"]
    selection_results[intergenic == T, gene := NA]
    selection_results$intergenic = NULL
    setattr(selection_results, "is_compound", FALSE)
    setcolorder(selection_results, c("variant_name", "variant_type", "gene"))
    setcolorder(selection_results, c(setdiff(names(selection_results), 'variant_id'), 'variant_id'))
  }


  if(hold_out_same_gene_samples == FALSE) {
    selection_results[, held_out := 0]
  } else {
    if('included_total' %in% names(selection_results)) {
      if (running_compound) {
        num_eligible_by_comp = sapply(compound_variants$definitions,
                                      function(x) variants[x, min(num_covered_and_in_samples)], USE.NAMES = TRUE)
        selection_results[, held_out := num_eligible_by_comp[variant_name] - included_total]
      } else {
        selection_results[variants, held_out := num_covered_and_in_samples - included_total, on = 'variant_id']
      }
    }
  }

  if(all(c('held_out', 'included_total') %in% names(selection_results))) {
    selection_results[, uncovered := samples[, .N] - included_total - held_out]
  }

  if(return_fit) {
    fits = lapply(fits, function(x) {
      x@call.orig = call('[not shown]')
      parent.env(environment(x)) = emptyenv() # otherwise, whole CESAnalysis will be embedded via parent
      return(x)
    })
    setattr(selection_results, 'fit', fits)
  }
  curr_results = list(selection_results)
  names(curr_results) = run_name


  cesa@selection_results = c(cesa@selection_results, curr_results)
  return(cesa)
}





#' Estimate variant effects with a logistic continuous-covariate model
#'
#' Estimates cancer effects for individual or compound variants when the
#' selection intensity varies with a continuous sample-level covariate,
#' according to a logistic curve:
#'
#' \deqn{\gamma(x) = \frac{\exp(L)}{1 + \exp[-k(x - m)]}.}
#'
#' Here, \eqn{\exp(L)} determines the upper asymptote, \eqn{k} determines the
#' steepness or growth/decay rate of the transition, and \eqn{m} is the
#' midpoint of the transition. A positive \eqn{k} gives an increasing curve and
#' a negative \eqn{k} gives a decreasing curve.
#'
#' Baseline mutation rates are obtained from the \code{CESAnalysis} object,
#' and model parameters are estimated separately for each variant by maximum
#' likelihood.
#'
#'
#' @inheritParams ces_variant
#'
#' @param model A likelihood-function factory defining the selection model.
#'   The package provides likelihood functions for continuous-covariate
#'   selection that can be used directly; for the logistic model, use
#'   \code{sswm_age_lik_logistic()}. A custom likelihood-function factory
#'   may also be supplied. See \code{ces_variant()} for the general custom-model
#'   interface.
#'
#' @param lik_args Named list of additional arguments passed to \code{model}.
#'   For continuous-covariate inference, this must include
#'   \code{covariate_data}, a table associating
#'   \code{Unique_Patient_Identifier} values with a numeric continuous
#'   covariate in \code{covariate_value}.
#'
#' @param optimizer Character scalar specifying the optimizer. Supported values
#'   are \code{"bbmle"} and \code{"COBYLA"}. \code{"bbmle"} fits the model
#'   using \code{bbmle::mle2()}, whereas \code{"COBYLA"} uses the constrained
#'   NLOPT COBYLA algorithm.
#'
#' @param optimizer_args Named list of additional arguments passed to
#'   \code{bbmle::mle2()} when \code{optimizer = "bbmle"}.
#'
#' @param constraint Positive numeric value giving the minimum permitted
#'   fitted selection intensity over the observed covariate range when using
#'   constrained optimization.
#'
#' @param conf Nominal confidence level. Set to \code{NULL} to skip confidence
#'   interval calculation. For inference on continuous-model parameters,
#'   confidence intervals can also be calculated separately using the
#'   continuous-model confidence-interval routines.
#'
#' @return A \code{CESAnalysis} object with a new entry appended to
#'   \code{selection_results}. For each analyzed variant, the result contains
#'   the fitted logistic-model parameters ordered as \code{L}, \code{k}, and
#'   \code{m} and maximized log-likelihood, together with variant annotations
#'   and applicable sample-count information.
#'
#' @seealso \code{\link{ces_variant}}
#'
#' @family continuous selection models
#'
#' @export
ces_variant_logistic <- function(cesa = NULL,
                                 variants = select_variants(cesa, min_freq = 2),
                                 samples = character(),
                                 model = "default",
                                 run_name = "auto",
                                 lik_args = list(),
                                 optimizer_args = if(identical(model, 'default')) list(method = 'L-BFGS-B', lower = 1e-3, upper = 1e9) else list(),
                                 return_fit = FALSE,
                                 hold_out_same_gene_samples = "auto",
                                 cores = 1,
                                 conf = .95,
                                 optimizer = c("bbmle", "COBYLA"),
                                 constraint = 1e-3)
{
  if(! is.numeric(cores) || length(cores) != 1 || cores - as.integer(cores) != 0 || cores < 1) {
    stop('cores should be 1-length positive integer')
  }

  if(! rlang::is_bool(return_fit)) {
    stop('return_fit should be TRUE/FALSE.')
  }

  validate_optimizer_args(optimizer_args)

  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (length(hold_out_same_gene_samples) == 1) {
    if (! is.logical(hold_out_same_gene_samples)) {
      if (identical(hold_out_same_gene_samples, "auto")) {
        hold_out_same_gene_samples = is(variants, "data.table")
      } else {
        stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
      }
    }
  } else {
    stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
  }

  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if(run_name %in% names(cesa@selection_results)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if (run_name == "auto") {
    # sequentially name results, handling nefarious run naming
    run_number = length(cesa@selection_results) + 1
    run_name = paste0('selection.', run_number)
    while(run_name %in% names(cesa@selection_results)) {
      run_number = run_number + 1
      run_name = paste0('variant_effects_', run_number)
    }
  }

  if(is(model, "character")) {
    # old names were basic, sswm-sequential (no one could remember if hyphen or underscore)
    model = tolower(model)
    model[model %in% c('sswm', 'default')] = 'basic'
    model[model %like% 'sswm[-_]sequential'] = 'sequential'
    if(length(model) != 1 || ! model %in% c("basic", "sequential")) {
      stop("model should specify a built-in selection model (i.e., \"default\") or a custom function factory.")
    } else {
      if (model == 'basic') {
        lik_factory = sswm_lik
      } else if(model == 'sequential') {
        lik_factory = sswm_sequential_lik
      } else {
        stop("Unrecognized model")
      }
    }
  } else if (! is(model, "function")) {
    stop("model should specify a built-in selection model (\"default\") or a custom function factory.")
  } else {
    lik_factory = model
  }

  if(! is(lik_args, "list")) {
    stop("lik args should be named list")
  }

  if(length(lik_args) != uniqueN(names(lik_args))) {
    stop('lik_args should be a named list without repeated names.')
  }

  optimizer <- match.arg(optimizer)

  samples = select_samples(cesa, samples)
  if(samples[, .N] < cesa@samples[, .N]) {
    num_excluded = cesa@samples[, .N] - samples[, .N]
    pretty_message(paste0("Note that ", format(num_excluded, big.mark = ','), " samples are being excluded from selection inference."))
  }

  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())

  # Set keys in case they've been lost
  mutations = cesa@mutations
  if(! is.null(conf)) {
    if(is(model, 'function')) {
      if(! rlang::is_scalar_double(conf) || conf != .95) {
        warning('conf is ignored when running a custom model.')
      }
      conf = NULL
    } else {
      if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
        stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
      }
    }
  }

  running_compound = FALSE

  # If an input variant table came directly from select_variants() and the variants are non-overlapping,
  # just accept the table. Otherwise, re-select the variants with the variant_id field.
  if (is(variants, "data.table")) {
    if(! "variant_id" %in% names(variants)) {
      stop("variants table is missing a variant_id column. Typically, variants is generated using select_variants().")
    }
    nonoverlapping = attr(variants, "nonoverlapping")
    if(is.null(nonoverlapping)) {
      if ('variant_id' %in% names(variants)) {
        pretty_message('Taking variants from variant_id column of input table....')
      }
    } else if(! identical(nonoverlapping, TRUE)) {
      stop("Input variants table may contain overlapping variants; re-run select_variants() to get a non-overlapping table.")
    }

    # re-select variants for maximum safety
    variants = select_variants(cesa, variant_ids = variants[, variant_id])

  } else if (is(variants, "CompoundVariantSet")) {
    running_compound = TRUE
    if (cesa@advanced$uid != variants@cesa_uid) {
      stop("Input CompoundVariantSet does not appear to derive from the input CESAnalysis.")
    }
    if (cesa@samples[, .N] != variants@cesa_num_samples) {
      stop("The number of samples in the CESAnalysis has changed since the CompoundVariantSet was created. ",
           "Please re-generate it.")
    }
    compound_variants = variants
    variants = select_variants(cesa, variant_ids = compound_variants@snvs$snv_id, include_subvariants = TRUE)

    # copy in compound variant names and overwrite covered_in with value of shared_cov
    variants = variants[compound_variants@snvs, compound_name := compound_name, on = c(variant_id = "snv_id")]
    if(variants[, .N] != compound_variants@snvs[, .N]) {
      stop("Internal error: select_variants() didn't return variant info 1-to-1 with compound input.")
    }
    variants[compound_variants@compounds, covered_in := shared_cov, on = "compound_name"]
  } else {
    stop("variants expected to be a variant table (from select_variants(), usually) or a CompoundVariantSet")
  }
  if (variants[, .N] == 0) {
    stop("There are no variants in the input!")
  }
  aac_ids = variants[variant_type == "aac", variant_id]

  # By noncoding, we just mean that SIs are calculated at the SNV site rather than at the AAC level,
  # regardless of whether there's a CDS annotation.
  noncoding_snv_ids = variants[variant_type == "snv", variant_id]

  if(length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", call. = F)
  }

  # identify mutations by nearest gene(s)
  maf = cesa@maf[samples$Unique_Patient_Identifier, on = "Unique_Patient_Identifier", nomatch = NULL]
  tmp = unique(maf[, .(gene = unlist(genes)), by = "Unique_Patient_Identifier"])[, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  tumors_with_variants_by_gene = list2env(tumors_with_variants_by_gene)

  # identify mutations by sample
  snv_aac_of_interest = cesa@mutations$aac_snv_key[aac_ids, on = 'aac_id']
  tmp = maf[snv_aac_of_interest, .(variant_id), on = c(variant_id = 'snv_id'),
            by = "Unique_Patient_Identifier", nomatch = NULL]
  tmp[snv_aac_of_interest, aac_id := aac_id, on = c(variant_id = 'snv_id')]
  tmp = tmp[, .(samples = list(unique(Unique_Patient_Identifier))), by = "aac_id"]
  samples_by_aac = setNames(tmp$samples, tmp$aac_id)

  setkey(maf, "variant_id")
  # need nomatch because some noncoding SNVs may not be present in the samples
  tmp = maf[noncoding_snv_ids, variant_id, by = "Unique_Patient_Identifier", nomatch = NULL][, .(samples = list(Unique_Patient_Identifier)), by = "variant_id"]
  samples_by_snv = tmp$samples
  names(samples_by_snv) = tmp$variant_id
  samples_by_variant = list2env(c(samples_by_aac, samples_by_snv))


  setkey(samples, "covered_regions")
  # These are WGS samples with purportedly whole-genome coverage.
  # That is, for better or worse, assuming that any variant can be found in these samples.
  # (Trimmed-interval WGS samples will have coverage = "genome" and covered_regions != "genome.")
  genome_wide_cov_samples = samples["genome", Unique_Patient_Identifier, nomatch = NULL]


  # Will process variants by coverage group (i.e., groups of variants that have the same tumors covering them)
  selection_results = NULL

  all_coverage = rbind(cesa@mutations$snv[, .(variant_id = snv_id, covered_in)],
                       cesa@mutations$amino_acid_change[, .(variant_id = aac_id, covered_in)])
  variants[all_coverage, covered_in := covered_in, on = 'variant_id']
  coverage_groups = unique(variants$covered_in)
  num_coverage_groups = length(coverage_groups)

  setkey(maf, "variant_id")
  setkey(variants, "variant_id")

  selection_results = lapply(1:length(coverage_groups), function(i) {
    coverage_group = coverage_groups[[i]]
    # one coverage group may be NA, for variants that are not covered by any specific covered_regions
    if(length(coverage_group) == 1 && is.na(coverage_group)) {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, NA_character_)))]
    } else {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, coverage_group)))]
    }
    message(sprintf("Preparing to calculate cancer effects (batch %i of %i)...", i, num_coverage_groups))
    if(is.null(coverage_group)) {
      covered_samples = genome_wide_cov_samples
    } else {
      covered_samples = c(samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL], genome_wide_cov_samples)
    }
    variants[curr_variants$variant_id, num_covered_and_in_samples := length(covered_samples), on = 'variant_id']

    # When not all samples are used, it's possible that no sampels will have coverage at the input variants.
    if(length(covered_samples) == 0) {
      # message("Skipped batch ", i, " because no samples had coverage at the variant sites in the batch.")
      return(list(data.table(), NULL))
    }
    # rough size of baseline rates data.table in bytes, if all included in one table
    work_size = length(covered_samples) * curr_variants[,.N] * 8

    # we divide into subgroups to cap baseline rates table at around 1 GB
    num_proc_groups = ceiling(work_size / 1e9)
    curr_variants[, subgroup := ceiling(num_proc_groups * 1:.N / .N)]

    if (running_compound) {
      # can't have subvariants of the same compound variant ending up in different subgroups
      curr_variants[, subgroup := rep.int(subgroup[1], .N), by = "compound_name"]
      num_proc_groups = max(curr_variants$subgroup) # rarely, last subgroup dropped by above
    }


    curr_results = lapply(1:num_proc_groups, function(j) {
      if (num_proc_groups > 1) {
        message(sprintf("Working on sub-batch %i of %i...", j, num_proc_groups))
      }
      curr_subgroup = curr_variants[subgroup == j]
      aac_ids = curr_subgroup[variant_type == "aac", variant_id]
      snv_ids = curr_subgroup[variant_type == "snv", variant_id]

      baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, snv_ids = snv_ids, samples = covered_samples)
      # put gene(s) by variant into env for quick access
      gene_lookup = curr_subgroup[, all_genes]
      names(gene_lookup) = curr_subgroup[, variant_id]
      gene_lookup = list2env(gene_lookup)

      # function to run MLE on given variant_id (vector of IDs for compound variants)
      process_variant = function(variant_id) {
        if(running_compound) {
          compound_id = variant_id
          # Important: sample_calls includes samples that are not in shared coverage
          tumors_with_variant = intersect(compound_variants@sample_calls[[compound_id]], covered_samples)
          current_snvs = compound_variants@snvs[compound_name == compound_id]
          all_genes = current_snvs[, unique(unlist(genes))]
          variant_id = current_snvs$snv_id
          rates = baseline_rates[, ..variant_id]
          # Sum Poisson rates across variants
          rates = rowSums(rates)
        } else {
          tumors_with_variant = samples_by_variant[[variant_id]]
          all_genes = gene_lookup[[variant_id]]
          rates = baseline_rates[, ..variant_id][[1]]
        }
        names(rates) = baseline_rates[, Unique_Patient_Identifier]

        # usually but not always just 1 gene when not compound (when compound, anything possible)
        if (hold_out_same_gene_samples) {
          if (length(all_genes) == 1) {
            if(is.na(all_genes)) {
              tumors_with_gene_mutated = tumors_with_variant
            } else {
              tumors_with_gene_mutated = tumors_with_variants_by_gene[[all_genes]]
            }
          } else {
            tumors_with_gene_mutated = unique(unlist(sapply(all_genes, function(x) tumors_with_variants_by_gene[[x]])))
          }
          tumors_without = setdiff(covered_samples, tumors_with_gene_mutated)
        } else {
          tumors_without = setdiff(covered_samples, tumors_with_variant)
        }
        rates_tumors_with = rates[tumors_with_variant]
        rates_tumors_without = rates[tumors_without]


        lik_args = c(list(rates_tumors_with = rates_tumors_with, rates_tumors_without = rates_tumors_without),
                     lik_args)
        fn = do.call(lik_factory, lik_args)

        # Logistic model can use bbmle:mle2(method = "Nelder-Mead"), since we don't
        # have any constraint when optimizing. So need to change ces_variant()

        if (optimizer == "bbmle"){
          par_init = formals(fn)[[1]]
          names(par_init) = bbmle::parnames(fn)
          # find optimized selection intensities
          # the selection intensity for any stage that has 0 variants will be on the lower boundary; will muffle the associated warning
          final_optimizer_args = c(list(minuslogl = fn, start = par_init, vecpar = T), optimizer_args)

          withCallingHandlers(
            {
              fit = do.call(bbmle::mle2, final_optimizer_args)
            },
            warning = function(w) {
              if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
                invokeRestart("muffleWarning")
              }
              if (grepl(x = conditionMessage(w), pattern = "convergence failure")) {
                # a little dangerous to muffle, but so far these warnings are
                # quite rare and have been harmless
                invokeRestart("muffleWarning")
              }
            }
          )
          selection_intensity = bbmle::coef(fit)
          loglikelihood = as.numeric(bbmle::logLik(fit))

        } else if (optimizer == "COBYLA"){
          min_covariate <- min(as.numeric(lik_args[["covariate_data"]]$covariate_value))
          max_covariate <- max(as.numeric(lik_args[["covariate_data"]]$covariate_value))

          eval_g0 <- function(params){
            L <- params[1]
            k <- params[2]
            m <- params[3]
            return(rbind(-(exp(L) / (1 + exp(-k * (min_covariate - m)))) + constraint,
                         -(exp(L) / (1 + exp(-k * (max_covariate - m)))) + constraint))
          }

          # Initial parameters
          par_init = formals(fn)[[1]]
          names(par_init) = bbmle::parnames(fn)

          print("Method used: nloptr COBYLA")
          fit <- nloptr(x0=par_init,
                        eval_f=fn,
                        lb = c(-100,-100,-1e3),
                        ub = c(100,100,1e3),
                        eval_g_ineq = eval_g0,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1.0e-8,
                                    "maxeval"=10000))

          selection_intensity <- fit$solution
          loglikelihood <- -fit$objective
        }

        names(selection_intensity) <- names(par_init)

        if (running_compound) {
          variant_id = compound_id
        }
        variant_output = c(list(variant_id = variant_id),
                           as.list(selection_intensity),
                           list(loglikelihood = loglikelihood))

        if(is.character(model) || is.null(lik_args$covariate_data)){

          if(is.character(model)){
            if (model == 'basic') {
              # Record counts of total samples included in inference and included samples with the variant.
              # This may vary from the naive output of variant_counts() due to issues of sample coverage and
              # (by default) the use of hold_out_same_gene_samples = TRUE.

              num_samples_with = length(tumors_with_variant)
              num_samples_total = num_samples_with + length(tumors_without)
              variant_output = c(variant_output, list(included_with_variant = num_samples_with,
                                                      included_total = num_samples_total))
            }
          }

          if(is(model, "function") && is.null(lik_args$covariate_data)){
            num_samples_with = length(tumors_with_variant)
            num_samples_total = num_samples_with + length(tumors_without)
            variant_output = c(variant_output, list(included_with_variant = num_samples_with,
                                                    included_total = num_samples_total))
          }
        }
        if(! is.null(conf)) {
          min_value = ifelse(is.null(final_optimizer_args$lower), -Inf, final_optimizer_args$lower)
          max_value = ifelse(is.null(final_optimizer_args$upper), Inf, final_optimizer_args$upper)
          variant_output = c(variant_output,
                             univariate_si_conf_ints(fit, fn, min_value, max_value, conf))
        }
        return(list(summary = variant_output, fit = if (return_fit) fit else NULL))
      }

      if (running_compound) {
        variants_to_run = curr_subgroup[, unique(compound_name)]
      } else {
        variants_to_run = curr_subgroup$variant_id
      }
      message("Calculating cancer effects...")
      subgroup_results = pbapply::pblapply(variants_to_run, process_variant, cl = cores)
      subgroup_selection = rbindlist(lapply(subgroup_results, '[[', 1))
      subgroup_fit = lapply(subgroup_results, '[[', 2)
      return(list(subgroup_selection, subgroup_fit))
    })
    group_selection = rbindlist(lapply(curr_results, '[[', 1))
    group_fit = lapply(curr_results, '[[', 2)
    return(list(group_selection, group_fit))
  })

  fits = unlist(lapply(selection_results, '[[', 2))
  selection_results = rbindlist(lapply(selection_results, '[[', 1))

  if(selection_results[, .N] == 0) {
    msg = paste0("No selection inference was performed, so returning CESAnalysis unaltered without selection output. ",
                 "Perhaps none of the variants had coverage in the specified samples?")
    pretty_message(msg)
    return(cesa)
  }

  if (running_compound) {
    selection_results[, variant_type := "compound"]
    setnames(selection_results, "variant_id", "variant_name")
    selection_results = selection_results[compound_variants@compounds, on = c(variant_name = "compound_name")]
    setattr(selection_results, "is_compound", TRUE)
    setcolorder(selection_results, c("variant_name", "variant_type"))
  } else {
    # Fill in top-priority gene for variants that are in coding regions or essential splice (or within 1 bp)
    # Other variants will get NA gene.
    selection_results[variants, c("variant_type", "variant_name", "gene", "intergenic") :=
                        list(variant_type, variant_name, gene, intergenic), on = "variant_id"]
    selection_results[intergenic == T, gene := NA]
    selection_results$intergenic = NULL
    setattr(selection_results, "is_compound", FALSE)
    setcolorder(selection_results, c("variant_name", "variant_type", "gene"))
    setcolorder(selection_results, c(setdiff(names(selection_results), 'variant_id'), 'variant_id'))
  }


  if(hold_out_same_gene_samples == FALSE) {
    selection_results[, held_out := 0]
  } else {
    if('included_total' %in% names(selection_results)) {
      if (running_compound) {
        num_eligible_by_comp = sapply(compound_variants$definitions,
                                      function(x) variants[x, min(num_covered_and_in_samples)], USE.NAMES = TRUE)
        selection_results[, held_out := num_eligible_by_comp[variant_name] - included_total]
      } else {
        selection_results[variants, held_out := num_covered_and_in_samples - included_total, on = 'variant_id']
      }
    }
  }

  if(all(c('held_out', 'included_total') %in% names(selection_results))) {
    selection_results[, uncovered := samples[, .N] - included_total - held_out]
  }

  if(return_fit) {
    fits = lapply(fits, function(x) {
      x@call.orig = call('[not shown]')
      parent.env(environment(x)) = emptyenv() # otherwise, whole CESAnalysis will be embedded via parent
      return(x)
    })
    setattr(selection_results, 'fit', fits)
  }
  curr_results = list(selection_results)
  names(curr_results) = run_name


  cesa@selection_results = c(cesa@selection_results, curr_results)
  return(cesa)
}




#' Estimate variant effects with a generalized sigmoid continuous-covariate model
#'
#' Estimates cancer effects for individual or compound variants when the
#' selection intensity varies with a continuous sample-level covariate,
#' according to a four-parameter generalized sigmoid curve:
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
#' Baseline mutation rates are obtained from the \code{CESAnalysis} object,
#' and model parameters are estimated separately for each variant by maximum
#' likelihood.
#'
#'
#' @inheritParams ces_variant
#'
#' @param model A likelihood-function factory defining the selection model.
#'   The package provides likelihood functions for continuous-covariate
#'   selection that can be used directly; for the generalized sigmoid model, use
#'   \code{sswm_age_lik_sigmoid()}. A custom likelihood-function factory may also
#'   be supplied. See \code{ces_variant()} for the general custom-model interface.
#'
#' @param lik_args Named list of additional arguments passed to \code{model}.
#'   For continuous-covariate inference, this must include
#'   \code{covariate_data}, a table associating
#'   \code{Unique_Patient_Identifier} values with a numeric continuous
#'   covariate in \code{covariate_value}.
#'
#' @param optimizer Character scalar specifying the optimization algorithm.
#'   Currently supported values are \code{"COBYLA"} and \code{"ISRES"},
#'   corresponding to the NLOPT algorithms \code{NLOPT_LN_COBYLA} and
#'   \code{NLOPT_GN_ISRES}, respectively.
#'
#' @param constraint Positive numeric value specifying the minimum permitted
#'   separation between the two asymptotes. The optimization imposes
#'   \code{exp(L) - exp(C) >= constraint}.
#'
#' @param conf Nominal confidence level. Set to \code{NULL} to skip confidence
#'   interval calculation. For inference on continuous-model parameters,
#'   confidence intervals can also be calculated separately using the
#'   continuous-model confidence-interval routines.
#'
#' @return A \code{CESAnalysis} object with a new entry appended to
#'   \code{selection_results}. For each analyzed variant, the result contains
#'   the fitted generalized sigmoid-model parameters ordered as \code{C}, \code{L},
#'   \code{s}, and \code{m} and maximized log-likelihood, together with
#'   variant annotations and applicable sample-count information.
#'
#' @seealso \code{\link{ces_variant}}
#'
#' @family continuous selection models
#'
#' @export
ces_variant_sigmoid <- function(cesa = NULL,
                                variants = select_variants(cesa, min_freq = 2),
                                samples = character(),
                                model = "default",
                                run_name = "auto",
                                lik_args = list(),
                                optimizer_args = if(identical(model, 'default')) list(method = 'L-BFGS-B', lower = 1e-3, upper = 1e9) else list(),
                                return_fit = FALSE,
                                hold_out_same_gene_samples = "auto",
                                cores = 1,
                                conf = .95,
                                optimizer = c("COBYLA", "ISRES"),
                                constraint = 1e-3)
{
  if(! is.numeric(cores) || length(cores) != 1 || cores - as.integer(cores) != 0 || cores < 1) {
    stop('cores should be 1-length positive integer')
  }

  if(! rlang::is_bool(return_fit)) {
    stop('return_fit should be TRUE/FALSE.')
  }

  validate_optimizer_args(optimizer_args)

  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (length(hold_out_same_gene_samples) == 1) {
    if (! is.logical(hold_out_same_gene_samples)) {
      if (identical(hold_out_same_gene_samples, "auto")) {
        hold_out_same_gene_samples = is(variants, "data.table")
      } else {
        stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
      }
    }
  } else {
    stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
  }

  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if(run_name %in% names(cesa@selection_results)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if (run_name == "auto") {
    # sequentially name results, handling nefarious run naming
    run_number = length(cesa@selection_results) + 1
    run_name = paste0('selection.', run_number)
    while(run_name %in% names(cesa@selection_results)) {
      run_number = run_number + 1
      run_name = paste0('variant_effects_', run_number)
    }
  }

  if(is(model, "character")) {
    # old names were basic, sswm-sequential (no one could remember if hyphen or underscore)
    model = tolower(model)
    model[model %in% c('sswm', 'default')] = 'basic'
    model[model %like% 'sswm[-_]sequential'] = 'sequential'
    if(length(model) != 1 || ! model %in% c("basic", "sequential")) {
      stop("model should specify a built-in selection model (i.e., \"default\") or a custom function factory.")
    } else {
      if (model == 'basic') {
        lik_factory = sswm_lik
      } else if(model == 'sequential') {
        lik_factory = sswm_sequential_lik
      } else {
        stop("Unrecognized model")
      }
    }
  } else if (! is(model, "function")) {
    stop("model should specify a built-in selection model (\"default\") or a custom function factory.")
  } else {
    lik_factory = model
  }

  if(! is(lik_args, "list")) {
    stop("lik args should be named list")
  }

  if(length(lik_args) != uniqueN(names(lik_args))) {
    stop('lik_args should be a named list without repeated names.')
  }

  optimizer <- match.arg(optimizer)

  samples = select_samples(cesa, samples)
  if(samples[, .N] < cesa@samples[, .N]) {
    num_excluded = cesa@samples[, .N] - samples[, .N]
    pretty_message(paste0("Note that ", format(num_excluded, big.mark = ','), " samples are being excluded from selection inference."))
  }

  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())

  # Set keys in case they've been lost
  mutations = cesa@mutations
  if(! is.null(conf)) {
    if(is(model, 'function')) {
      if(! rlang::is_scalar_double(conf) || conf != .95) {
        warning('conf is ignored when running a custom model.')
      }
      conf = NULL
    } else {
      if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
        stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
      }
    }
  }

  running_compound = FALSE

  # If an input variant table came directly from select_variants() and the variants are non-overlapping,
  # just accept the table. Otherwise, re-select the variants with the variant_id field.
  if (is(variants, "data.table")) {
    if(! "variant_id" %in% names(variants)) {
      stop("variants table is missing a variant_id column. Typically, variants is generated using select_variants().")
    }
    nonoverlapping = attr(variants, "nonoverlapping")
    if(is.null(nonoverlapping)) {
      if ('variant_id' %in% names(variants)) {
        pretty_message('Taking variants from variant_id column of input table....')
      }
    } else if(! identical(nonoverlapping, TRUE)) {
      stop("Input variants table may contain overlapping variants; re-run select_variants() to get a non-overlapping table.")
    }

    # re-select variants for maximum safety
    variants = select_variants(cesa, variant_ids = variants[, variant_id])

  } else if (is(variants, "CompoundVariantSet")) {
    running_compound = TRUE
    if (cesa@advanced$uid != variants@cesa_uid) {
      stop("Input CompoundVariantSet does not appear to derive from the input CESAnalysis.")
    }
    if (cesa@samples[, .N] != variants@cesa_num_samples) {
      stop("The number of samples in the CESAnalysis has changed since the CompoundVariantSet was created. ",
           "Please re-generate it.")
    }
    compound_variants = variants
    variants = select_variants(cesa, variant_ids = compound_variants@snvs$snv_id, include_subvariants = TRUE)

    # copy in compound variant names and overwrite covered_in with value of shared_cov
    variants = variants[compound_variants@snvs, compound_name := compound_name, on = c(variant_id = "snv_id")]
    if(variants[, .N] != compound_variants@snvs[, .N]) {
      stop("Internal error: select_variants() didn't return variant info 1-to-1 with compound input.")
    }
    variants[compound_variants@compounds, covered_in := shared_cov, on = "compound_name"]
  } else {
    stop("variants expected to be a variant table (from select_variants(), usually) or a CompoundVariantSet")
  }
  if (variants[, .N] == 0) {
    stop("There are no variants in the input!")
  }
  aac_ids = variants[variant_type == "aac", variant_id]

  # By noncoding, we just mean that SIs are calculated at the SNV site rather than at the AAC level,
  # regardless of whether there's a CDS annotation.
  noncoding_snv_ids = variants[variant_type == "snv", variant_id]

  if(length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", call. = F)
  }

  # identify mutations by nearest gene(s)
  maf = cesa@maf[samples$Unique_Patient_Identifier, on = "Unique_Patient_Identifier", nomatch = NULL]
  tmp = unique(maf[, .(gene = unlist(genes)), by = "Unique_Patient_Identifier"])[, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  tumors_with_variants_by_gene = list2env(tumors_with_variants_by_gene)

  # identify mutations by sample
  snv_aac_of_interest = cesa@mutations$aac_snv_key[aac_ids, on = 'aac_id']
  tmp = maf[snv_aac_of_interest, .(variant_id), on = c(variant_id = 'snv_id'),
            by = "Unique_Patient_Identifier", nomatch = NULL]
  tmp[snv_aac_of_interest, aac_id := aac_id, on = c(variant_id = 'snv_id')]
  tmp = tmp[, .(samples = list(unique(Unique_Patient_Identifier))), by = "aac_id"]
  samples_by_aac = setNames(tmp$samples, tmp$aac_id)

  setkey(maf, "variant_id")
  # need nomatch because some noncoding SNVs may not be present in the samples
  tmp = maf[noncoding_snv_ids, variant_id, by = "Unique_Patient_Identifier", nomatch = NULL][, .(samples = list(Unique_Patient_Identifier)), by = "variant_id"]
  samples_by_snv = tmp$samples
  names(samples_by_snv) = tmp$variant_id
  samples_by_variant = list2env(c(samples_by_aac, samples_by_snv))


  setkey(samples, "covered_regions")
  # These are WGS samples with purportedly whole-genome coverage.
  # That is, for better or worse, assuming that any variant can be found in these samples.
  # (Trimmed-interval WGS samples will have coverage = "genome" and covered_regions != "genome.")
  genome_wide_cov_samples = samples["genome", Unique_Patient_Identifier, nomatch = NULL]


  # Will process variants by coverage group (i.e., groups of variants that have the same tumors covering them)
  selection_results = NULL

  all_coverage = rbind(cesa@mutations$snv[, .(variant_id = snv_id, covered_in)],
                       cesa@mutations$amino_acid_change[, .(variant_id = aac_id, covered_in)])
  variants[all_coverage, covered_in := covered_in, on = 'variant_id']
  coverage_groups = unique(variants$covered_in)
  num_coverage_groups = length(coverage_groups)

  setkey(maf, "variant_id")
  setkey(variants, "variant_id")

  selection_results = lapply(1:length(coverage_groups), function(i) {
    coverage_group = coverage_groups[[i]]
    # one coverage group may be NA, for variants that are not covered by any specific covered_regions
    if(length(coverage_group) == 1 && is.na(coverage_group)) {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, NA_character_)))]
    } else {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, coverage_group)))]
    }
    message(sprintf("Preparing to calculate cancer effects (batch %i of %i)...", i, num_coverage_groups))
    if(is.null(coverage_group)) {
      covered_samples = genome_wide_cov_samples
    } else {
      covered_samples = c(samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL], genome_wide_cov_samples)
    }
    variants[curr_variants$variant_id, num_covered_and_in_samples := length(covered_samples), on = 'variant_id']

    # When not all samples are used, it's possible that no sampels will have coverage at the input variants.
    if(length(covered_samples) == 0) {
      # message("Skipped batch ", i, " because no samples had coverage at the variant sites in the batch.")
      return(list(data.table(), NULL))
    }
    # rough size of baseline rates data.table in bytes, if all included in one table
    work_size = length(covered_samples) * curr_variants[,.N] * 8

    # we divide into subgroups to cap baseline rates table at around 1 GB
    num_proc_groups = ceiling(work_size / 1e9)
    curr_variants[, subgroup := ceiling(num_proc_groups * 1:.N / .N)]

    if (running_compound) {
      # can't have subvariants of the same compound variant ending up in different subgroups
      curr_variants[, subgroup := rep.int(subgroup[1], .N), by = "compound_name"]
      num_proc_groups = max(curr_variants$subgroup) # rarely, last subgroup dropped by above
    }


    curr_results = lapply(1:num_proc_groups, function(j) {
      if (num_proc_groups > 1) {
        message(sprintf("Working on sub-batch %i of %i...", j, num_proc_groups))
      }
      curr_subgroup = curr_variants[subgroup == j]
      aac_ids = curr_subgroup[variant_type == "aac", variant_id]
      snv_ids = curr_subgroup[variant_type == "snv", variant_id]

      baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, snv_ids = snv_ids, samples = covered_samples)
      # put gene(s) by variant into env for quick access
      gene_lookup = curr_subgroup[, all_genes]
      names(gene_lookup) = curr_subgroup[, variant_id]
      gene_lookup = list2env(gene_lookup)

      # function to run MLE on given variant_id (vector of IDs for compound variants)
      process_variant = function(variant_id) {
        if(running_compound) {
          compound_id = variant_id
          # Important: sample_calls includes samples that are not in shared coverage
          tumors_with_variant = intersect(compound_variants@sample_calls[[compound_id]], covered_samples)
          current_snvs = compound_variants@snvs[compound_name == compound_id]
          all_genes = current_snvs[, unique(unlist(genes))]
          variant_id = current_snvs$snv_id
          rates = baseline_rates[, ..variant_id]
          # Sum Poisson rates across variants
          rates = rowSums(rates)
        } else {
          tumors_with_variant = samples_by_variant[[variant_id]]
          all_genes = gene_lookup[[variant_id]]
          rates = baseline_rates[, ..variant_id][[1]]
        }
        names(rates) = baseline_rates[, Unique_Patient_Identifier]

        # usually but not always just 1 gene when not compound (when compound, anything possible)
        if (hold_out_same_gene_samples) {
          if (length(all_genes) == 1) {
            if(is.na(all_genes)) {
              tumors_with_gene_mutated = tumors_with_variant
            } else {
              tumors_with_gene_mutated = tumors_with_variants_by_gene[[all_genes]]
            }
          } else {
            tumors_with_gene_mutated = unique(unlist(sapply(all_genes, function(x) tumors_with_variants_by_gene[[x]])))
          }
          tumors_without = setdiff(covered_samples, tumors_with_gene_mutated)
        } else {
          tumors_without = setdiff(covered_samples, tumors_with_variant)
        }
        rates_tumors_with = rates[tumors_with_variant]
        rates_tumors_without = rates[tumors_without]


        lik_args = c(list(rates_tumors_with = rates_tumors_with, rates_tumors_without = rates_tumors_without),
                     lik_args)
        fn = do.call(lik_factory, lik_args)

        # Generalized sigmoid model, constrain the upper asymptote to exceed the lower asymptote
        eval_g0 <- function(params){
          C <- params[1]
          L <- params[2]
          s <- params[3]
          m <- params[4]
          return(C - L + constraint)
        }

        # Initial parameters
        par_init = formals(fn)[[1]]
        names(par_init) = bbmle::parnames(fn)

        if (optimizer == "COBYLA") {
          print("Method used: nloptr COBYLA")
          fit <- nloptr(x0=par_init,
                        eval_f=fn,
                        lb = c(-100,-100,-100,1e-6),
                        ub = c(100,100,100,1e3),
                        eval_g_ineq = eval_g0,
                        opts = list("algorithm"="NLOPT_LN_COBYLA",
                                    "xtol_rel"=1.0e-8,
                                    "maxeval"=10000))

          selection_intensity <- fit$solution
          loglikelihood <- -fit$objective
          print(paste0("loglik: ", loglikelihood))

        } else if (optimizer == "ISRES"){
          print("Method used: nloptr ISRES")
          fit <- nloptr(x0=par_init,
                        eval_f=fn,
                        lb = c(-100,-100,-100,1e-6),
                        ub = c(100,100,100,1e3),
                        eval_g_ineq = eval_g0,
                        opts = list("algorithm"="NLOPT_GN_ISRES",
                                    "xtol_rel"=1.0e-9,
                                    "maxeval"=50000))

          selection_intensity <- fit$solution
          loglikelihood <- -fit$objective
          print(paste0("loglik: ", loglikelihood))
        }

        names(selection_intensity) <- names(par_init)

        if (running_compound) {
          variant_id = compound_id
        }
        variant_output = c(list(variant_id = variant_id),
                           as.list(selection_intensity),
                           list(loglikelihood = loglikelihood))

        if(is.character(model) || is.null(lik_args$covariate_data)){

          if(is.character(model)){
            if (model == 'basic') {
              # Record counts of total samples included in inference and included samples with the variant.
              # This may vary from the naive output of variant_counts() due to issues of sample coverage and
              # (by default) the use of hold_out_same_gene_samples = TRUE.

              num_samples_with = length(tumors_with_variant)
              num_samples_total = num_samples_with + length(tumors_without)
              variant_output = c(variant_output, list(included_with_variant = num_samples_with,
                                                      included_total = num_samples_total))
            }
          }

          if(is(model, "function") && is.null(lik_args$covariate_data)){
            num_samples_with = length(tumors_with_variant)
            num_samples_total = num_samples_with + length(tumors_without)
            variant_output = c(variant_output, list(included_with_variant = num_samples_with,
                                                    included_total = num_samples_total))
          }
        }
        if(! is.null(conf)) {
          min_value = ifelse(is.null(final_optimizer_args$lower), -Inf, final_optimizer_args$lower)
          max_value = ifelse(is.null(final_optimizer_args$upper), Inf, final_optimizer_args$upper)
          variant_output = c(variant_output,
                             univariate_si_conf_ints(fit, fn, min_value, max_value, conf))
        }
        return(list(summary = variant_output, fit = if (return_fit) fit else NULL))
      }

      if (running_compound) {
        variants_to_run = curr_subgroup[, unique(compound_name)]
      } else {
        variants_to_run = curr_subgroup$variant_id
      }
      message("Calculating cancer effects...")
      subgroup_results = pbapply::pblapply(variants_to_run, process_variant, cl = cores)
      subgroup_selection = rbindlist(lapply(subgroup_results, '[[', 1))
      subgroup_fit = lapply(subgroup_results, '[[', 2)
      return(list(subgroup_selection, subgroup_fit))
    })
    group_selection = rbindlist(lapply(curr_results, '[[', 1))
    group_fit = lapply(curr_results, '[[', 2)
    return(list(group_selection, group_fit))
  })

  fits = unlist(lapply(selection_results, '[[', 2))
  selection_results = rbindlist(lapply(selection_results, '[[', 1))

  if(selection_results[, .N] == 0) {
    msg = paste0("No selection inference was performed, so returning CESAnalysis unaltered without selection output. ",
                 "Perhaps none of the variants had coverage in the specified samples?")
    pretty_message(msg)
    return(cesa)
  }

  if (running_compound) {
    selection_results[, variant_type := "compound"]
    setnames(selection_results, "variant_id", "variant_name")
    selection_results = selection_results[compound_variants@compounds, on = c(variant_name = "compound_name")]
    setattr(selection_results, "is_compound", TRUE)
    setcolorder(selection_results, c("variant_name", "variant_type"))
  } else {
    # Fill in top-priority gene for variants that are in coding regions or essential splice (or within 1 bp)
    # Other variants will get NA gene.
    selection_results[variants, c("variant_type", "variant_name", "gene", "intergenic") :=
                        list(variant_type, variant_name, gene, intergenic), on = "variant_id"]
    selection_results[intergenic == T, gene := NA]
    selection_results$intergenic = NULL
    setattr(selection_results, "is_compound", FALSE)
    setcolorder(selection_results, c("variant_name", "variant_type", "gene"))
    setcolorder(selection_results, c(setdiff(names(selection_results), 'variant_id'), 'variant_id'))
  }


  if(hold_out_same_gene_samples == FALSE) {
    selection_results[, held_out := 0]
  } else {
    if('included_total' %in% names(selection_results)) {
      if (running_compound) {
        num_eligible_by_comp = sapply(compound_variants$definitions,
                                      function(x) variants[x, min(num_covered_and_in_samples)], USE.NAMES = TRUE)
        selection_results[, held_out := num_eligible_by_comp[variant_name] - included_total]
      } else {
        selection_results[variants, held_out := num_covered_and_in_samples - included_total, on = 'variant_id']
      }
    }
  }

  if(all(c('held_out', 'included_total') %in% names(selection_results))) {
    selection_results[, uncovered := samples[, .N] - included_total - held_out]
  }

  if(return_fit) {
    fits = lapply(fits, function(x) {
      x@call.orig = call('[not shown]')
      parent.env(environment(x)) = emptyenv() # otherwise, whole CESAnalysis will be embedded via parent
      return(x)
    })
    setattr(selection_results, 'fit', fits)
  }
  curr_results = list(selection_results)
  names(curr_results) = run_name


  cesa@selection_results = c(cesa@selection_results, curr_results)
  return(cesa)
}
