#' Calculate cancer effects of variants
#'
#' This function calculates variant effect sizes under the chosen model of selection. Under the
#' default model, a variant is assumed to have a consistent scaled selection coefficient (cancer
#' effect) across all included samples.
#' 
#' Definitions of the sample count columns in the effects output:
#' \itemize{
#'   \item included_with_variant: Number of samples that have the variant and were included in the inference.
#'   \item included_total: Number of samples that have coverage at the site and were included in the inference.
#'   \item held_out: Samples that have coverage at the site, but were held out of the inference due to \code{hold_out_same_gene_samples = TRUE}.
#'   \item uncovered: Samples that were not included in the inference because their sequencing did not cover the variant site.
#'  }
#' Note that if a table of samples to include in the inference is specified with \code{samples}, any
#' CESAnalysis samples not present in the table will not be included in any of the above accounts.
#' 
#' It's possible to pass in your own selection model. You'll need to create a "function factory"
#' that, for any variant, produces a likelihood function that can be evaluated on the data. The
#' first two arguments must be \code{rates_tumors_with} and \code{rates_tumors_without}, which give the baseline
#' site mutation rates in samples with and without the variant. The third argument must be
#' \code{sample_index}, a data.table that associates \code{Unique_Patient_Identifier} with group names and
#' indices. (Your function factory must accept this argument, but it doesn't have to use its value.)
#' Values for all three of these arguments will be calculated by ces_variant and passed to your
#' function factory automatically. Your function can take whatever additional arguments you like,
#' and you can pass in values using \code{lik_args}. The likelihood function parameters that
#' ces_variant will optimize should be named and have default values. See the source code of
#' \code{sswm_lik()} for an example.
#' 
#' 
#' @param cesa CESAnalysis object
#' @param cores Number of cores to use for processing variants in parallel (not useful for Windows
#'   systems).
#' @param conf Cancer effect confidence interval width (NULL skips calculation, speeds runtime). Ignored
#' when running custom models.
#' @param samples Which samples to include in inference. Defaults to all samples. Can be a vector of
#'   Unique_Patient_Identifiers, or a data.table containing rows from the CESAnalysis sample table.
#' @param variants Which variants to estimate effects for, specified with a variant table such as
#'   from \code{[CESAnalysis]$variants} or \code{select_variants()}, or a \code{CompoundVariantSet}
#'   from \code{define_compound_variants()}. Defaults to all recurrent mutations; that is,
#'   \code{[CESAnalysis]$variants[maf_prevalence > 1]}. To include all variants, set to
#'   \code{[CESAnalysis]$variants}.
#' @param model Set to "basic" (default) or "sequential" (not yet available) to use built-in
#'   models of selection, or supply a custom function factory (see details).
#' @param run_name Optionally, a name to identify the current run.
#' @param lik_args Extra arguments, given as a list, to pass to custom likelihood functions.
#' @param optimizer_args Named list of arguments to pass to the optimizer, bbmle::mle2. Use, for example,
#' to choose optimization algorithm or parameter boundaries on custom models.
#' @param return_fit TRUE/FALSE (default FALSE): Embed model fit for each variant in a "fit"
#'   attribute of the selection results table. Use \code{attr(selection_table, 'fit')} to access the
#'   list of fitted models. Defaults to FALSE to save memory. Model fit objects can be of moderate
#'   or large size. If you run thousands of variants at once, you may exhaust your system memory.
#' @param hold_out_same_gene_samples When finding likelihood of each variant, hold out samples that
#'   lack the variant but have any other mutations in the same gene. By default, TRUE when running
#'   with single variants, FALSE with a CompoundVariantSet.
#' @return CESAnalysis object with selection results appended to the selection output list
#' @export
ces_variant <- function(cesa = NULL,
                        variants = select_variants(cesa, min_freq = 2),
                        samples = character(),
                        model = "default",
                        run_name = "auto",
                        lik_args = list(),
                        optimizer_args = if(identical(model, 'default')) list(method = 'L-BFGS-B', lower = 1e-3, upper = 1e9) else list(),
                        return_fit = FALSE,
                        hold_out_same_gene_samples = "auto",
                        cores = 1,
                        conf = .95)
{
  if(! is.numeric(cores) || length(cores) != 1 || cores - as.integer(cores) != 0 || cores < 1) {
    stop('cores should be 1-length positive integer')
  }
  
  if(! rlang::is_bool(return_fit)) {
    stop('return_fit should be TRUE/FALSE.')
  }
  
  
  if (! is(optimizer_args, "list") || uniqueN(names(optimizer_args)) != length(optimizer_args)) {
    stop("optimizer_args should a named list of arguments to pass.")
  }
  
  reserved_args = c('minuslogl', 'start', 'vecpar')
  if(any(reserved_args %in% names(optimizer_args))) {
    msg = paste0('Optimizer arguments start, vecpar, and minuslogl cannot be changed here. ', 
                 'If you are using a custom model, your likelihood function can declare these ',
                 'values directly (see docs).')
    stop(pretty_message(msg, emit = F))
  }
  
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
  

  
  samples = select_samples(cesa, samples)
  if(samples[, .N] < cesa@samples[, .N]) {
    num_excluded = cesa@samples[, .N] - samples[, .N]
    pretty_message(paste0("Note that ", num_excluded, " samples are being excluded from selection inference."))
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
  
  i = 0
  setkey(maf, "variant_id")
  setkey(variants, "variant_id")
  
  selection_results = lapply(coverage_groups, function(coverage_group) {
    # one coverage group may be NA, for variants that are not covered by any specific covered_regions
    if(length(coverage_group) == 1 && is.na(coverage_group)) {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, NA_character_)))]
    } else {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, coverage_group)))]
    }
    i = i+1
    message(sprintf("Preparing to calculate cancer effects (batch %i of %i)...", i, num_coverage_groups))
    if(is.null(coverage_group)) {
      covered_samples = genome_wide_cov_samples
    } else {
      covered_samples = c(samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL], genome_wide_cov_samples)
    }
    variants[curr_variants$variant_id, num_covered_and_in_samples := length(covered_samples), on = 'variant_id']
    
    if(length(covered_samples) == 0) {
      warning("Skipped batch ", i, " because no samples had coverage at the variant sites in the batch.")
      next
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
        
        if (running_compound) {
          variant_id = compound_id
        }
        variant_output = c(list(variant_id = variant_id), 
                           as.list(selection_intensity),
                           list(loglikelihood = loglikelihood))
        
        if(is.character(model) || is.null(lik_args$sample_index)){
          
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
          
          if(is(model, "function") && is.null(lik_args$sample_index)){
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

#' Clear variant effect output
#' 
#' Remove output from previous ces_variant() runs from CESAnalysis
#' @param cesa CESAnalysis
#' @param run_names Which previous runs to remove; defaults to removing all.
#' @export
clear_effect_output = function(cesa, run_names = names(cesa$selection)) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  cesa = copy_cesa(cesa)
  
  if(length(cesa@selection_results) == 0) {
    warning("No ces_variant() runs to clear, so returning CESAnalysis unaltered.")
    return(cesa)
  }
  if(! is.character(run_names) || length(run_names) < 1) {
    stop("run_name should be type character.")
  }
  
  if(! all(run_names %in% names(cesa@selection_results))) {
    stop("Input run_names are not all present in ces_variant() output.")
  }
  cesa@selection_results[run_names] = NULL
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}

#' Clear epistasis output
#' 
#' Remove previous epistatic effect estimations from CESAnalysis.
#' @param cesa CESAnalysis.
#' @param run_names Which previous runs to remove; defaults to removing all.
#' @export
clear_epistasis_output = function(cesa, run_names = names(cesa$epistasis)) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  cesa = copy_cesa(cesa)
  
  if(! is.character(run_names) || length(run_names) < 1) {
    stop("run_name should be type character.")
  }
  
  if(! all(run_names %in% names(cesa@epistasis))) {
    stop("Input run_names are not all present in epistasis output.")
  }
  cesa@epistasis[run_names] = NULL
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}




