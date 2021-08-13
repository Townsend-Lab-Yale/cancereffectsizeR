#' Calculate variant effect size
#'
#' This function calculates variant effect sizes under the chosen model of selection. By
#' default, a variant is assumed to have a consistent selection intensity across all
#' samples. Set `model = "sswm_sequential"` to allow selection intensity to vary among
#' sequential sample groups (e.g., stages 1-4; local/distant metastases). Use `groups` to
#' define group ordering or to restrict which groups are considered under either built-in
#' model. By default, only variants with MAF frequency > 1 (i.e., recurrent variants) are
#' tested. To include all variants, or to otherwise customize which variants to include,
#' call select_variants() with desired parameters.
#' 
#' It's possible to pass in your own selection model. You'll need to create a "function
#' factory" that, for any variant, produces a likelihood function that can be evaluated on
#' the data. The first two arguments must be rates_tumors_with and rates_tumors_without,
#' which take the baseline site mutation rates in samples with and without the variant.
#' The third argument must be sample_index, which associates Unique_Patient_Identifiers
#' with their sample groups. Values for all three of these arguments will be calculated by
#' ces_variant and passed to your function factory automatically. Your function can take
#' whatever additional arguments you like, and you can pass in values using
#' `custom_lik_args`. The likelihood function parameters that ces_variant will optimize
#' should be named and have default values. See the source code of `sswm_sequential_lik()`
#' for an example.
#' 
#' 
#' @param cesa CESAnalysis object
#' @param cores number of cores to use
#' @param conf selection intensity confidence interval width (NULL skips calculation,
#'   speeds runtime)
#' @param run_name Optionally, a name to identify the current run.
#' @param variants Variant table from \code{select_variants()}, or a \code{CompoundVariantSet} from
#'   \code{define_compound_variants()}. Defaults to all recurrent noncoding SNVs and
#'   (SNV-derived) coding mutations, where recurrent means appearing in at least two
#'   samples in the MAF data set.
#' @param model "sswm" or "sswm_sequential" to use built-in models of selection, or supply a
#' custom function factory (see details).
#' @param groups Which sample groups to include (default all); data for outgroup samples
#'   will not inform selection calculation. For models (like sswm_sequential) that assume
#'   ordered groups of samples, use a list to indicate group ordering.
#'   Examples: \code{c("group1", "group2")} includes groups 1 and 2, but doesn't indicate
#'   an ordering, so is invalid for ordered-group models.
#'   \code{list("Primary", "Metastatic")} indicates two ordered groups, and
#'   \code{list("A", c("B", "C"), "D")} means that group A is first, groups B and C share
#'   an intermediate state, and group D is last.
#' @param custom_lik_args Extra arguments, given as a list, to pass to custom likelihood
#'   functions.
#' @param hold_out_same_gene_samples When finding likelihood of each variant, hold out
#'   samples that lack the variant but have any other mutations in the same gene. By default,
#'   TRUE when running with single variants, FALSE with a CompoundVariantSet.
#' @return CESAnalysis object with selection results appended to the selection output list
#' @export

ces_variant <- function(cesa = NULL,
                        variants = select_variants(cesa, min_freq = 2),
                        run_name = "auto",
                        model = "sswm",
                        groups = NULL,
                        custom_lik_args = list(),
                        hold_out_same_gene_samples = "auto",
                        cores = 1,
                        conf = .95) 
{
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis", call. = F)
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
    if(length(model) != 1 || ! model %in% c("sswm", "sswm_sequential")) {
      stop("model should specify a built-in selection model (sswm, sswm_sequential) or a custom function factory.")
    } else {
      if (model == "sswm") {
        lik_factory = sswm_lik
      } else if(model == "sswm_sequential") {
        if (is.null(groups)) {
          stop("an ordered list of groups must be supplied with the \"groups\" option to use the sswm_sequential model")
        }
        lik_factory = sswm_sequential_lik
      } else {
        stop("Unrecognized model")
      }
    }
  } else if (! is(model, "function")) {
    stop("model should specify a built-in selection model (sswm, sswm_sequential) or a custom function factory.")
  } else {
    lik_factory = model
  }
  
  
  sample_index = numeric() # named vector giving position of each tumor in sequential ordering (unused in some models)
  group_ordering = groups
  if (is.null(group_ordering)) {
    group_ordering = cesa@groups
  }

  running_builtin_sswm = FALSE
  if (identical(lik_factory, sswm_lik)) {
    running_builtin_sswm = TRUE
    if (! is(group_ordering, "character")) {
      stop("\"groups\" should be supplied as a character vector for the sswm model.")
    }
  }
  
  if(is(group_ordering, "character")) {
    group_ordering = as.list(group_ordering)
  }
  if(! is(group_ordering, "list")) {
    stop("group_ordering should be character vector or list")
  }
  if (length(group_ordering) < 2 && identical(lik_factory, sswm_sequential_lik)) {
    stop("\"groups\" should be a list with length at least two for the chosen model")
  }
  
  used_groups = character()
  for (i in 1:length(group_ordering)) {
    curr_state = group_ordering[[i]]
    if (! is(curr_state, "character")) {
      stop("Each element of group_ordering should be type character")
    }
    if (length(curr_state) != length(unique(curr_state))) {
      stop("Double-check your group_ordering")
    }
    if (! all(curr_state %in% cesa@groups)) {
      invalid_groups = setdiff(curr_state, cesa@groups)
      stop("Group not declared in CESAnalysis: ", paste0(invalid_groups, collapse = ", "), ".")
    }
    curr_samples = cesa@samples[group %in% curr_state, Unique_Patient_Identifier]
    if (length(curr_samples) == 0) {
      stop("Some groupings given by group_ordering have no associated samples.")
    }
    if (any(curr_state %in% used_groups)) {
      stop("CESAnalysis groups are re-used in group_ordering")
    }
    used_groups = c(curr_state, used_groups)
    curr_grouping = rep(i, length(curr_samples))
    names(curr_grouping) = curr_samples
    sample_index = c(sample_index, curr_grouping)
  }
  cesa = copy_cesa(cesa)
  unused_groups = setdiff(cesa@groups, used_groups)
  if (length(unused_groups) > 0) {
    pretty_message(paste0("The following CESAnalysis groups were not included in \"groups\", so they are not informing effect size:\n",
            paste(unused_groups, collapse = ", "), "."), black = F)
    samples = cesa@samples[used_groups, on = "group"]
    maf = cesa@maf[samples$Unique_Patient_Identifier, on = "Unique_Patient_Identifier"]
  } else {
    samples = cesa@samples
    maf = cesa@maf
  }
  
  cesa = update_cesa_history(cesa, match.call())
  
  # Set keys in case they've been lost
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  mutations = cesa@mutations
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }
  
  running_compound = FALSE
  if (is(variants, "data.table")) {
    nonoverlapping = attr(variants, "nonoverlapping")
    if(is.null(nonoverlapping)) {
      stop("Input variants table is missing attribute nonoverlapping (probably, it wasn't generated by select_variants()")
    } else if(! identical(nonoverlapping, TRUE)) {
      stop("Input variants table may contain overlapping variants; re-run select_variants() to get a non-overlapping table.")
    }
    
    # use input table without calling select_variants again if it came from the CESAnalysis
    if (! identical(attr(variants, "cesa_id"), cesa@advanced$uid)) {
      if(! "variant_id" %in% names(variants)) {
        stop("variants is missing a variant_id column. Typically, variants is generated using select_variants().")
      }
      variants = select_variants(cesa, variant_passlist = variants[, variant_id])
    }
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
    variants = select_variants(cesa, variant_passlist = compound_variants@snvs$snv_id, include_subvariants = TRUE)
    
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
  # regardless of whether there's a CDS annotation. We're not going to check to see if the user's
  # AACs and SNVs overlap (that's their problem if they decided to supply overlapping variants)
  noncoding_snv_ids = variants[variant_type == "snv", variant_id]
  
  if(length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", call. = F)
  }

  # identify mutations by nearest gene(s)
  tmp = unique(maf[, .(gene = unlist(genes)), by = "Unique_Patient_Identifier"])[, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  tumors_with_variants_by_gene = list2env(tumors_with_variants_by_gene)
  
  # identify mutations by sample
  tmp = maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "aac_id"]
  tmp = tmp[aac_ids, on = "aac_id"] # subset to coding mutations in use
  aacs_by_tumor = tmp$samples
  names(aacs_by_tumor) = tmp$aac_id
  
  setkey(maf, "variant_id")
  # need nomatch because some noncoding SNVs may not be present in the samples
  tmp = maf[noncoding_snv_ids, variant_id, by = "Unique_Patient_Identifier", nomatch = NULL][, .(samples = list(Unique_Patient_Identifier)), by = "variant_id"]
  snvs_by_tumor = tmp$samples
  names(snvs_by_tumor) = tmp$variant_id
  variants_by_tumor = list2env(c(aacs_by_tumor, snvs_by_tumor))
  
  
  setkey(samples, "covered_regions")
  # These are WGS samples with purportedly whole-genome coverage.
  # That is, for better or worse, assuming that any variant can be found in these samples.
  # (Trimmed-interval WGS samples will have coverage = "genome" and covered_regions != "genome.")
  genome_wide_cov_samples = samples["genome", Unique_Patient_Identifier, nomatch = NULL]
  

  # Will process variants by coverage group (i.e., groups of variants that have the same tumors covering them)
  selection_results = NULL
  coverage_groups = unique(variants$covered_in)
  num_coverage_groups = length(coverage_groups)
  
  i = 0
  setkey(maf, "variant_id")
  setkey(variants, "variant_id")
  for (coverage_group in coverage_groups) {
    # one coverage group may be NA, for variants that are not covered by any specific covered_regions
    if(length(coverage_group) == 1 && is.na(coverage_group)) {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, NA_character_)))]
    } else {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, coverage_group)))]
    }
    i = i+1
    message(sprintf("Preparing to calculate selection intensities (batch %i of %i)...", i, num_coverage_groups))
    covered_samples = c(samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL], genome_wide_cov_samples)
    
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
    
    
    for (j in 1:num_proc_groups) {
      if (num_proc_groups > 1) {
        message(sprintf("Working on sub-batch %i of %i...", j, num_proc_groups))
      }
      curr_subgroup = curr_variants[subgroup == j]
      aac_ids = curr_subgroup[variant_type == "aac", variant_id]
      snv_ids = curr_subgroup[variant_type == "snv", variant_id]
      
      baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, snv_ids = snv_ids, samples = covered_samples, cores = cores)
      
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
          # binomial probability of having 1 or more of the mutations (assuming independence)
          rates = apply(rates, 1, function(x) 1 - prod(1 - x))
        } else {
          tumors_with_variant = variants_by_tumor[[variant_id]]
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
            tumors_with_gene_mutated = unique(sapply(all_genes, function(x) tumors_with_variants_by_gene[[x]]))
          }
          tumors_without = setdiff(covered_samples, tumors_with_gene_mutated)
        } else {
          tumors_without = setdiff(covered_samples, tumors_with_variant)
        }
        rates_tumors_with = rates[tumors_with_variant]
        rates_tumors_without = rates[tumors_without]
        
        lik_args = c(list(rates_tumors_with = rates_tumors_with, rates_tumors_without = rates_tumors_without), 
                     custom_lik_args)
        
        # sample_index supplied if defined and not running our SSWM (it doensn't use it)
        if(! running_builtin_sswm) {
          lik_args = c(lik_args, list(sample_index = sample_index))
        }
        fn = do.call(lik_factory, lik_args)
        
        par_init = formals(fn)[[1]]
        names(par_init) = bbmle::parnames(fn)
        

        # find optimized selection intensities
        # the selection intensity for any stage that has 0 variants will be on the lower boundary; will muffle the associated warning
        withCallingHandlers(
          {
            fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = T, lower=1e-3, upper=1e9, control=list(fnscale=1e-12))
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
        
        #dndscv_q = sapply(cesa@dndscv_out_list, function(x) x$sel_cv[x$sel_cv$gene_name == mut_record$gene, "qallsubs_cv"])
        if (running_compound) {
          variant_id = compound_id
        }
        variant_output = c(list(variant_id = variant_id), 
                           as.list(selection_intensity),
                           list(loglikelihood = loglikelihood))
        
        if(! is.null(conf)) {
          variant_output = c(variant_output, univariate_si_conf_ints(fit, fn, .001, 1e20, conf))
        }
        return(variant_output)
      }
      
      if (running_compound) {
        variants_to_run = curr_subgroup[, unique(compound_name)]
      } else {
        variants_to_run = curr_subgroup$variant_id
      }
      message("Calculating SIs...")
      selection_results = rbind(selection_results, rbindlist(pbapply::pblapply(variants_to_run, process_variant, cl = cores)))
    }
  }
  
  if (running_compound) {
    selection_results[, variant_type := "compound"]
    setnames(selection_results, "variant_id", "variant_name")
    selection_results = selection_results[compound_variants@compounds, on = c(variant_name = "compound_name")]
    setattr(selection_results, "is_compound", TRUE)
    setcolorder(selection_results, c("variant_name", "variant_type"))
  } else {
    selection_results[variants, variant_type := variant_type, on = "variant_id"]
    setattr(selection_results, "is_compound", FALSE)
    setcolorder(selection_results, c("variant_id", "variant_type"))
  }
  
  # isolate SI columns for plotting functions (and maybe more, eventually)
  ll_col_num = which(colnames(selection_results) == "loglikelihood")
  si_cols = colnames(selection_results)[3:(ll_col_num - 1)]
  setattr(selection_results, "si_cols", si_cols)
  
  curr_results = list(selection_results)
  names(curr_results) = run_name
  cesa@selection_results = c(cesa@selection_results, curr_results)
  
  return(cesa)
}






