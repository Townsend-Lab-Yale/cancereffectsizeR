#' Calculate selection intensity for single-nucleotide variants and amino acid changes
#' 
#' @param cesa CESAnalysis object
#' @param cores number of cores to use
#' @param conf selection intensity confidence interval width (NULL skips calculation,
#'   speeds runtime)
#' @param variants Variant table from select_variants(), or a CompoundVariantSet from
#'   define_compound_variants(). Defaults to all recurrent noncoding SNVs and
#'   (SNV-derived) coding mutations, where recurrent means appearing in at least two
#'   samples in the MAF data set.
#' @param lik_fn "sswm" or "sswm_sequential" to use built-in models of selection, or supply a
#' custom function factory (see details).
#' @param group_ordering for models (like sswm_sequential) that assume different groups of
#'   samples are from sequential tumor progression states, a vector or list giving the
#'   ordering of CESAnalysis groups. Examples: c("Primary", "Metastatic") is a
#'   two-state ordering. Lists allow groups to share progression states: 
#'   list("A", c("B", "C"), "D") means that group A is first, groups B and C 
#'   belong to the same intermediate state, and group D is last.
#' @param custom_lik_args Extra arguments, given as a list, to pass to custom likelihood
#'   functions.
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export

ces_variant <- function(cesa = NULL,
                        variants = select_variants(cesa, min_freq = 2),
                        lik_fn = "sswm",
                        group_ordering = NULL,
                        custom_lik_args = list(),
                        cores = 1,
                        conf = .95) 
{
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis", call. = F)
  }
  
  if(is(lik_fn, "character")) {
    if(length(lik_fn) != 1 || ! lik_fn %in% c("sswm", "sswm_sequential")) {
      stop("lik_fn should specify a built-in selection model (sswm, sswm_sequential) or a custom function factory.")
    } else {
      if (lik_fn == "sswm") {
        lik_factory = sswm_lik
      } else if(lik_fn == "sswm_sequential") {
        if (is.null(group_ordering)) {
          stop("group_ordering must be supplied with the sswm_sequential model")
        }
        lik_factory = sswm_sequential_lik
      } else {
        stop("Unrecognized lik_fn")
      }
    }
  } else if (! is(lik_fn, "function")) {
    stop("lik_fn should specify a built-in selection model (sswm, sswm_sequential) or a custom function factory.")
  } else {
    lik_factory = lik_fn
  }
  
  sample_index = numeric() # named vector giving position of each tumor in sequential ordering (unused in some models)
  if (! is.null(group_ordering)) {
    if (length(cesa@groups) == 1) {
      stop("Your CESAnalysis does not have defined sample groups (see ?CESAnalysis).")
    }
    if (identical(lik_factory, sswm_lik)) {
      stop("group_ordering is not used by the sswm model (it assumes one selection intensity per variant).")
    }
    
    if(is(group_ordering, "character")) {
      group_ordering = as.list(group_ordering)
    }
    if(! is(group_ordering, "list")) {
      stop("group_ordering should be character vector or list")
    }
    if (length(group_ordering) < 2) {
      stop("group_ordering must have length of at least 2")
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
        stop("Not all group_ordering groups declared in the CESAnalysis.")
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
    unused_groups = setdiff(cesa@groups, used_groups)
    if (length(unused_groups) > 0) {
      warning("The following groups were not included in group_ordering, so they are not informing effect size:\n",
              paste(unused_groups, collapse = ", "), ".")
    }
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
  
  if (is(variants, "data.table")) {
    # use input table without calling select_variants again if it came from the CESAnalysis
    if (! identical(attr(variants, "cesa_id"), cesa@advanced$uid)) {
      if(! "variant_id" %in% names(variants)) {
        stop("variants is missing a variant_id column. Typically, variants is generated using select_variants().")
      }
      variants = select_variants(cesa, variant_passlist = variants[, variant_id])
    }
  } else if (is(variants, "CompoundVariantSet")) {
    compound_variants = variants
    variants = select_variants(cesa, variant_passlist = unlist(compound_variants@snv_id)) 
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
  tmp = unique(cesa@maf[variant_type == "snv", .(gene = unlist(genes)), by = "Unique_Patient_Identifier"])[, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  tumors_with_variants_by_gene = list2env(tumors_with_variants_by_gene)
  
  # identify mutations by sample
  tmp = cesa@maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "aac_id"]
  tmp = tmp[aac_ids, on = "aac_id"] # subset to coding mutations in use
  aacs_by_tumor = tmp$samples
  names(aacs_by_tumor) = tmp$aac_id
  
  setkey(cesa@maf, "variant_id")
  tmp = cesa@maf[noncoding_snv_ids, variant_id, by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "variant_id"]
  snvs_by_tumor = tmp$samples
  names(snvs_by_tumor) = tmp$variant_id
  variants_by_tumor = list2env(c(aacs_by_tumor, snvs_by_tumor))
  
  
  setkey(cesa@samples, "covered_regions")
  
  # these are WGS samples without any coverage limitation
  # that is, for better or worse, assuming that any variant can be found in these samples
  genome_wide_cov_samples = cesa@samples["genome", Unique_Patient_Identifier, nomatch = NULL]
  

  # Will process variants by coverage group (i.e., groups of variants that have the same tumors covering them)
  selection_results = NULL
  coverage_groups = unique(variants$covered_in)
  num_coverage_groups = length(coverage_groups)
  
  i = 0
  setkey(cesa@maf, "variant_id")
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
    covered_samples = c(cesa@samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL], genome_wide_cov_samples)
    
    # rough size of baseline rates data.table in bytes, if all included in one table
    work_size = length(covered_samples) * curr_variants[,.N] * 8
    
    # we divide into subgroups to cap baseline rates table at around 1 GB
    num_proc_groups = ceiling(work_size / 1e9)
    curr_variants[, subgroup := ceiling(num_proc_groups * 1:.N / .N)]
    
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
      
      process_variant = function(variant_id) {
        tumors_with_variant = variants_by_tumor[[variant_id]]
        curr_id = variant_id
        all_genes = gene_lookup[[variant_id]]
        
        # usually just 1 gene
        if (length(all_genes) == 1) {
          tumors_with_gene_mutated = tumors_with_variants_by_gene[[all_genes]]
        } else {
          tumors_with_gene_mutated = unique(sapply(all_genes, function(x) tumors_with_variants_by_gene[[x]]))
        }
        tumors_without_gene_mutated = setdiff(covered_samples, tumors_with_gene_mutated)
        
        rates = baseline_rates[, ..variant_id][[1]]
        names(rates) = baseline_rates[, Unique_Patient_Identifier]
        rates_tumors_with = rates[tumors_with_variant]
        rates_tumors_without = rates[tumors_without_gene_mutated]
        
        lik_args = c(list(rates_tumors_with = rates_tumors_with, rates_tumors_without = rates_tumors_without), 
                     custom_lik_args)
        if(length(sample_index) > 0) {
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
          }
        )
        
        selection_intensity = bbmle::coef(fit)
        single_stage = length(cesa@groups) == 1
        if (length(selection_intensity) == 1) {
          names(selection_intensity) = "selection_intensity"
        }
        loglikelihood = as.numeric(bbmle::logLik(fit))
        
        #dndscv_q = sapply(cesa@dndscv_out_list, function(x) x$sel_cv[x$sel_cv$gene_name == mut_record$gene, "qallsubs_cv"])
        variant_output = c(list(variant_id = variant_id), 
                           as.list(selection_intensity),
                           list(loglikelihood = loglikelihood))
        
        if(! is.null(conf)) {
          variant_output = c(variant_output, univariate_si_conf_ints(fit, fn, .001, 1e20, conf))
        }
        return(variant_output)
      }
      
      message("Calculating SIs...")
      selection_results = rbind(selection_results, rbindlist(pbapply::pblapply(curr_subgroup$variant_id, process_variant, cl = cores)))
    }
  }
  selection_results[variants, variant_type := variant_type, on = "variant_id"]
  setcolorder(selection_results, c("variant_id", "variant_type"))
  cesa@selection_results = c(cesa@selection_results, list(selection_results))
  
  return(cesa)
}






