#' Calculate selection intensity for single-nucleotide variants and amino acid changes
#' 
#' @param cesa CESAnalysis object
#' @param cores number of cores to use
#' @param conf selection intensity confidence interval width (NULL skips calculation,
#'   speeds runtime)
#' @param variants Variants to include, obtained from calling select_variants() on the
#'   CESAnalysis. Defaults to selecting all recurrent amino-acid-changing SNVs plus all
#'   recurrent noncoding SNVs, where recurrent means appearing at least twice across the
#'   full MAF data set. Calling select_variants with `min_freq` set to 0 or 1 can be
#'   useful for comparing sample groups or analyzing SIs across genomic regions.
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
  setkey(cesa@samples, "Unique_Patient_Identifier")
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  mutations = cesa@mutations
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }
  
  selected = variants
  if (is.null(selected) || ! is.data.table(selected)) {
    stop("variants expected to be a data.table (from select_variants()).")
  }
  if (selected[, .N] == 0) {
    stop("There are no variants in the input!")
  }
  aac_ids = selected[variant_type == "aac", variant_id]
  
  # By noncoding, we just mean that SIs are calculated at the SNV site rather than at the AAC level,
  # regardless of whether there's a CDS annotation. We're not going to check to see if the user's
  # AACs and SNVs overlap (that's their problem if they decided to supply overlapping variants)
  noncoding_snv_ids = selected[variant_type == "snv", variant_id]
  
  if(length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", call. = F)
  }

  # get SNVs and AACs of interest
  noncoding_table = mutations$snv[noncoding_snv_ids]
  setkey(noncoding_table, "snv_id")
  coding_table = mutations$amino_acid_change[aac_ids]
  setkey(coding_table, "aac_id")
  
  # identify mutations by nearest gene(s)
  tmp = cesa@maf[variant_type == "snv", .(gene = unlist(genes)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  
  # identify mutations by sample
  tmp = cesa@maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "aac_id"]
  tmp = tmp[coding_table[, aac_id], on = "aac_id"]
  aacs_by_tumor = tmp$samples
  names(aacs_by_tumor) = tmp$aac_id
  
  wgs_samples = cesa@samples[covered_regions == "genome", Unique_Patient_Identifier]
  

  
  # function takes in AAC or SNV ID, returns SI table output
  process_variant = function(mut_id, snv_or_aac) {
    if(snv_or_aac == "aac") {
      mut_record = coding_table[mut_id]
      tumors_with_variant = aacs_by_tumor[[mut_id]]
      tumors_with_gene_mutated = tumors_with_variants_by_gene[[mut_record$gene]]
    } else {
      mut_record = noncoding_table[mut_id]
      tumors_with_variant = cesa@maf[variant_id == mut_id, Unique_Patient_Identifier]
      tumors_with_gene_mutated = unlist(tumors_with_variants_by_gene[unlist(mut_record$genes)], use.names = F) # for rare case of multiple gene hits, take all
    }
    
    eligible_tumors = cesa@samples[covered_regions %in% unlist(mut_record$covered_in), Unique_Patient_Identifier]
    eligible_tumors = union(eligible_tumors, wgs_samples)
    tumors_without_gene_mutated = setdiff(eligible_tumors, tumors_with_gene_mutated)

    rates = baseline_rates[, ..mut_id][[1]]
    names(rates) = baseline_rates[, Unique_Patient_Identifier]
    
    lik_args = c(list(rates_tumors_with = rates[tumors_with_variant], rates_tumors_without = rates[tumors_without_gene_mutated]), 
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
    variant_output = c(list(variant_id = mut_id, variant_type = snv_or_aac), 
                       as.list(selection_intensity),
                       list(loglikelihood = loglikelihood))
    
    if(! is.null(conf)) {
      variant_output = c(variant_output, univariate_si_conf_ints(fit, fn, .001, 1e20, conf))
    }
    return(variant_output)
  }
  

  # Will process variants by coverage group (i.e., groups of variants that have the same tumors covering them)
  selection_results = NULL
  coverage_groups = unique(c(coding_table$covered_in, noncoding_table$covered_in))
  num_coverage_groups = length(coverage_groups)
  
  for (i in 1:num_coverage_groups) {
    message(sprintf("Preparing to calculate selection intensities (batch %i of %i)...", i, num_coverage_groups))
    aac_ids = coding_table[, identical(unlist(covered_in),coverage_groups[[i]]), by = "aac_id"][V1 == T, aac_id]
    noncoding_snv_ids = noncoding_table[, identical(unlist(covered_in),coverage_groups[[i]]), by = "snv_id"][V1 == T, snv_id]
    covered_samples = cesa@samples[covered_regions %in% coverage_groups[[i]], Unique_Patient_Identifier]
    
    muts_in_group = data.table(mut_id = c(aac_ids, noncoding_snv_ids), snv_or_aac = c(rep.int("aac", length(aac_ids)), 
                                                            rep.int("snv", length(noncoding_snv_ids))))
    
    # rough size of baseline rates data.table in bytes, if all included in one table
    work_size = length(covered_samples) * nrow(muts_in_group) * 8
    
    # we divide into subgroups to cap basline rates table size at approx. 1 GB
    num_proc_groups = ceiling(work_size / 1e9)
    muts_in_group[, subgroup := ceiling(num_proc_groups * 1:.N / .N)]
    
    for (j in 1:num_proc_groups) {
      if (num_proc_groups > 1) {
        message(sprintf("Working on sub-batch %i of %i...", j, num_proc_groups))
      }
      muts_in_subgroup = muts_in_group[subgroup == j]
      aac_ids = muts_in_subgroup[snv_or_aac == "aac", mut_id]
      snv_ids = muts_in_subgroup[snv_or_aac == "snv", mut_id]

      baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, snv_ids = snv_ids, samples = covered_samples, cores = cores)
      message("Calculating SIs for coding mutations...")
      selection_results = rbind(selection_results, rbindlist(pbapply::pblapply(aac_ids, process_variant, snv_or_aac = "aac", cl = cores)))
      message("Calculating SIs for noncoding SNVs...")
      selection_results = rbind(selection_results, rbindlist(pbapply::pblapply(snv_ids, process_variant, snv_or_aac = "snv", cl = cores)))
    }
  }
  cesa@selection_results = c(cesa@selection_results, list(selection_results))
  
  return(cesa)
}








