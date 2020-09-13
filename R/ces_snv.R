#' Calculate selection intensity for single-nucleotide variants and amino acid changes
#' @param cesa CESAnalysis object
#' @param genes limit analysis to just these genes (otherwise, all genes, unless variants is specified)
#' @param cores number of cores to use
#' @param conf selection intensity confidence interval width (NULL skips calculation, speeds runtime)
#' @param variant_ids When specified, CES-style variant IDs of variants to include, or
#'   short names like "KRAS G12C" (or KRAS_G12C). If "genes" is also specified, the union
#'   is taken. If you want to use this parameter to calculate SI confidence intervals at
#'   sites with no MAF variants, set min_freq = 0.
#' @param min_freq default 2; setting to 0 or 1 can be useful for
#'   comparing sample groups or analyzing SIs across genomic regions
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export

ces_snv <- function(cesa = NULL,
                    genes = NULL,
                    min_freq = 2,
                    variant_ids = NULL,
                    cores = 1,
                    conf = .95) 
{
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis", call. = F)
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
  

  
  selected = select_variants(cesa, genes = genes, variant_ids = variant_ids, min_freq = min_freq)
  if (is.null(selected)) {
    stop("No variants passed your filters!")
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
  
  # identify mutations by variant
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
    tumor_stage_indices = cesa@samples[eligible_tumors, group_index]
    names(tumor_stage_indices) = eligible_tumors
    rates = baseline_rates[, ..mut_id][[1]]
    names(rates) = baseline_rates[, Unique_Patient_Identifier]
    
    fn = ml_objective(tumor_stages = tumor_stage_indices, tumors_without_gene_mutated = tumors_without_gene_mutated,
                      tumors_with_variant = tumors_with_variant, baseline_rates = rates)
    
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
  cesa@selection_results = selection_results
  
  return(cesa)
}








