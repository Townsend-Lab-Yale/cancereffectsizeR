#' Calculate selection intensity for single-nucleotide variants and amino acid changes
#' @param cesa CESAnalysis object
#' @param genes which genes to calculate effect sizes within; defaults to all genes in data set
#' @param cores number of cores to use
#' @param conf selection intensity confidence interval width (NULL skips calculation, speeds runtime)
#' @param variants When specified, list of specific mutations to include, excluding all others; 
#'                 e.g., list(aac_id = c(...), snv_id = c(...)). If "genes" is also specified, the intersection is taken. If
#'                 you want to use this parameter to calculate SI at sites with no MAF variants, set include_nonrecurrent_variants = T.
#' @param include_nonrecurrent_variants default false; will increase runtime and SIs at non-recurrent sites aren't very informative
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export

ces_snv <- function(cesa = NULL,
                    genes = "all",
                    cores = 1,
                    include_nonrecurrent_variants = F,
                    variants = NULL,
                    conf = .95) 
{
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
  
  if (! "character" %in% class(genes)) {
    stop("Expected argument \"genes\" to take a character vector.", call. = F)
  }
  
  if (length(mutations) == 0) {
    stop("There are no annotated mutations in the analysis!", call. = F)
  }
  
  # If no list of variants is specified, take all AACs in data set and all SNVs that are not covered in AACs
  if (is.null(variants)) {
    aac_ids = unique(na.omit(unlist(cesa@maf$assoc_aa_mut)))
    
    # tak all SNVs in MAF, then subtract those that have AAC annotation
    noncoding_snv_ids = setdiff(cesa@maf[! is.na(snv_id), snv_id], mutations$amino_acid_change[aac_ids, unlist(constituent_snvs)])
  } else {
    if(! is(variants, "list")) {
      stop("Expected variants to be a list", call. = F)
    }
    if (! all(names(variants) %in% c("snv_id", "aac_id"))) {
      stop("variants must be a named list containing snv_id and/or aac_id")
    }
    if(length(variants) != length(unique(names(variants)))) {
      stop("variants must be a named list containing snv_id and/or aac_id")
    }
    
    if(! unique(sapply(variants, class)) == "character") {
      stop("Elements of variants list must be character vectors of variant IDs")
    }
    
    # By noncoding, we just mean that SI will be calculate at the SNV site rather than at the AAC level,
    # regardless of whether there's a CDS annotation. We're not going to check to see if the user's
    # AACs and SNVs overlap (that's their problem to deal with!)
    noncoding_snv_ids = character()
    if (! is.null(variants$snv_id)) {
      noncoding_snv_ids = unique(variants$snv_id)
      num_present = mutations$snv[noncoding_snv_ids, .N, nomatch = NULL]
      num_missing = length(noncoding_snv_ids) - num_present
      if (num_missing > 0) {
        stop(paste0("There are no annotations (in $mutations$snv) for ", num_missing, " of your requested SNVs."), call. = F)
      }
    }
    
    aac_ids = character()
    if (! is.null(variants$aac_id)) {
      aac_ids = unique(variants$aac_id)
      num_present = mutations$amino_acid_change[aac_ids, .N, nomatch = NULL]
      num_missing = length(aac_ids) - num_present
      if (num_missing > 0) {
        stop(paste0("There are no annotations (in $mutations$amino_acid_change) for ", num_missing, " of your requested coding mutations."), call. = F)
      }
    }
  }
  
  
  # filter variants based on recurrency
  if(include_nonrecurrent_variants == F) {
    recurrent_aac = cesa@maf[! is.na(assoc_aa_mut), .(aac_id = unlist(assoc_aa_mut))][, .N, by = aac_id][N > 1, aac_id]
    aac_ids = aac_ids[aac_ids %in% recurrent_aac]
    
    recurrent_snv = cesa@maf[! is.na(snv_id), .(snv_id)][, .N, by = "snv_id"][N > 1, snv_id]
    noncoding_snv_ids = noncoding_snv_ids[noncoding_snv_ids %in% recurrent_snv]
  }
  
  # filter based on gene
  genes_in_dataset = unlist(cesa@maf$genes)
  if(genes[1] != "all") {
    genes = unique(genes)
    genes_in_analysis <- genes[genes %in% genes_in_dataset]
    missing_genes = genes[! genes %in% genes_in_dataset]
    gene_names = get_genome_data(cesa, "gene_names")
    invalid_genes = missing_genes[! missing_genes %in% gene_names]
    num_invalid = length(invalid_genes)
    if(num_invalid > 0) {
      additional_msg = ""
      if(num_invalid > 50) {
        invalid_genes = invalid_genes[1:40]
        additional_msg = paste0(" (and ", num_invalid - 40, " more)")
      }
      list_of_invalid = paste(invalid_genes, collapse = ", ")
      stop(paste0("Note: The following requested genes were not found in reference data for your genome build:\n\t",
                  list_of_invalid, additional_msg, "\n"))
    }
    if (length(genes_in_analysis) == 0) {
      stop("None of the requested genes have eligible mutations in the SNV data set.")
    }

    num_missing = length(missing_genes)
    if(num_missing > 0) {
      additional_msg = ""
      if(num_missing > 50) {
        missing_genes = missing_genes[1:40]
        additional_msg = paste0(" (and ", num_missing - 40, " more)")
      }
      list_of_missing = paste(missing_genes, collapse = ", ")

      message(paste0("The following requested genes have no eligible mutations in the SNV data set:\n\t",
        list_of_missing, additional_msg))
    }
    aac_ids = aac_ids[mutations$amino_acid_change[aac_ids, gene %in% genes_in_analysis]]
    
    # determine which noncoding SNVs have a gene annotation containing the analysis genes, and include just those
    noncoding_snv_ids = mutations$snv[noncoding_snv_ids, .(include = any(unlist(genes) %in% genes_in_analysis)), by = snv_id][include == T, snv_id]
  }
  
  if(length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", call. = F)
  }


  # get SNVs and AACs of interest
  noncoding_table = mutations$snv[noncoding_snv_ids]
  setkey(noncoding_table, "snv_id")
  coding_table = mutations$amino_acid_change[aac_ids]
  setkey(coding_table, "aac_id")
  
  # identify mutations by gene
  tmp = cesa@maf[Variant_Type == "SNV", .(gene = unlist(genes)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  
  # identify mutations by variant
  tmp = cesa@maf[! is.na(assoc_aa_mut), .(aac_id = unlist(assoc_aa_mut)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "aac_id"]
  tmp = tmp[coding_table[, aac_id], on = "aac_id"]
  aacs_by_tumor = tmp$samples
  names(aacs_by_tumor) = tmp$aac_id

  # include upper/lower CI in column name
  if(! is.null(conf)) {
    ci_high_colname = paste0("ci_high_", conf * 100)
    ci_low_colname = paste0("ci_low_", conf * 100)
  }
  
  # function takes in AAC or SNV ID, returns SI table output
  process_variant = function(mut_id, snv_or_aac) {
    if(snv_or_aac == "aac") {
      mut_record = coding_table[mut_id]
      tumors_with_variant = aacs_by_tumor[[mut_id]]
      tumors_with_gene_mutated = tumors_with_variants_by_gene[[mut_record$gene]]
    } else {
      mut_record = noncoding_table[mut_id]
      tumors_with_variant = cesa@maf[snv_id == mut_id, Unique_Patient_Identifier]
      tumors_with_gene_mutated = unlist(tumors_with_variants_by_gene[unlist(mut_record$genes)], use.names = F) # for rare case of multiple gene hits, take all
    }
    eligible_tumors = cesa@samples[covered_regions %in% unlist(mut_record$covered_in), Unique_Patient_Identifier]
    
    tumors_without_gene_mutated = setdiff(eligible_tumors, tumors_with_gene_mutated)
    tumor_stage_indices = cesa@samples[eligible_tumors, progression_index]
    names(tumor_stage_indices) = eligible_tumors
    rates = baseline_rates[, ..mut_id][[1]]
    names(rates) = baseline_rates[, Unique_Patient_Identifier]
    fn = ml_objective(tumor_stages = tumor_stage_indices, tumors_without_gene_mutated = tumors_without_gene_mutated,
                      tumors_with_variant = tumors_with_variant, baseline_rates = rates)
    
    # initialize all gamma (SI) values at 1000; bbmle requires a parnames attribute be set to name each gamma (here, g1, g2, etc.)
    par_init = rep(1000, length(cesa@progressions))
    names(par_init) <- bbmle::parnames(fn) <- paste0("g", 1:length(cesa@progressions))
    
    # find optimized selection intensities
    # the selection intensity for any stage that has 0 variants will be on the lower boundary; will muffle the associated warning
    withCallingHandlers(
      {
        fit = bbmle::mle2(fn, start = par_init, method="L-BFGS-B", vecpar = T, lower=1e-3, upper=1e9, control=list(fnscale=1e-12), hessian.opts = list(method = "complex"))
      },
      warning = function(w) {
        if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
          invokeRestart("muffleWarning")
        }
      }
    )
    
    selection_intensity =  bbmle::coef(fit)
    loglikelihood = as.numeric(bbmle::logLik(fit))
    loglikelihood = rep(loglikelihood, length(selection_intensity))
    progression_name = cesa@progressions
    if (length(cesa@progressions) == 1) {
      progression_name = "Not applicable"
    }
    
    # Get number of tumors of each named stage with the variant (in proper progression order)
    stages = cesa@samples[tumors_with_variant, progression_name]
    tumors_with_variant = as.numeric(table(factor(stages, levels = cesa@progressions)))
    
    # Also get number of eligible tumors per stage
    stages = cesa@samples[eligible_tumors, progression_name]
    tumors_with_coverage = as.numeric(table(factor(stages, levels = cesa@progressions)))
    
    # if a tumor has no coverage in a given stage, set selection intensity to NA
    uncovered_stages = which(tumors_with_coverage == 0)
    selection_intensity[uncovered_stages] = NA
    
    #dndscv_q = sapply(cesa@dndscv_out_list, function(x) x$sel_cv[x$sel_cv$gene_name == mut_record$gene, "qallsubs_cv"])
    
    variant_output = data.table(variant = mut_id, variant_type = snv_or_aac, selection_intensity, loglikelihood, 
                                progression = progression_name, tumors_with_variant, tumors_with_coverage)
    if(! is.null(conf)) {
      # can't get confidence intervals for progression states that have no tumors with the variant
      offset = qchisq(conf, 1)/2
      max_ll = -1 * loglikelihood[1]
      
      # for each SI, use uniroot to get a single-parameter confidence interval
      for (i in 1:length(cesa@progressions)) {
        if(is.na(selection_intensity[i])) {
          lower = NA_real_
          upper = NA_real_
        } else {
          # univariate likelihood function freezes all but one SI at MLE
          # offset makes output at MLE negative; function should be positive at lower/upper boundaries (.001, 1e20),
          # and uniroot will find the zeroes, which should represent the lower/uppper CIs
          ulik = function(x) { 
            pars = selection_intensity
            pars[i] = x
            return(fn(pars) - max_ll - offset)
          }
          # if ulik(.001) is negative, no root on (.001, SI), so assigning .001 as lower bound
          if(ulik(.001) < 0) {
            lower = .001
          } else {
            # Enforcing an SI floor of .001, as in optimization
            lower = max(uniroot(ulik, lower = .001, upper = selection_intensity[i])$root, .001)
          }
          if(ulik(1e20) < 0){
            # this really shouldn't happen
            upper = NA_real_
          } else {
            upper = uniroot(ulik, lower = selection_intensity[i], upper = 1e20)$root
          }
        }
        variant_output[i, c(ci_low_colname) := lower]
        variant_output[i, c(ci_high_colname) := upper]
      }
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








