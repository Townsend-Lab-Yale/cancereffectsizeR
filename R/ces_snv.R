#' Calculate selection intensity for single-nucleotide variants and amino acid changes
#' @param cesa CESAnalysis object
#' @param gene which genes to calculate effect sizes within; defaults to all genes in data set
#' @param cores number of cores to use
#' @param conf selection intensity confidence interval width (NULL skips calculation, speeds runtime)
#' @param include_nonrecurrent_variants default false; will increase runtime and SIs at non-recurrent sites aren't very informative
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export

ces_snv <- function(cesa = NULL,
                            genes = "all",
                            cores = 1,
                            include_nonrecurrent_variants = F,
                            conf = .95) 
{
  setkey(cesa@samples, "Unique_Patient_Identifier") # in case dt has forgotten its key
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }
  
  if (! "character" %in% class(genes)) {
    stop("Expected argument \"genes\" to take a character vector.", call. = F)
  }

  mutations = cesa@mutations
  if (length(mutations) == 0) {
    stop("There are no annotated mutations in the analysis!", call. = F)
  }
  
  
  # take all SNVs in data set, then subtract any that are a part of aac mutations
  setkey(mutations$amino_acid_change, "aac_id")
  setkey(mutations$snv, "snv_id")
  
  
  # filter variants based on recurrency
  if(include_nonrecurrent_variants == T) {
    aac_ids = unique(mutations$amino_acid_change$aac_id)
    noncoding_snv_ids = setdiff(cesa@maf[! is.na(snv_id), snv_id], mutations$amino_acid_change[aac_ids, unlist(all_snv_ids)])
    
  } else {
    aac_ids = cesa@maf[! is.na(assoc_aa_mut), .(aac_id = unlist(assoc_aa_mut))][, .N, by = aac_id][N > 1, aac_id]
    recurrent_snvs = cesa@maf[! is.na(snv_id), .(snv_id)][, .N, by = "snv_id"][N > 1, snv_id]
    noncoding_snv_ids = setdiff(recurrent_snvs, mutations$amino_acid_change[aac_ids, unlist(all_snv_ids)])
  }
  
  # filter based on gene
  genes_in_dataset = union(mutations$amino_acid_change[aac_ids, gene], mutations$snv[noncoding_snv_ids, unlist(genes)])
  if(genes[1] =="all") {
    genes_in_analysis <- genes_in_dataset
  } else{
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
    
    # determine which noncoding SNVs have a gene annotation containg the analysis genes, and include just those
    noncoding_snv_ids = mutations$snv[noncoding_snv_ids, .(include = any(unlist(genes) %in% genes_in_analysis)), by = snv_id][include == T, snv_id]
  }

  message("Collecting variants and determining baseline mutation rates...")
  
  baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, snv_ids = noncoding_snv_ids)

  
  # get SNVs and AACs of interest
  noncoding_table = mutations$snv[snv_id %in% noncoding_snv_ids]
  coding_table = mutations$amino_acid_change[aac_id %in% aac_ids]
  
  # identify mutations by gene
  tmp = cesa@maf[Variant_Type == "SNV", .(gene = unlist(genes)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  
  # identify mutations by variant
  tmp = cesa@maf[! is.na(assoc_aa_mut), .(aac_id = unlist(assoc_aa_mut)), by = "Unique_Patient_Identifier"][, .(samples = list(Unique_Patient_Identifier)), by = "aac_id"]
  tmp = tmp[aac_id %in% coding_table$aac_id]
  aacs_by_tumor = tmp$samples
  names(aacs_by_tumor) = tmp$aac_id

  if(! is.null(conf)) {
    ci_high_colname = paste0("ci_high_", conf * 100)
    ci_low_colname = paste0("ci_low_", conf * 100)
  }
  
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
      fn <<- fit@minuslogl
      
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
  
  
  
  
  # start with aac SIs
  message("Estimating selection intensities for amino acid changes...")
  aac_results = rbindlist(pbapply::pblapply(aac_ids, process_variant, snv_or_aac = "aac", cl = cores))
  message("Estimating selection intensities for noncoding SNVs...")
  snv_results = rbindlist(pbapply::pblapply(noncoding_snv_ids, process_variant, snv_or_aac = "snv", cl = cores))
  cesa@selection_results = rbind(aac_results, snv_results)
  
  cesa@status[["SNV selection"]] = "view effect sizes with snv_results()"
  return(cesa)
}








