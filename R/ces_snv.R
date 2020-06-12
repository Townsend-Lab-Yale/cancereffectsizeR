#' Calculate selection intensity for single-nucleotide variants and amino acid changes
#' @param cesa CESAnalysis object
#' @param gene which genes to calculate effect sizes within; defaults to all genes with recurrent mutations in data set
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
  setkey(mutations$amino_acid_change, "aa_mut_id")
  setkey(mutations$snv, "snv_id")
  
  if(length(genes_in_dataset) == 0) {
    stop("There are no annotated mutations in the data set!", call. = F)
  }
  
  # filter variants based on recurrency
  if(include_nonrecurrent_variants == T) {
    aac_ids = unique(mutations$amino_acid_change$aa_mut_id)
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

  snvs_in_analysis = union(noncoding_snv_ids, mutations$amino_acid_change[aac_ids, unlist(all_snv_ids)])
  
  sample_gene_rates = as.data.table(expand.grid(gene = genes_in_analysis, Unique_Patient_Identifier = cesa@samples$Unique_Patient_Identifier, stringsAsFactors = F),
                                    key = "Unique_Patient_Identifier")
  
  sample_gene_rates = sample_gene_rates[cesa@samples[, .(Unique_Patient_Identifier, progression_name)]]
  
  melted_mutrates = melt.data.table(cesa@mutrates, id.vars = c("gene"))
  setnames(melted_mutrates, c("variable", "value"), c("progression_name", "raw_rate"))
  
  sample_gene_rates = melted_mutrates[sample_gene_rates, , on = c("gene", "progression_name")]
  dt = as.data.table(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix, keep.rownames = "Unique_Patient_Identifier")
  setkey(dt, "Unique_Patient_Identifier")
  
  gene_trinuc_comp = get_genome_data(cesa, "gene_trinuc_comp")
  
  sample_gene_rates[, trinuc_comp := lapply(gene, function(x) gene_trinuc_comp[[x]])]
  N = sample_gene_rates[, .N]
  for (i in 1:N) {
    # Performance note: If you don't explicitly convert the dt row with as.numeric, runs much slower
    sample_gene_rates[i, aggregate_rate := raw_rate / sum(gene_trinuc_comp[[gene]] * as.numeric(dt[Unique_Patient_Identifier, 2:97]))]
  }
  
  baseline_rates = dt[, .(Unique_Patient_Identifier)]
  for (aac in aac_ids) {
    curr_gene = mutations$amino_acid_change[aac, gene]
    
    # Get SNV IDs associated with the aac
    snv_ids = mutations$amino_acid_change[aac, unlist(all_snv_ids)]
    trinuc_mut = mutations$snv[snv_ids, trinuc_mut]
    tmp = dt[, sum(as.numeric(.SD)), .SDcols = trinuc_mut, by = "Unique_Patient_Identifier"]
    tmp = sample_gene_rates[gene == curr_gene][tmp, , on = "Unique_Patient_Identifier"]
    baseline_rates[, (aac) := tmp$aggregate_rate * tmp$V1]
  }
  
  for (snv in noncoding_snv_ids) {
    # occasionally more than one gene annotated; in this case take the average
    curr_genes = mutations$snv[snv, unlist(genes)]
    trinuc_mut = mutations$snv[snv, trinuc_mut]
    tmp = dt[, sum(as.numeric(.SD)), .SDcols = trinuc_mut, by = "Unique_Patient_Identifier"]
    tmp =  sample_gene_rates[gene %in% curr_genes][tmp, , on = "Unique_Patient_Identifier"]
    tmp[, per_gene_rate := tmp$aggregate_rate * tmp$V1]
    averaged_over_genes = tmp[, mean(per_gene_rate), by = "Unique_Patient_Identifier"]
    baseline_rates[, (snv) := averaged_over_genes$V1]
  }

  return(baseline_rates) # for now
  selection_results <- rbindlist(pbapply::pblapply(genes_in_analysis, get_gene_results, cesa = cesa, conf = conf,
                                             gene_trinuc_comp = gene_trinuc_comp, cl = cores))
  cesa@selection_results = selection_results
  cesa@status[["SNV selection"]] = "view effect sizes with snv_results()"
  return(cesa)
}


#' Single-stage SNV effect size analysis (gets called by ces_snv)
#' @keywords internal
get_gene_results <- function(gene, cesa, conf, gene_trinuc_comp) {
  if(! is.null(conf)) {
    ci_high_colname = paste0("ci_high_", conf * 100)
    ci_low_colname = paste0("ci_low_", conf * 100)
  }
  snv.maf = cesa@maf[Variant_Type == "SNV"]
  current_gene_maf = snv.maf[Gene_name == gene]
  these_mutation_rates <-
    mutation_rate_calc(
      this_MAF = current_gene_maf,
      gene = gene,
      gene_mut_rate = cesa@mutrates,
      trinuc_proportion_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix,
      gene_trinuc_comp = gene_trinuc_comp,
      samples = cesa@samples)

  variants = colnames(these_mutation_rates)
  process_variant = function(variant) {
    # use the first matching record as the locus 
    # (will assume that for amino acid variants, coverage at one site in codon implies coverage for whole codon)
    variant_maf = current_gene_maf[nt_mut_id == variant] # no need to subset further because already dealing with a gene-specific MAF
    # covered_in is a 1-item list with a character vector of coverage_grs that cover the variant site
    site_coverage = unlist(variant_maf[1, covered_in])
    eligible_tumors = cesa@samples[covered_regions %in% site_coverage, Unique_Patient_Identifier]
    eligible_tumors[eligible_tumors %in% rownames(these_mutation_rates)] # To-do: This check isn't necessary if these_mutation_rates covers all tumors

    
    # given the tumors with coverage, their mutation rates at the variant sites, and their mutation status,
    # find most likely selection intensities (by stage if applicable)
    tumors_with_pos_mutated <- variant_maf[nt_mut_id==variant, Unique_Patient_Identifier]
    tumors_without_gene_mutated <- eligible_tumors[! eligible_tumors %in% current_gene_maf$Unique_Patient_Identifier]
    tumor_stage_indices = cesa@samples[eligible_tumors, progression_index]
    names(tumor_stage_indices) = eligible_tumors
    fn = ml_objective(tumor_stages = tumor_stage_indices, tumors_without_gene_mutated = tumors_without_gene_mutated,
                      tumors_with_pos_mutated = tumors_with_pos_mutated, variant=variant, specific_mut_rates=these_mutation_rates)
    
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
    unsure_gene_name = rep(variant_maf$unsure_gene_name[1], length(selection_intensity))
    progression_name = cesa@progressions
    if (length(cesa@progressions) == 1) {
      progression_name = "Not applicable"
    }
    
    # Get number of tumors of each named stage with the variant (in proper progression order)
    stages = cesa@samples[tumors_with_pos_mutated, progression_name]
    tumors_with_variant = as.numeric(table(factor(stages, levels = cesa@progressions)))
    
    # Also get number of eligible tumors per stage
    stages = cesa@samples[eligible_tumors, progression_name]
    tumors_with_coverage = as.numeric(table(factor(stages, levels = cesa@progressions)))
    
    # if a tumor has no coverage in a given stage, set selection intensity to NA
    uncovered_stages = which(tumors_with_coverage == 0)
    selection_intensity[uncovered_stages] = NA
    
    dndscv_q = sapply(cesa@dndscv_out_list, function(x) x$sel_cv[x$sel_cv$gene_name == gene, "qallsubs_cv"])
    
    # legacy_id = variant_maf$unique_variant_ID[1]
    # if (variant_maf$is_coding[1] == TRUE) {
    #   legacy_id = variant_maf[1, paste(Gene_name, unique_variant_ID_AA)]
    # }

    variant_id = ifelse(is.na(variant_maf$aa_mut_id[1]), variant, variant_maf$aa_mut_id[1])
    
    
    
    variant_output = data.table(variant = variant_id, selection_intensity, unsure_gene_name, loglikelihood, gene, 
                         progression = progression_name, tumors_with_variant, tumors_with_coverage, dndscv_q,
                         legacy_id)
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
  return(data.table::rbindlist(lapply(variants, process_variant)))
}





