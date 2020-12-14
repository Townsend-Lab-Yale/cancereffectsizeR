#' Gene-level epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise gene-level epistasis.
#' More specifically, selection at the gene level is assumed to act through all recurrent
#' variants (variants that appear in just one sample are ignored), with selection
#' intensity of each gene assumed to vary based on mutational status of the other gene.
#' 
#' Only samples that have coverage at all recurrent sites in both genes can be included in
#' analysis since samples lacking full coverage may or may not have mutations at the
#' uncovered sites.
#' 
#' @param cesa CESAnalysis object
#' @param genes Vector of gene names; SIs will be calculated for all gene pairs
#' @param run_name Optionally, a name to identify the current run.
#' @param conf Confidence interval size from 0 to 1 (.95 -> 95%). NULL skips calculation,
#'   which reduces runtime.
#' @param cores number of cores for parallel processing of gene pairs
#' @return CESAnalysis with epistasis analysis results added
#' @export
ces_gene_epistasis = function(cesa = NULL, genes = character(), run_name = "auto", cores = 1, conf = NULL)
{
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  cesa = update_cesa_history(cesa, match.call())
  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if(run_name %in% names(cesa@epistasis)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if (run_name == "auto") {
    # sequentially name results, handling nefarious run naming
    run_number = length(cesa@epistasis) + 1
    run_name = paste0('gene_epistasis_', run_number)
    while(run_name %in% names(cesa@epistasis)) {
      run_number = run_number + 1
      run_name = paste0('gene_epistasis_', run_number)
    }
  }
  
  setkey(cesa@samples, "Unique_Patient_Identifier") # in case dt has forgotten its key
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }

  if (length(genes) < 2) {
    stop("Supply at least two genes to analyze.")
  }
  
  maf = cesa@maf[variant_type == "snv"]
	
	
	# subset input genes to those with data
  genes = unique(genes)
	genes_to_analyze = intersect(unique(unlist(maf$genes)), genes)
	recurrent_variants = select_variants(cesa = cesa, genes = genes_to_analyze, min_freq = 2, collapse_lists = F)
	recurrent_snv_id = recurrent_variants[variant_type == "snv", variant_id]
	recurrent_aac_id = recurrent_variants[variant_type == "aac", variant_id]
	
	# use "all_genes" column because SNVs might have multiple genes, and only one appears in the "gene" column
	genes_with_recurrent_variants = unique(unlist(recurrent_variants$all_genes)) 

	genes_to_analyze = intersect(genes_with_recurrent_variants, genes_to_analyze)
	num_passing_genes = length(genes_to_analyze)
  if(num_passing_genes == 0) {
    stop("None of your requested genes have any recurrent variants (see docs for why this is required).", call. = F)
  }
	if(num_passing_genes != length(genes)) {
	  num_missing = length(genes) - num_passing_genes
	  pretty_message(paste0(num_missing, " of your requested genes had no recurrent variants, so will not be included (see docs for why this is required)."))
	}

	if (length(genes_to_analyze) == 1) {
	  stop(sprintf("Only 1 requested gene (%s) has recurrent variants, so epistasis can't be run.", genes_to_analyze), call. = F)
	}
	
	gr_table = as.data.table(get_ref_data(cesa@ref_data_dir, "gr_genes"))
	setkey(gr_table, "names")
	passing_gene_ranges = gr_table[genes_to_analyze]
	setkey(passing_gene_ranges, "seqnames", "start", "end")
	
  gene_pairs <- utils::combn(sort(genes_to_analyze), 2, simplify = F)
  bad_pairs = numeric()
  for (i in 1:length(gene_pairs)) {
    if (foverlaps(passing_gene_ranges[names == gene_pairs[[i]][1]], 
              passing_gene_ranges[names == gene_pairs[[i]][2]], nomatch = NULL)[, .N] > 0) {
      bad_pairs = c(bad_pairs, i)
    }
  }
  if (length(bad_pairs) > 0) {
    msg = paste0("The following pairs of genes were skipped because their ranges overlap in reference data:\n",
                paste(sapply(gene_pairs[bad_pairs], function(x) paste(x, collapse = "/")), collapse = ", "))
    pretty_message(msg, black = F)
    gene_pairs = gene_pairs[- bad_pairs]
    if (length(gene_pairs) == 0) {
      stop("There are no remaining gene pairs to test for epistasis.")
    }
  }
  
  selection_results = pbapply::pblapply(X = gene_pairs, FUN = pairwise_gene_epistasis, cesa=cesa, conf = conf, cl = cores)
  results = data.table::rbindlist(selection_results)
  
  # pairwise epistasis function uses v1/v2 in parameter names (as in, variants); sub in g1/g2 for gene
  colnames(results) = gsub("_v1", "_g1", colnames(results))
  colnames(results) = gsub("_v2", "_g2", colnames(results))
  
  if (results[(joint_cov_samples_just_g1_mut == 0 | joint_cov_samples_just_g2_mut == 0) & 
              joint_cov_samples_both_mut == 0, .N] > 0) {
    pretty_message(paste0("Some gene pairs had no recurrent variants in one or both genes of jointly-covering samples. ",
                   "Epistatic selection intensities are all NAs for these pairs."))
  }
  
  results = list(results)
  names(results) = run_name
  cesa@epistasis = c(cesa@epistasis, results)
  
  return(cesa)
}

#' Variant-level pairwise epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise epistasis between pairs of variants.
#' CompoundVariantSets are supported.
#' 
#' @param cesa CESAnalysis
#' @param variants To test pairs of variants, supply a list where each element is a
#'   2-length vector of CES-style variant IDs. Alternatively (and often more usefully),
#'   supply a CompoundVariantSet (see \code{define_compound_variants}) to test all pairs
#'   of compound variants in the set.
#' @param run_name Optionally, a name to identify the current run.
#' @param conf confidence interval size from 0 to 1 (.95 -> 95%); NULL skips calculation,
#'   reducing runtime.
#' @param cores number of cores for parallel processing of gene pairs
#' @return a data table with pairwise-epistatic selection intensities and variant
#'   frequencies for in tumors that have coverage at both variants in each pair
#' @export
ces_epistasis = function(cesa = NULL, variants = NULL, run_name = "auto", cores = 1, conf = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should specify a CESAnalysis object", call. = F)
  }
  cesa = update_cesa_history(cesa, match.call())
  
  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if(run_name %in% names(cesa@epistasis)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (run_name == "auto") {
    # sequentially name results, handling nefarious run naming
    run_number = length(cesa@epistasis) + 1
    run_name = paste0('variant_epistasis_', run_number)
    while(run_name %in% names(cesa@epistasis)) {
      run_number = run_number + 1
      run_name = paste0('variant_epistasis_', run_number)
    }
  }
  
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }
  
  if (is.null(cesa@mutations$snv)) {
    stop("The CESAnalysis must have mutation annotations.", call. = F)
  } 
  
  if (! is.null(variants) && length(variants) < 1) {
    stop("variants is 0-length.", call. = F)
  }
  
  if(is(variants, "CompoundVariantSet")) {
    index_pairs = utils::combn(1:length(variants), 2, simplify = F)
    results = rbindlist(pbapply::pblapply(X = index_pairs, FUN = pairwise_variant_epistasis, compound_variants = variants, cesa=cesa, conf = conf, cl = cores))
  } else if (is(variants, "list")) {
    setkey(cesa@samples, "Unique_Patient_Identifier") # in case dt has forgotten its key
    setkey(cesa@mutations$amino_acid_change, "aac_id")
    setkey(cesa@mutations$snv, "snv_id")
    
    # validate that all variant IDs appear in mutation annotations
    aac_names = cesa@mutations$amino_acid_change[, .(aac_id, variant_name = paste(gene, aachange, sep = "_"))]
    notify_name_conversion = F
    for (i in 1:length(variants)) {
      pair = variants[[i]]
      if (! is(pair, "character") || length(pair) != 2) {
        stop("All elements of variants list must be character vectors of length 2 (for two variant IDs)")
      }
      
      for (j in 1:2) {
        variant = pair[j]
        if (cesa@mutations$snv[variant, .N, nomatch = NULL] == 0 && aac_names[variant == aac_id, .N, nomatch = NULL] == 0) {
          variant = gsub(' ', '_', variant)
          record_by_name = aac_names[variant, on = "variant_name", nomatch = NULL]
          num_records = record_by_name[, .N]
          if (num_records == 0) {
            stop("Variant ", j, " could not be found. Either it's not present in the CESAnalysis, or it's a bad ID.")
          } else if (num_records > 1) {
            stop("More than one variant in the CESAnalysis has the name \"", j, "\". Use a full variant_id to specify which to test.")
          } else {
            notify_name_conversion = T
            variants[[i]][j] = record_by_name[, aac_id]
          }
        }
      }
    }
    if (notify_name_conversion) {
      pretty_message("Shorthand variant names were recognized and uniquely paired with cancereffectsizeR's aac_ids.")
    }
    results = rbindlist(pbapply::pblapply(X = variants, FUN = pairwise_variant_epistasis, cesa=cesa, conf = conf, cl = cores))
  } else {
    stop("variants should be of type list or CompoundVariantSet")
  }
  
  if(results[(joint_cov_samples_just_v1 == 0 | joint_cov_samples_just_v2 == 0) & 
             joint_cov_samples_with_both == 0, .N] > 0) {
    pretty_message(paste0("Some variant pairs had frequency of 0 in one or both variants across jointly-covering samples. ",
                          "Epistatic selection intensities are all NAs for these pairs."))
  }
  results = list(results)
  names(results) = run_name
  cesa@epistasis = c(cesa@epistasis, results)
  return(cesa)
  
}

#' Calculate SIs at variant level under pairwise epistasis model
#' @param cesa CESAnalysis
#' @param variant_pair 2-length character of variant IDs, or 2-length numeric giving
#'   indices of CompoundVariantSet for the current two compound variants
#' @keywords internal
pairwise_variant_epistasis = function(cesa, variant_pair, conf, compound_variants = NULL) {
  running_compound = FALSE
  if (is(compound_variants, "CompoundVariantSet")) {
    running_compound = TRUE
    compound_variants = compound_variants[variant_pair]
    joint_coverage = c("genome", intersect(compound_variants@compounds$shared_cov[[1]], compound_variants@compounds$shared_cov[[2]]))
    v1 = compound_variants@compounds$compound_name[[1]]
    v2 = compound_variants@compounds$compound_name[[2]]
    v1_ids = compound_variants@snvs[compound_name == v1, snv_id]
    v2_ids = compound_variants@snvs[compound_name == v2, snv_id]
    variant_ids = c(v1_ids, v2_ids)
  } else {
    v1 = variant_pair[1]
    v2 = variant_pair[2]
    
    v1_coverage = c(cesa@mutations$amino_acid_change[v1, unlist(covered_in), nomatch = NULL], 
                    cesa@mutations$snv[v1, unlist(covered_in), nomatch = NULL])
    v2_coverage = c(cesa@mutations$amino_acid_change[v2, unlist(covered_in), nomatch = NULL], 
                    cesa@mutations$snv[v2, unlist(covered_in), nomatch = NULL])
    
    # Samples have to have v1 and v2 coverage (and samples with covered_regions == "genome" always have coverage)
    joint_coverage = c("genome", intersect(v1_coverage, v2_coverage))
    variant_ids = c(v1, v2)
  }

  eligible_tumors = cesa@samples[covered_regions %in% joint_coverage, Unique_Patient_Identifier]
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
    v1_rates = apply(v1_rates, 1, function(x) 1 - prod(1 - x))
    names(v1_rates) = all_rates$Unique_Patient_Identifier
    v2_rates = all_rates[, ..v2_ids]
    v2_rates = apply(v2_rates, 1, function(x) 1 - prod(1 - x))
    names(v2_rates) = all_rates$Unique_Patient_Identifier
    
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
    setcolorder(all_rates, c("Unique_Patient_Identifier", v1, v2))
    setkey(all_rates, "Unique_Patient_Identifier")
    covered_maf = cesa@maf[eligible_tumors, on = "Unique_Patient_Identifier", nomatch = NULL]
    tumors_with_v1 = covered_maf[, v1 %in% unlist(c(variant_id, assoc_aac)), by = "Unique_Patient_Identifier"][V1 == T, Unique_Patient_Identifier]
    tumors_with_v2 = covered_maf[, v2 %in% unlist(c(variant_id, assoc_aac)), by = "Unique_Patient_Identifier"][V1 == T, Unique_Patient_Identifier]
    tumors_with_both = intersect(tumors_with_v1, tumors_with_v2)
    tumors_just_v1 = setdiff(tumors_with_v1, tumors_with_both)
    tumors_just_v2 = setdiff(tumors_with_v2, tumors_with_both)
    tumors_with_neither = setdiff(eligible_tumors, c(tumors_with_v1, tumors_with_v2))
    
    # 2-item lists: first item is baseline rates for v1; second for v2
    with_just_1 = as.list(all_rates[tumors_just_v1])[2:3]
    with_just_2 = as.list(all_rates[tumors_just_v2])[2:3]
    with_both = as.list(all_rates[tumors_with_both])[2:3]
    with_neither = as.list(all_rates[tumors_with_neither])[2:3]
  }

  # call factory function to get variant-specific likelihood function
  fn = pairwise_epistasis_lik(with_just_1, with_just_2, with_both, with_neither)
  par_init = formals(fn)[[1]]
  names(par_init) = bbmle::parnames(fn)
  
  # No point of testing epistasis if either variant doesn't appear
  if (length(tumors_with_v1) == 0 || length(tumors_with_v2) == 0) {
    early_output = list(variant1 = v1, variant2 = v2, ces_v1 = NA_real_, ces_v2 = NA_real_,
                        ces_v1_after_v2 = NA_real_, ces_v2_after_v1 = NA_real_, 
                        joint_cov_samples_just_v1 = length(tumors_just_v1),
                        joint_cov_samples_just_v2 = length(tumors_just_v2),
                        joint_cov_samples_with_both = length(tumors_with_both),
                        joint_cov_samples_with_neither = length(tumors_with_neither))
    if (! is.null(conf)) {
     si_names = names(par_init)
     ci_high_colname = paste0("ci_high_", conf * 100)
     ci_low_colname = paste0("ci_low_", conf * 100)
     low_colnames = paste(ci_low_colname, si_names, sep = "_")
     high_colnames = paste(ci_high_colname, si_names, sep = "_")
     ci_colnames = unlist(S4Vectors::zipup(low_colnames, high_colnames))
     ci_output = as.list(rep.int(NA_real_, 8))
     names(ci_output) = ci_colnames
     early_output = c(early_output, ci_output) # low/high for 4 parameters
    }
    return(early_output)
  }
  
  # find optimized selection intensities
  # the selection intensity when some group has 0 variants will be on the lower boundary; will muffle the associated warning
  withCallingHandlers(
    {
      fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = T, lower=1e-3, upper=1e9)
    },
    warning = function(w) {
      if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  params = bbmle::coef(fit)
  
  variant_ep_results = list(variant1 = v1, variant2 = v2, ces_v1 = params[1], ces_v2 = params[2],
                     ces_v1_after_v2 = params[3], ces_v2_after_v1 = params[4], 
                     joint_cov_samples_just_v1 = length(tumors_just_v1),
                     joint_cov_samples_just_v2 = length(tumors_just_v2),
                     joint_cov_samples_with_both = length(tumors_with_both),
                     joint_cov_samples_with_neither = length(tumors_with_neither))
  if(! is.null(conf)) {
    variant_ep_results = c(variant_ep_results, univariate_si_conf_ints(fit, fn, .001, 1e9, conf))
  }
  return(variant_ep_results)
}


#' Calculate SIs at gene level under pairwise epistasis model
#' 
#' The genes are assumed to not overlap in any ranges (the calling
#' function checks for this).
#' @keywords internal
pairwise_gene_epistasis = function(cesa, genes, conf) {
  gene1 = genes[1]
  gene2 = genes[2]

  # When including targeted gene sequencing data, we need to throw out any samples
  # that don't have coverage at ALL recurrent variant sites in these genes. This
  # could be a problem if some of the "exome" samples are actually whole-genome
  # if they haven't been trimmed strictly enough, especially if some of those
  # sites are low-complexity and filled with calling errors

  # Get the IDs of recurrent variants in the genes
  # Note that intergenic SNVs won't be included, even if they're close to a gene
  # Doing each gene individually to facilitate finding t
  rec_variants_g1 = select_variants(cesa, genes = gene1, min_freq = 2, collapse_lists = F)
  aac_g1 = rec_variants_g1[variant_type == "aac", variant_id]
  noncoding_g1 = rec_variants_g1[variant_type == "snv", variant_id]
  all_g1_snv = unique(c(noncoding_g1, na.omit(unlist(rec_variants_g1$constituent_snvs))))
  
  rec_variants_g2 = select_variants(cesa, genes = gene2, min_freq = 2, collapse_lists = F)
  aac_g2 = rec_variants_g2[variant_type == "aac", variant_id]
  noncoding_g2 = rec_variants_g2[variant_type == "snv", variant_id]
  all_g2_snv = unique(c(noncoding_g2, na.omit(unlist(rec_variants_g2$constituent_snvs))))
  
  rec_variants = rbind(rec_variants_g1, rec_variants_g2)
  aac = c(aac_g1, aac_g2)
  noncoding_snv = c(noncoding_g1, noncoding_g2)
  
  # Can only use samples that have coverage at all sites
  eligible_regions = c(Reduce(intersect, rec_variants$covered_in), "genome")
  eligible_tumors = cesa@samples[covered_regions %in% eligible_regions, Unique_Patient_Identifier]
  
  if(length(eligible_tumors) == 0) {
    no_tumors_with_coverage_warning = sprintf("There are no tumors that have coverage at all recurrent variant sites in genes %s and %s, so the pair must be skipped.",
                                              gene1, gene2)
    warning(no_tumors_with_coverage_warning)
    return(NULL)
  }

  # calculate rates across all sites and sum
  baseline_rates_g1 = baseline_mutation_rates(cesa, aac_ids = aac_g1, snv_ids = noncoding_g1, samples = eligible_tumors)
  g1_tmp = as.matrix(baseline_rates_g1[, -"Unique_Patient_Identifier"], rownames = baseline_rates_g1$Unique_Patient_Identifier)
  baseline_rates_g1 = apply(g1_tmp, 1, function(x) 1 - prod(1 - x))
  
  baseline_rates_g2 = baseline_mutation_rates(cesa, aac_ids = aac_g2, snv_ids = noncoding_g2, samples = eligible_tumors)
  g2_tmp = as.matrix(baseline_rates_g2[, -"Unique_Patient_Identifier"], rownames = baseline_rates_g2$Unique_Patient_Identifier)
  baseline_rates_g2 = apply(g2_tmp, 1, function(x) 1 - prod(1 - x))

  # Find recurrent mutation status of all tumors in the two genes
  maf = cesa@maf[eligible_tumors, on = "Unique_Patient_Identifier", nomatch = NULL]
  setkey(maf, "variant_id")
  tumors_with_gene1_mutated = maf[all_g1_snv, unique(Unique_Patient_Identifier), nomatch = NULL]
  tumors_with_gene2_mutated = maf[all_g2_snv, unique(Unique_Patient_Identifier), nomatch = NULL]
  
  tumors_with_both_mutated = intersect(tumors_with_gene1_mutated,tumors_with_gene2_mutated)
  tumors_with_ONLY_gene1_mutated <- setdiff(tumors_with_gene1_mutated, tumors_with_both_mutated)
  tumors_with_ONLY_gene2_mutated <- setdiff(tumors_with_gene2_mutated, tumors_with_both_mutated)
  tumors_with_neither_mutated <- setdiff(eligible_tumors, c(tumors_with_gene1_mutated, tumors_with_gene2_mutated))
  
  get_summed_mut_rates = function(tumors) {
    if(length(tumors) == 0) {
      return(NULL)
    }
    rates1 = unname(baseline_rates_g1[tumors])
    rates2 = unname(baseline_rates_g2[tumors])
    return(list(rates1, rates2))
  }
  
  with_just_1 = get_summed_mut_rates(tumors_with_ONLY_gene1_mutated)
  with_just_2 = get_summed_mut_rates(tumors_with_ONLY_gene2_mutated)
  with_both = get_summed_mut_rates(tumors_with_both_mutated)
  with_neither = get_summed_mut_rates(tumors_with_neither_mutated)
  
  fn = pairwise_epistasis_lik(with_just_1, with_just_2, with_both, with_neither)
  par_init = formals(fn)[[1]]
  names(par_init) = bbmle::parnames(fn)
  
  if(length(tumors_with_gene1_mutated) == 0 || length(tumors_with_gene2_mutated) == 0) {
    early_output = list(gene_1 = gene1, gene_2 = gene2, ces_g1 = NA_real_, ces_g2 = NA_real_,
                           ces_g1_after_g2 = NA_real_, ces_g2_after_g1 = NA_real_, 
                           joint_cov_samples_just_g1_mut = length(tumors_with_ONLY_gene1_mutated),
                           joint_cov_samples_just_g2_mut = length(tumors_with_ONLY_gene2_mutated),
                           joint_cov_samples_both_mut = length(tumors_with_both_mutated),
                           joint_cov_samples_no_mut = length(tumors_with_neither_mutated))
    if(! is.null(conf)) {
      si_names = names(par_init)
      ci_high_colname = paste0("ci_high_", conf * 100)
      ci_low_colname = paste0("ci_low_", conf * 100)
      low_colnames = paste(ci_low_colname, si_names, sep = "_")
      high_colnames = paste(ci_high_colname, si_names, sep = "_")
      ci_colnames = unlist(S4Vectors::zipup(low_colnames, high_colnames))
      ci_output = as.list(rep.int(NA_real_, 8))
      names(ci_output) = ci_colnames
      early_output = c(early_output, ci_output) # low/high for 4 parameters
    }
    return(early_output)
  }
  
  
  
  # find optimized selection intensities
  # the selection intensity when some group has 0 variants will be on the lower boundary; will muffle the associated warning
  withCallingHandlers(
    {
      fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = T, lower=1e-3, upper=1e9)
    },
    warning = function(w) {
      if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  
  params = bbmle::coef(fit)
  gene_ep_results = list(gene_1 = gene1, gene_2 = gene2, ces_g1 = params[1], ces_g2 = params[2],
              ces_g1_after_g2 = params[3], ces_g2_after_g1 = params[4], 
              joint_cov_samples_just_g1_mut = length(tumors_with_ONLY_gene1_mutated),
              joint_cov_samples_just_g2_mut = length(tumors_with_ONLY_gene2_mutated),
              joint_cov_samples_both_mut = length(tumors_with_both_mutated),
              joint_cov_samples_no_mut = length(tumors_with_neither_mutated))
  
  
  if(! is.null(conf)) {
     gene_ep_results = c(gene_ep_results, univariate_si_conf_ints(fit, fn, .001, 1e9, conf))
  }
  
  return(gene_ep_results)
}
