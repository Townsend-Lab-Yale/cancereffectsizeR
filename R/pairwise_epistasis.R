#' Gene-level epistasis
#' 
#' Calculate selection intensity under an assumption of pairwise gene-level epistasis.
#' Selection at the gene level is assumed to act through all 
#' specified variants (see options), with the selection for each gene's variants allowed
#' to vary based on the mutational status of the other gene's variants.
#' 
#' Only samples that have coverage at all included sites in both genes can be included in the
#' inference since samples lacking full coverage may or may not have mutations at the
#' uncovered sites.
#' 
#' @param cesa CESAnalysis object
#' @param genes Vector of gene names; SIs will be calculated for all gene pairs. Alternatively,
#' a list of gene pairs (2-length character vectors) to run just the given pairings.
#' @param variants Which variants to include in inference for each gene. Either "recurrent" for all
#'   variants present in two or more samples (across all MAF data), "nonsilent" for nonsynonymous
#'   coding variants and variants in essential splice sites, or a data.table containing all variants
#'   to include (as returned by \code{select_variants()} or by
#'   subsetting\code{[CESAnalysis\]$variants}). For noncoding variants with multiple gene
#'   annotations, the one listed in the "gene" column is used. In the recurrent method, nearby
#'   noncoding variants may be included.
#' @param samples Which samples to include in inference. Can be a vector of
#'   Unique_Patient_Identifiers, or a data.table containing rows from the CESAnalysis sample table.
#'   Samples that do not have coverage at all variant sites in a given inference will be set aside.
#' @param run_name Optionally, a name to identify the current run.
#' @param conf Confidence interval size from 0 to 1 (.95 -> 95\%). NULL skips calculation,
#'   which may be helpful to reduce runtime when analyzing many gene pairs.
#' @param cores Number of cores for parallel processing of gene pairs.
#' @param return_fit TRUE/FALSE (default FALSE): Embed epistatic model fits for each gene pair in a "fit" attribute
#' of the epistasis results table. Use \code{attr(my_results, 'fit')} to access the list of fitted models.
#' @return CESAnalysis with a table of epistatic inferences appended to list \code{[CESAnalysis]$epistasis}. Some column definitions:
#'   \itemize{
#'    \item variant_A, variant_B: Gene names. Specifically, A and B refer to the merged sets of included variants from each gene.
#'    \item ces_A0: Cancer effect (scaled selection coefficient) of variant A that acts in the absence of variant B.
#'    \item ces_B0: Cancer effect of variant B that acts in the absence of variant A.
#'    \item ces_A_on_B: Cancer effect of variant A that acts when a sample already has variant B.
#'    \item ces_B_on_A: Cancer effect of variant B that acts when a sample already has variant A.
#'    \item p_A_change: P-value of likelihood ratio test (LRT) that informs whether selection for
#'    variant A significantly changes after acquiring variant B. The LRT compares the likelihood of
#'    the full epistatic model to that of a reduced model in which ces_A0 and ces_A_on_B are set
#'    equal. The p-value is the probability, under the reduced model, of the likelihood ratio being
#'    greater than or equal to the ratio observed.
#'    \item p_B_change:  P-value of likelihood ratio test (LRT) that informs whether selection for
#'    variant B significantly changes after acquiring variant A. The LRT compares the likelihood of
#'    the full epistatic model to that of a reduced model in which ces_B0 and ces_B_on_A are set
#'    equal. The p-value is the probability, under the reduced model, of the likelihood ratio being
#'    greater than or equal to the ratio observed.
#'   \item p_epistasis: P-value of likelihood ratio test that informs whether the epistatic model
#'   better explains the mutation data than a non-epistatic model in which selection for mutations
#'   in each gene are independent of the mutation status in the other gene. Quite often, p_epistasis
#'   will suggest a significant epistatic effect even though p_A_change and p_B_change do not
#'   suggest significant changes in selection for either gene individually. This is because the
#'   degree of co-occurrence can often be explained equally well by a strong change in selection for
#'   either gene.
#'   \item expected_nAB_epistasis: The expected number of samples with both A and B mutated under the fitted epistatic model.
#'   Typically, this will be very close to the actual number of AB samples (nAB).
#'   \item expected_nAB_null: The expected number of samples with both A and B mutated under a no-epistasis model.
#'   \item AB_epistatic_ratio: The ratio \code{expected_nAB_epistasis/expected_nAB_null}. Useful to gauge the overall
#'   impact of epistatic interactions on the co-occurrence of variants A and B. Since the expectations take mutation rates into account,
#'   this ratio is a better indicator than the relative frequencies of A0, B0, AB, 00 in the data set.
#'   \item nA0, nB0, nAB, n00: Number of (included) samples with mutations in just A, just B, both A and B, and neither.
#'   \item ces_A_null, ces_B_null: Cancer effects of A and B when estimated independently; that is,
#'   effects under a no-epistasis model.
#'  }
#' @export
ces_gene_epistasis = function(cesa = NULL, genes = NULL, variants = NULL,
                              samples = character(), run_name = "auto", cores = 1, conf = .95,
                              return_fit = FALSE)
{
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  
  samples = select_samples(cesa, samples)
  num_excluded = cesa@samples[, .N] - samples[, .N]
  if(num_excluded > 0) {
    pretty_message(paste0("Note that ", num_excluded, " samples are being excluded from selection inference."))
  }
  
  cesa = copy_cesa(cesa)
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
  
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }

  # genes can be list of character pairs or a character vector, from which all possible pairs will be taken
  gene_pairs = NULL
  if(is(genes, "list") && all(sapply(genes, is.character)) && all(sapply(genes, length) ==2)) {
    gene_pairs = unique(genes)
    genes = unique(unlist(genes))
  } else if (is(genes, "character")) {
    if (length(genes) < 2) {
      stop("Supply at least two genes to analyze.")
    }
  } else{
    stop("genes should be character (at least two gene names) or a list of gene pairs.")
  }
  
  invalid_genes = setdiff(genes, .ces_ref_data[[cesa@ref_key]]$gene_names)
  num_invalid = length(invalid_genes)
  if(num_invalid > 0) {
    if(num_invalid > 19) {
      invalid_genes = c(invalid_genes[1:15], paste0("and ", num_invalid - 15, " more"))
    }
    stop("Input genes are not present in reference data: ", paste(invalid_genes, collapse = ", "))
  } 
  
  maf = cesa@maf[variant_type == "snv"]
  
  if(is.null(variants)) {
    variants = 'recurrent'
    pretty_message(paste0("As of cancereffectsizeR 2.4.0, there is a \"variants\" argument that lets you control ",
      "which variants from each gene are used in epistasis inference. Defaulting to \"recurrent\" strategy to match ",
      "previous behavior. See docs for other options."))
  }
  if(is.character(variants) && length(variants) == 1) {
    if(variants == 'recurrent') {
      variants_to_use = select_variants(cesa = cesa, genes = genes, min_freq = 2)
    } else if (variants == 'nonsilent') {
      # Use variants that appear in MAF data and that are nonsynonymous coding or at essential
      # splice sites (coding or not)
      variants_to_use = select_variants(cesa, genes = genes, min_freq = 1)
      variants_to_use = variants_to_use[essential_splice == TRUE | 
                                          (variant_type == 'aac' & (aa_ref != aa_alt | essential_splice == TRUE))]
    } else {
      stop("variants should be \"recurrent\", \"nonsilent\", or a data.table.")
    }
  } else if (is.data.table(variants)) {
    if (! 'variant_id' %in% names(variants)) {
      stop("Input variants table lacks variant_id column")
    }
    message("Verifying input variant table....", appendLF = F)
    variants_to_use = select_variants(cesa = cesa, variant_ids = variants$variant_id)
    if(attr(variants_to_use, 'nonoverlapping') == FALSE) {
      msg = paste0("Input variants table contains overlapping variant records (e.g., amino acid substitutions that ",
           "contain the same constituent SNVs). To avoid confusion, ",
           "all variants in the inference should be nonoverlapping.")
      stop(pretty_message(msg, emit = F))
    }
    message(" Done.")
  } else if (is.null(variants)) {
    stop("Argument variants should be set to \"recurrent\", \"nonsilent\", or a data.table of variants to include for each gene.")
  } else {
    stop("variants should be \"recurrent\", \"nonsilent\", or a data.table.")
  }
	
	
	# use "all_genes" column because SNVs might have multiple genes, and only one appears in the "gene" column
	genes_present = unique(variants_to_use$gene)
	genes_to_analyze = intersect(genes_present, genes)
	num_passing_genes = length(genes_to_analyze)
  if(num_passing_genes == 0) {
    stop("None of your requested genes have any eligible variants to include in selection inference. ",
         "Changing the \"variants\" argument might help.")
  }
	if(num_passing_genes != length(genes)) {
	  num_missing = length(genes) - num_passing_genes
	  pretty_message(paste0(num_missing, " of your input genes had no eligible variants in input, so will not be tested. ",
	                        "You may want to verify that your \"variants\" input is what you intended."))
	}

	if (num_passing_genes == 1) {
	  msg = paste0("Only 1 requested gene (", genes_to_analyze, ") has eligible variants  in input, so epistasis can't be tested. ",
	               "You may want to verify that your \"variants\" input is what you intended.")
	  stop(pretty_message(msg, emit = F))
	}
	
	if(is.null(gene_pairs)) {
	  gene_pairs = utils::combn(sort(genes_to_analyze), 2, simplify = F)
	} else {
	  gene_pairs = gene_pairs[sapply(gene_pairs, function(x) all(x %in% genes_to_analyze))]
	}
	
	gr_table = as.data.table(.ces_ref_data[[cesa@ref_key]]$gr_genes)
	# Handle pid-based gr_genes
	if('gene' %in% names(gr_table)) {
	  gr_table[, names := gene]
	  gr_table[, gene := NULL]
	}
	setkey(gr_table, "names")
	passing_gene_ranges = gr_table[genes_to_analyze]
	setkey(passing_gene_ranges, "seqnames", "start", "end")
  
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
    gene_pairs = gene_pairs[-bad_pairs]
    if (length(gene_pairs) == 0) {
      stop("There are no remaining gene pairs to test for epistasis.")
    }
  }
  
  # use separate CompoundVariantSets for each pair to avoid possible gene-overlap issues
  setkey(variants_to_use, 'gene')
  results = pbapply::pblapply(X = gene_pairs, FUN = 
    function(x) {
      variant_ids = setNames(list(variants_to_use[x[1], variant_id], variants_to_use[x[2], variant_id]), x)
      comp = CompoundVariantSet(cesa, variant_ids)
      tmp = pairwise_variant_epistasis(cesa = cesa, variant_pair = c(1, 2), samples = samples, compound_variants = comp,
                                           conf = conf)
    }, cl = cores
  )
  fits = lapply(results, '[[', 'fit')
  results = rbindlist(lapply(results, '[[', 'summary'))
  
  if(return_fit) {
    setattr(results, 'fit', fits)
  }
  
  if (results[(nA0 == 0 | nB0 == 0) & 
              nAB == 0, .N] > 0) {
    pretty_message(paste0("Some gene pairs had no eligible variants in one or both genes of jointly-covering samples. ",
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
#'   supply a CompoundVariantSet (see \code{define_compound_variants()}) to test all pairs
#'   of compound variants in the set.
#' @param samples Which samples to include in inference. Defaults to all samples.
#'   Can be a vector of Unique_Patient_Identifiers, or a data.table containing rows from
#'   the CESAnalysis sample table.
#' @param run_name Optionally, a name to identify the current run.
#' @param conf confidence interval size from 0 to 1 (.95 -> 95\%); NULL skips calculation,
#'   reducing runtime.
#' @param cores number of cores for parallel processing of variant pairs
#' @param return_fit TRUE/FALSE (default FALSE): Embed epistatic model fits for each variant pair in
#'   a "fit" attribute of the epistasis results table. Use \code{attr(my_results, 'fit')} to access
#'   the list of fitted models.
#' @return CESAnalysis with a table of epistatic inferences appended to list \code{[CESAnalysis]$epistasis}. Some column definitions:
#'   \itemize{
#'    \item variant_A, variant_B: Names for the two variants or merged sets of variants in each
#'    epistatic inference. For brevity in the case of merged variant sets, we say that a sample with
#'    any variant in variant set A "has variant A."
#'    \item ces_A0: Cancer effect (scaled selection coefficient) of variant A that acts in the absence of variant B.
#'    \item ces_B0: Cancer effect of variant B that acts in the absence of variant A.
#'    \item ces_A_on_B: Cancer effect of variant A that acts when a sample already has variant B.
#'    \item ces_B_on_A: Cancer effect of variant B that acts when a sample already has variant A.
#'    \item p_A_change: P-value of likelihood ratio test (LRT) that informs whether selection for
#'    variant A significantly changes after acquiring variant B. The LRT compares the likelihood of
#'    the full epistatic model to that of a reduced model in which ces_A0 and ces_A_on_B are set
#'    equal. The p-value is the probability, under the reduced model, of the likelihood ratio being
#'    greater than or equal to the ratio observed.
#'    \item p_B_change:  P-value of likelihood ratio test (LRT) that informs whether selection for
#'    variant B significantly changes after acquiring variant A. The LRT compares the likelihood of
#'    the full epistatic model to that of a reduced model in which ces_B0 and ces_B_on_A are set
#'    equal. The p-value is the probability, under the reduced model, of the likelihood ratio being
#'    greater than or equal to the ratio observed.
#'   \item p_epistasis: P-value of likelihood ratio test that informs whether the epistatic model
#'   better explains the mutation data than a non-epistatic model in which selection for mutations
#'   in each variant are independent of the mutation status of the other variant. Quite often, p_epistasis
#'   will suggest a significant epistatic effect even though p_A_change and p_B_change do not suggest
#'   significant changes in selection for either variant individually. This is because the degree of
#'   co-occurrence can often be explained equally well by a strong change in selection for either variant.
#'   \item expected_nAB_epistasis: The expected number of samples with both A and B mutated under the fitted epistatic model.
#'   Typically, this will be very close to the actual number of AB samples (nAB).
#'   \item expected_nAB_null: The expected number of samples with both A and B mutated under a no-epistasis model.
#'   \item AB_epistatic_ratio: The ratio \code{expected_nAB_epistasis/expected_nAB_null}. Useful to gauge the overall
#'   impact of epistatic interactions on the co-occurrence of variants A and B. Since the expectations take mutation rates into account,
#'   this ratio is a better indicator than the relative frequencies of A0, B0, AB, 00 in the data set.
#'   \item nA0, nB0, nAB, n00: Number of (included) samples with mutations in just A, just B, both A and B, and neither.
#'   \item ces_A_null, ces_B_null: Cancer effects of A and B when estimated independently; that is,
#'   effects under a no-epistasis model.
#'  }
#' @export
ces_epistasis = function(cesa = NULL, variants = NULL, samples = character(), run_name = "auto", cores = 1, 
                         conf = .95, return_fit = FALSE) {
  
  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if(! is.null(conf)) {
    if(! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals)", call. = F)
    }
  }
  
  if (! is.null(variants) && length(variants) < 1) {
    stop("variants is 0-length.", call. = F)
  }
  
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should specify a CESAnalysis object", call. = F)
  }
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  
  samples = select_samples(cesa, samples)
  num_excluded_samples = cesa@samples[, .N] - samples[, .N]
  if(num_excluded_samples > 0) {
    message("Note that ", num_excluded_samples, " samples are being excluded from selection inference.")
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
  
  if(is(variants, "CompoundVariantSet")) {
    index_pairs = utils::combn(1:length(variants), 2, simplify = F)
    results = pbapply::pblapply(X = index_pairs, FUN = pairwise_variant_epistasis, samples = samples, compound_variants = variants, cesa=cesa, conf = conf, cl = cores)
  } else if (is(variants, "list")) {
    setkey(cesa@mutations$amino_acid_change, "aac_id")
    setkey(cesa@mutations$snv, "snv_id")
    
    # Validate that all variant IDs appear in mutation annotations
    aac_names = cesa@mutations$amino_acid_change[, .(aac_id, variant_name)]
    
    notify_name_conversion = F
    for (i in 1:length(variants)) {
      pair = variants[[i]]
      if (! is(pair, "character") || length(pair) != 2) {
        stop("All elements of variants list must be character vectors of length 2 (for two variant IDs)")
      }
      
      for (j in 1:2) {
        variant = pair[j]
        if (cesa@mutations$snv[variant, .N, nomatch = NULL] == 0 && aac_names[variant == aac_id, .N, nomatch = NULL] == 0) {
          record_by_name = aac_names[variant, on = "variant_name", nomatch = NULL]
          num_records = record_by_name[, .N]
          if (num_records == 0) {
            stop("Variant ", variant, " could not be found. Either it's not present in the CESAnalysis, or it's a bad ID.")
          } else if (num_records > 1) {
            stop("More than one variant in the CESAnalysis has the name \"", variant, "\". Use a full variant_id to specify which to test.")
          } else {
            # This should never happen with ces.refset.hg38 1.3 and later, because variant_name is uniquely identifying.
            notify_name_conversion = T
            variants[[i]][j] = record_by_name[, aac_id]
          }
        }
      }
    }
    if (notify_name_conversion) {
      pretty_message("Shorthand variant names were recognized and uniquely paired with cancereffectsizeR's aac_ids.")
    }
    results = pbapply::pblapply(X = variants, FUN = pairwise_variant_epistasis, cesa=cesa, samples = samples, conf = conf, cl = cores)
  } else {
    stop("variants should be of type list or CompoundVariantSet")
  }
  fits = lapply(results, '[[', 'fit')
  results = rbindlist(lapply(results, '[[', 'summary'))
  if(results[(nA0 == 0 | nB0 == 0) & 
             nAB == 0, .N] > 0) {
    pretty_message(paste0("Some variant pairs had prevalence of 0 in one or both variants across jointly-covering input samples. ",
                          "Epistatic selection intensities are all NAs for these pairs."))
  }
  if(return_fit) {
    setattr(results, 'fit', fits)
  }
  
  results = list(results)
  names(results) = run_name
  cesa@epistasis = c(cesa@epistasis, results)
  return(cesa)
  
}

#' Calculate SIs at variant level under pairwise epistasis model
#' @param cesa CESAnalysis
#' @param samples Validated samples data.table (as from select_samples())
#' @param variant_pair 2-length character of variant IDs, or 2-length numeric giving
#'   indices of CompoundVariantSet for the current two compound variants
#' @param compound_variants If testing a pair of compound variants, the CompoundVariantSet defining them
#' @keywords internal
pairwise_variant_epistasis = function(cesa, variant_pair, samples, conf, compound_variants = NULL) {
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

  eligible_tumors = samples[covered_regions %in% joint_coverage, Unique_Patient_Identifier]
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
    v1_rates = rowSums(v1_rates)
    names(v1_rates) = all_rates$Unique_Patient_Identifier
    v2_rates = all_rates[, ..v2_ids]
    v2_rates = rowSums(v2_rates)
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
    tumors_with_v1 = intersect(samples_with(cesa, v1), eligible_tumors)
    tumors_with_v2 = intersect(samples_with(cesa, v2), eligible_tumors)
    tumors_with_both = intersect(tumors_with_v1, tumors_with_v2)
    tumors_just_v1 = setdiff(tumors_with_v1, tumors_with_both)
    tumors_just_v2 = setdiff(tumors_with_v2, tumors_with_both)
    tumors_with_neither = setdiff(eligible_tumors, c(tumors_with_v1, tumors_with_v2))
    
    # 2-item lists: first item is baseline rates for v1; second for v2
    with_just_1 = lapply(as.list(all_rates[tumors_just_v1])[2:3], setNames, tumors_just_v1)
    with_just_2 = lapply(as.list(all_rates[tumors_just_v2])[2:3], setNames, tumors_just_v2)
    with_both = lapply(as.list(all_rates[tumors_with_both])[2:3], setNames, tumors_with_both)
    with_neither = lapply(as.list(all_rates[tumors_with_neither])[2:3], setNames, tumors_with_neither)
  }

  # call factory function to get variant-specific likelihood function
  epistasis_lik_fn = pairwise_epistasis_lik(with_just_1, with_just_2, with_both, with_neither)
  par_init = formals(epistasis_lik_fn)[[1]]
  names(par_init) = bbmle::parnames(epistasis_lik_fn)
  
  # No point of testing epistasis if either variant doesn't appear
  if (length(tumors_with_v1) == 0 || length(tumors_with_v2) == 0) {
    n_total = length(tumors_just_v1) + length(tumors_just_v2) + length(tumors_with_neither) # none have both
    early_output = list(variant_A = v1, variant_B = v2, ces_A0 = NA_real_, ces_B0 = NA_real_,
                        ces_A_on_B = NA_real_, ces_B_on_A = NA_real_, 
                        p_A_change = NA_real_, p_B_change = NA_real_, p_epistasis = NA_real_, 
                        expected_nAB_epistasis = NA_real_, expected_nAB_null = NA_real_,
                        AB_epistatic_ratio = NA_real_,
                        nA0 = length(tumors_just_v1),
                        nB0 = length(tumors_just_v2),
                        nAB = length(tumors_with_both),
                        n00 = length(tumors_with_neither),
                        n_total = n_total,
                        ces_A_null = NA_real_, ces_B_null = NA_real_)
    if (! is.null(conf)) {
     si_names = c('ces_A0', 'ces_B0', 'ces_A_on_B', 'ces_B_on_A')
     ci_high_colname = paste0("ci_high_", conf * 100)
     ci_low_colname = paste0("ci_low_", conf * 100)
     low_colnames = paste(ci_low_colname, si_names, sep = "_")
     high_colnames = paste(ci_high_colname, si_names, sep = "_")
     ci_colnames = unlist(S4Vectors::zipup(low_colnames, high_colnames))
     ci_output = as.list(rep.int(NA_real_, 8))
     names(ci_output) = ci_colnames
     early_output = c(early_output, ci_output) # low/high for 4 parameters
    }
    return(list(summary = early_output, fit = NULL)) # fit list will have a NULL entry
  }
  
  # find optimized selection intensities
  # the selection intensity when some group has 0 variants will be on the lower boundary; will muffle the associated warning
  withCallingHandlers(
    {
      fit = bbmle::mle2(epistasis_lik_fn, method="L-BFGS-B", start = par_init, control = list(ndeps = c(.1, .1, .1, .1)),
                        vecpar = TRUE, lower=1e-3, upper=1e9)
    },
    warning = function(w) {
      if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  
  # Also remove the whole copy of CESAnalysis in there.
  # To-do: probably shouldn't be there in the first place.
  rm('cesa', envir = environment(fit))
  params = bbmle::coef(fit)
  
  
  # Get non-epistatic effects for both v1 and v2
  v1_rates_samples_with = c(with_just_1[[1]], with_both[[1]])
  v1_rates_samples_without = c(with_just_2[[1]], with_neither[[1]])
  fn = sswm_lik(v1_rates_samples_with, v1_rates_samples_without)
  par_init = formals(fn)[[1]]
  names(par_init) = bbmle::parnames(fn)
  withCallingHandlers(
    {
      v1_fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = TRUE, lower=1e-3, upper=1e9)
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
  
  v2_rates_samples_with = c(with_just_2[[2]], with_both[[2]])
  v2_rates_samples_without = c(with_just_1[[2]], with_neither[[2]])
  fn = sswm_lik(v2_rates_samples_with, v2_rates_samples_without)
  par_init = formals(fn)[[1]]
  names(par_init) = bbmle::parnames(fn)
  withCallingHandlers(
    {
      v2_fit = bbmle::mle2(fn, method="L-BFGS-B", start = par_init, vecpar = TRUE, lower=1e-3, upper=1e9)
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
  
  v1_simple_estimate = bbmle::coef(v1_fit) 
  v2_simple_estimate = bbmle::coef(v2_fit)
  ll_check = as.numeric(bbmle::logLik(bbmle::mle2(epistasis_lik_fn, start = list(ces_v1 = v1_simple_estimate, ces_v1_after_v2 = v1_simple_estimate,
                                                 ces_v2 = v2_simple_estimate, ces_v2_after_v1 = v2_simple_estimate),
                                                 control = list(ndeps = c(.1, .1, .1, .1)), vecpar = TRUE, eval.only = T)))
  if(! isTRUE(all.equal(ll_check, as.numeric(bbmle::logLik(v1_fit)) + as.numeric(bbmle::logLik(v2_fit)), tolerance = 1e-6))) {
    msg = paste0('An internal check failed during significance testing on ', v1, '/', v2, '. P-values cannot be reported ', 
                 'for this pair. We would appreciate it if you could report the issue on GitHub.')
    warning(pretty_message(msg, emit = F))
    p_v1 = NA
    p_v2 = NA
    p_epistasis = NA
  } else {
    param_init_v1 = list(ces_v1 = v1_simple_estimate, ces_v2 = 1e3, 
                         ces_v1_after_v2 = v1_simple_estimate, ces_v2_after_v1 = 1e3)
    
    withCallingHandlers(
      {
        refit_v1 = bbmle::mle2(epistasis_lik_fn, method="L-BFGS-B", start = param_init_v1, 
                               control = list(ndeps = c(.1, .1)),
                               vecpar = TRUE, lower = 1e-3, upper= 1e9,
                               fixed = list(ces_v1 = v1_simple_estimate, ces_v1_after_v2 = v1_simple_estimate))
      },
      warning = function(w) {
        if (conditionMessage(w) %ilike% "some parameters are on the boundary") {
          invokeRestart("muffleWarning")
        }
      }
    )
    v1_chisquared = as.numeric(-2 * (bbmle::logLik(refit_v1) - bbmle::logLik(fit)))
    p_v1 =  pchisq(v1_chisquared, df = 2, lower.tail = F)
    
    
    param_init_v2 = list(ces_v1 = 1e3, ces_v2 = v2_simple_estimate,
                         ces_v1_after_v2 = 1e3, ces_v2_after_v1 = v2_simple_estimate)
    
    withCallingHandlers(
      {
        refit_v2 = bbmle::mle2(epistasis_lik_fn, method="L-BFGS-B", start = param_init_v2, control = list(ndeps = c(.1, .1)),
                               vecpar = TRUE, lower = 1e-3, upper = 1e9, 
                               fixed = list(ces_v2 = v2_simple_estimate, ces_v2_after_v1 = v2_simple_estimate))
      },
      warning = function(w) {
        if (conditionMessage(w) %ilike% "some parameters are on the boundary") {
          invokeRestart("muffleWarning")
        }
      }
    )
    v2_chisquared = as.numeric(-2 * (bbmle::logLik(refit_v2) - bbmle::logLik(fit)))
    p_v2 = pchisq(v2_chisquared, df = 2, lower.tail = F)
    
    chisquared = as.numeric(-2 * (bbmle::logLik(v1_fit) + bbmle::logLik(v2_fit) - bbmle::logLik(fit)))
    p_epistasis = pchisq(chisquared, df = 2, lower.tail = F)
  }
  
  # Collect v1, v2 rates in the same sample order
  v1_rates = c(v1_rates_samples_with, v1_rates_samples_without)[eligible_tumors]
  v2_rates = c(v2_rates_samples_with, v2_rates_samples_without)[eligible_tumors]
  
  # parameter estimates from epistatic model
  ces_v1 = params[1]
  ces_v2 = params[2]
  ces_v1_after_v2 = params[3]
  ces_v2_after_v1 = params[4]
  
  # Variables named for parallelism with pairwise_epistasis_lik().
  A = params[1] * v1_rates
  B = params[2] * v2_rates
  A_on_B = params[3] * v1_rates
  B_on_A = params[4] * v2_rates
  
  # See pairwise_epistasis_lik()
  p_wt = exp(-1 * (A+B))
  p_A = (A / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))
  p_B = (B / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))
  p_AB = 1 - p_wt - p_A - p_B
  expected_nAB = sum(p_AB)
  
  # Under no-epistasis model, p(AB) = p(A) * p(B)
  A = v1_simple_estimate * v1_rates
  B = v2_simple_estimate * v2_rates
  expected_nAB_null = sum((1 - exp(-1 * A)) * (1 - exp(-1 * B)))
  
  variant_ep_results = list(variant_A = v1, variant_B = v2, ces_A0 = params[1], ces_B0 = params[2],
                            ces_A_on_B = params[3], ces_B_on_A = params[4], 
                            p_A_change = p_v1, p_B_change = p_v2, p_epistasis = p_epistasis, 
                            expected_nAB_epistasis = expected_nAB, expected_nAB_null = expected_nAB_null,
                            AB_epistatic_ratio = expected_nAB / expected_nAB_null,
                            nA0 = length(tumors_just_v1),
                            nB0 = length(tumors_just_v2),
                            nAB = length(tumors_with_both),
                            n00 = length(tumors_with_neither),
                            n_total = length(v1_rates),
                            ces_A_null = v1_simple_estimate, ces_B_null = v2_simple_estimate)
  if(! is.null(conf)) {
    bbmle::parnames(epistasis_lik_fn) = c('ces_A0', 'ces_B0', 'ces_A_on_B', 'ces_B_on_A')
    variant_ep_results = c(variant_ep_results, univariate_si_conf_ints(fit, epistasis_lik_fn, .001, 1e9, conf))
  }
  return(list(summary = variant_ep_results, fit = fit))
}

