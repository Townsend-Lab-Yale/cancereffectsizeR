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
#'   patient_ids, or a data.table containing rows from the CESAnalysis sample table.
#'   Samples that do not have coverage at all variant sites in a given inference will be set aside.
#' @param run_name Optionally, a name to identify the current run.
#' @param model Set to "default"  to use built-in
#'   model of epistatic selection, or supply a custom function factory (see details).
#' @param lik_args Extra arguments, given as a list, to pass to custom likelihood functions.
#' @param optimizer_args optimizer_args List of arguments to pass to the optimizer (bbmle::mle2).
#' @param pval_calc_fn For use with custom models; optional. A function that takes an epistasis
#'   model fit as input and returns p-values and other descriptives.
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
ces_gene_epistasis = function(cesa = NULL, 
                              genes = NULL, 
                              variants = NULL,
                              samples = character(), 
                              model = "default",
                              run_name = "auto", 
                              cores = 1, 
                              conf = .95,
                              lik_args = list(),
                              optimizer_args = list(),
                              pval_calc_fn = NULL,
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
  setkey(cesa@mutations$sbs, "sbs_id")
  
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
  
  maf = cesa@maf[variant_type == "sbs"]
  
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
    
    # Check if SBS are shared between AAC
    aac_to_check = variants_to_use[variant_type == 'aac', variant_id]
    nonoverlapping = cesa@mutations$aac_sbs_key[aac_to_check, .(is_unique = uniqueN(aac_id) == 1), on = 'aac_id', by = 'sbs_id'][, all(is_unique)]
    if(! nonoverlapping) {
      msg = paste0("The input variants table contains overlapping variants (typically, amino acid substitutions with different protein IDs that ",
           "are caused by the same SBS variants). To avoid confusion, ",
           "all variants in the inference should be nonoverlapping.")
      stop(pretty_message(msg, emit = F))
    }
    message(" Done.")
  } else if (is.null(variants)) {
    stop("Argument variants should be set to \"recurrent\", \"nonsilent\", or a data.table of variants to include for each gene.")
  } else {
    stop("variants should be \"recurrent\", \"nonsilent\", or a data.table.")
  }
	
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
      tmp = pairwise_variant_epistasis(cesa = cesa, variant_pair = c(1, 2), samples = samples, compound_variants = comp, conf = conf, 
                                       model = model, lik_args = lik_args, pval_calc_fn = pval_calc_fn,
                                       optimizer_args = optimizer_args)
    }, cl = cores
  )
  fits = lapply(results, '[[', 'fit')
  results = rbindlist(lapply(results, '[[', 'summary'), fill = TRUE)
  
  if(return_fit) {
    fits = lapply(fits, function(x) {
      x@call.orig = call('[not shown]')
      parent.env(environment(x)) = emptyenv()
      return(x)
    })
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
#'   Can be a vector of patient_ids, or a data.table containing rows from
#'   the CESAnalysis sample table.
#' @param run_name Optionally, a name to identify the current run.
#' @param conf confidence interval size from 0 to 1 (.95 -> 95\%); NULL skips calculation,
#'   reducing runtime.
#' @param model Set to "default"  to use built-in
#'   model of epistatic selection, or supply a custom function factory (see details).
#' @param lik_args Extra arguments, given as a list, to pass to custom likelihood functions.
#' @param optimizer_args List of arguments to pass to the optimizer (bbmle::mle2).
#' @param pval_calc_fn For use with custom models; optional. A function that takes an epistasis
#'   model fit as input and returns p-values and other descriptives.
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
                         model = 'default', lik_args = list(), pval_calc_fn = NULL, optimizer_args = list(),
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
    results = pbapply::pblapply(X = index_pairs, FUN = pairwise_variant_epistasis, samples = samples, compound_variants = variants, cesa=cesa, conf = conf, cl = cores, model = model, lik_args = lik_args)
  } else if (is(variants, "list")) {
    setkey(cesa@mutations$amino_acid_change, "aac_id")
    setkey(cesa@mutations$sbs, "sbs_id")
    
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
        if (cesa@mutations$sbs[variant, .N, nomatch = NULL] == 0 && aac_names[variant == aac_id, .N, nomatch = NULL] == 0) {
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
    results = pbapply::pblapply(X = variants, FUN = pairwise_variant_epistasis, cesa=cesa, samples = samples, conf = conf, 
                                cl = cores, model = model, lik_args = lik_args, optimizer_args = optimizer_args,
                                pval_calc_fn = pval_calc_fn)
  } else {
    stop("variants should be of type list or CompoundVariantSet")
  }
  fits = lapply(results, '[[', 'fit')
  results = rbindlist(lapply(results, '[[', 'summary'), fill = TRUE)
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


