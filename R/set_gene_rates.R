#' Assign pre-calculated neutral gene mutation rates
#'
#' This function allows you to assign gene mutation rates calculated however you like to
#' samples in your CESAnalysis. Rates can be assigned to all samples or to samples
#' with the specified sample_group labels (see \code{?CESAnalysis}).
#' 
#' Provide rates in a gene_rates data.table with two columns: gene (character) and rate
#' (numeric on \code{[0, 1]}). Gene names must match the CESAnalysis reference data; to
#' see all genes in your CESAnalysis \code{cesa}, call
#' \code{unique(cesa$reference_data$gene_ranges$names)}. If you don't want to supply rates
#' for every gene, set \code{missing_genes_take_nearest = T} to have each missing gene
#' take the rate of the nearest non-missing gene.
#' 
#' 
#' @param cesa CESAnalysis object
#' @param gene_rates A two-column data.table with genes in column 1 and rate in column 2
#' @param sample_group Character vector giving sample group(s) that the rates apply to. If
#'   unspecified, the rates are assigned to all samples.
#' @param missing_genes_take_nearest Set to TRUE
#' to have each gene missing from gene_rates take the rate of the nearest non-missing gene.
#' @export
set_gene_rates = function(cesa = NULL, gene_rates = NULL, sample_group = NULL, missing_genes_take_nearest = FALSE) {
  stopifnot(is(cesa, "CESAnalysis"),
            is(gene_rates, "data.table"),
            is(missing_genes_take_nearest, "logical"),
            length(missing_genes_take_nearest) == 1)
  if(is.null(sample_group)) {
    sample_group = cesa@groups
  }
  if(! is(sample_group, "character")) {
    stop("sample_group should be character")
  }
  sample_group = unique(sample_group)
  if (! all(sample_group %in% cesa@groups)) {
    possible_groups = setdiff(cesa@groups, "stageless") # historical
    if (length(possible_groups) == 0) {
      stop("The CESAnalysis has no user-defined sample groups, so you can't use the sample_group parameter.")
    }
    stop("Unrecognized sample groups. Those defined in the CESAnalysis are ", paste(possible_groups, sep = ", "), ".")
  }
  
  curr_samples = cesa@samples[group %in% sample_group, "Unique_Patient_Identifier"]
  if (length(curr_samples) == 0) {
    stop("No samples are in the specified group(s)")
  }
  
  # Don't allow rates to be overwritten (too confusing for now to risk overlapping sample group specification in sequential calls)
  if (! is.null(cesa@samples$gene_rate_grp)) {
    if (! all(is.na(cesa@samples[group %in% sample_group, gene_rate_grp]))) {
      stop("Gene mutation rates have already been calculated/assigned for some or all samples specified.")
    }
  }
  
  # validate gene_rates
  if(ncol(gene_rates) != 2) {
    stop("gene_rates should be a two-column table giving genes and their associated rates.")
  }
  colnames(gene_rates) = c("gene", "rate") # will rename "rate" to current group after more validation
  if (! is(gene_rates$gene, "character") | ! is(gene_rates$rate, "numeric")) {
    stop("In gene_rates, gene and rate should be character and numeric, respectively")
  }
  if(anyNA(gene_rates)) {
    stop("There are NA entries in gene_rates.")
  }
  if(any(duplicated(gene_rates$gene))) {
    stop("There are duplicate gene entries in gene_rates.")
  }
  
  if(gene_rates[, any(rate < 0 | rate > 1)]) {
    stop("All rates must be on [0, 1]")
  }
  
  # Load ref set genes and compare to gene_rates
  gr_genes = get_ref_data(cesa, "gr_genes")
  all_genes = unique(gr_genes$names)
  
  if(length(setdiff(gene_rates$gene, all_genes)) > 0) {
    stop("There are gene names in gene_rates that are absent from the CESAnalysis's reference data set")
  }
  
  missing_genes = setdiff(all_genes, gene_rates$gene)
  if (length(missing_genes) > 0) {
    if (missing_genes_take_nearest == FALSE) {
      stop("Not all reference genes present in gene_rates. To have each gene take the rate of the nearest ",
           "non-missing gene, re-run with missing_genes_take_nearest = TRUE")
    } else {
      gene_rates = fill_missing_genes(gr_genes, gene_rates)
    }
  }
  if (ncol(cesa@mutrates) == 0) {
    curr_rate_group = 'rate_grp_1'
  } else {
    curr_rate_group = paste0('rate_grp_', ncol(cesa@mutrates))
  }
  
  setnames(gene_rates, "rate", curr_rate_group)
  
  cesa@samples = cesa@samples[group %in% sample_group, gene_rate_grp := curr_rate_group]
  
  if (ncol(cesa@mutrates) > 0) {
    cesa@mutrates = cesa@mutrates[gene_rates, on = "gene"]
  } else {
    cesa@mutrates = gene_rates
  }
  
  if (! any(is.na(cesa@samples$gene_rate_grp))) {
    cesa@advanced$gene_rates_done = T
  }
  
  return(cesa)
}