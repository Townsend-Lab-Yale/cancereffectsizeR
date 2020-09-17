#' Use dNdScv with tissue-specific covariates to calculate gene-level mutation rates
#'
#' This function calculates gene-level neutral mutation rates based on counts
#' of nonsynonymous and synonymous mutations per gene under the dNdScv package's model, 
#' as described in Martincorena et al. (https://doi.org/10.1016/j.cell.2017.09.042).
#'
#' @param cesa CESAnalysis object
#' @param covariates Name of gene mutation rate covariates to use (run
#'   list_ces_covariates() to see choices). For hg19 only, you can also use "default" for
#'   dNdScv's non-tissue-specific covariates. If no covariates are available, set NULL to
#'   run without.
#' @param save_all_dndscv_output default false; when true, saves all dndscv output, not
#'   just what's needed by cancereffectsizeR (will make object size very large)
#' @return CESAnalysis object with gene-level mutation rates calculated
#' @export
# don't change this function at all without being sure you're not messing up tests
gene_mutation_rates <- function(cesa, covariates = NULL, save_all_dndscv_output = FALSE){
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa expected to be a CESAnalysis object", call. = F)
  }
  cesa = update_cesa_history(cesa, match.call())
  if (! cesa@ref_key %in% ls(.ces_ref_data)) {
    preload_ref_data(cesa@ref_data_dir)
  }
  RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS
  gr_genes = .ces_ref_data[[cesa@ref_key]]$gr_genes
  
  # hacky way of forcing an object of name gr_genes into the dndscv::dndscv function environment,
  # since the object is required by dndscv but there's no argument to supply your own copy of it
  our_env = new.env(parent = environment(dndscv::dndscv))
  our_env$gr_genes = gr_genes
  f = dndscv::dndscv
  environment(f) = our_env
  dndscv_input = dndscv_preprocess(cesa = cesa, covariates = covariates)
  dndscv_input = lapply(dndscv_input, function(x) { x[["refdb"]] = RefCDS; return(x)})
  message("Running dNdScv...")
  withCallingHandlers(
    {
      dndscv_raw_output = lapply(dndscv_input, function(x) do.call(f, x))
    }, error = function(e) {
      if (startsWith(conditionMessage(e), "bad 'file' argument"))  {
        stop("You need to update dNdScv. Try running \"devtools::update_packages(packages = \"dndscv\")\".")
      }
    }, warning = function(w) {
      dndscv_msg = conditionMessage(w)
      # Recurrent mutations are presumed real in CES
      # As of 09/03/20, dNdScv's MNV check is flawed (and CES has already handled them in load_maf())
      if (startsWith(dndscv_msg, "Same mutations observed in different sampleIDs") || 
          startsWith(dndscv_msg, "Mutations observed in contiguous sites")) {
        invokeRestart("muffleWarning")
      }
    }
  )
  cesa = dndscv_postprocess(cesa = cesa, dndscv_raw_output = dndscv_raw_output, save_all_dndscv_output = save_all_dndscv_output)
  return(cesa)
}


#' Internal function to prepare for running dNdScv
#' @keywords internal
dndscv_preprocess = function(cesa, covariates = "default") {
  groups_with_data = cesa@samples[coverage %in% c("exome", "genome"), unique(group)]
  groups_without_data = setdiff(cesa@groups, groups_with_data)
  if(length(groups_without_data) > 0) {
    stop(paste0("Cannot run dNdScv by sample group because the following groups have no whole-exome/whole-genome samples,\n",
                "which are needed for gene mutation rate calculation: ", paste(groups_without_data, collapse = ", ")))
  }
  if(is.null(covariates)){
    cv = NULL
    genes_in_pca = NULL
    warning("Calculating gene mutation rates with no covariate data; stop and re-run with covariates if available.\n",
            "(Check with list_ces_covariates(); for hg19 only, you can also specify \"default\" for dNdScv default\n",
            "covariates.)", call. = F, immediate. = T)
  } else if (is(covariates, "character") && covariates[1] == "default") {
    if(cesa@ref_key == "ces_hg19_v1") {
      message("Loading dNdScv default covariates for hg19 (stop and re-run with tissue-specific covariates if available)...")
      data("covariates_hg19", package = "dndscv", envir = environment())
      genes_in_pca <- rownames(covs)
      cv = "hg19"
    } else {
      stop("There is no default covariates data for this genome build, so you'll need to supply your own\n",
           "or run without by setting covariates = NULL.", call. = F)
    }
  } else if (! is(covariates, "character") || length(covariates) != 1) {
    stop("covariates expected to be 1-length character. Check available covariates with list_ces_covariates()")
  } else {
    covariates = paste0("covariates/", covariates)
    if(! check_for_ref_data(cesa, covariates)) {
      stop("Covariates could not be found. Check available covariates with list_ces_covariates().")
    }
    this_cov_pca <- get_ref_data(cesa, covariates) 
    cv = this_cov_pca$rotation
    genes_in_pca = rownames(cv)
  }

  # select MAF records to use for gene mutation rate calculation (via dNdScv)
  # for now, records from TGS samples are kept out; in the future, we could check out dNdScv's panel sequencing features
  dndscv_samples = cesa@samples[coverage == "exome" | coverage == "genome", Unique_Patient_Identifier]
  
  dndscv_maf = cesa@maf[Unique_Patient_Identifier %in% dndscv_samples,]

  dndscv_input = list()
  for (curr_group in cesa@groups) {
    current_subset_tumors = cesa@samples[group == curr_group, Unique_Patient_Identifier]
    mutations = dndscv_maf[dndscv_maf$Unique_Patient_Identifier %in% current_subset_tumors,]
    if(nrow(mutations) == 0) {
      stop(paste0("Can't run dNdScv because sample group ", curr_group, " has no usable SNV mutations."))
    }
    dndscv_input[[curr_group]] = list(mutations = mutations, gene_list = genes_in_pca, cv = cv)
  }
  return(dndscv_input)
}

#' Internal function to calculate gene-level mutation rates from dNdScv output
#' @keywords internal
dndscv_postprocess = function(cesa, dndscv_raw_output, save_all_dndscv_output = FALSE) {
  RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS
  gr_genes = .ces_ref_data[[cesa@ref_key]]$gr_genes
  dndscv_out_list = dndscv_raw_output
  names(dndscv_out_list) = cesa@groups
   # Get RefCDS data on number of synonymous mutations possible at each site
  # Per dNdScv docs, L matrices list "number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context"
  num_syn = sapply(RefCDS, function(x) colSums(x$L)[1])
  names(num_syn) = sapply(RefCDS, function(x) x$gene_name)

  dndscv_genes = dndscv_out_list[[1]]$genemuts$gene_name # dndscv uses same set of genes for each stage
  num_syn = num_syn[names(num_syn) %in% dndscv_genes]

  mutrates_list <- vector(mode = "list",length = length(cesa@groups))
  names(mutrates_list) <- cesa@groups

  message("Using dNdScv output to calculate gene-level mutation rates...")
  for(this_subset in 1:length(mutrates_list)){
    dndscv_output = dndscv_out_list[[this_subset]] 
    number_of_tumors_in_this_subset <- length(unique(dndscv_output$annotmuts$sampleID))
    if(dndscv_output$nbreg$theta>1){
      # see page e4 of dNdScv paper (Martincorena 2017, Cell)
      mutrates_vec <- ((dndscv_output$genemuts$n_syn +
                          dndscv_output$nbreg$theta -
                          1) /
                         (1 +
                            (dndscv_output$nbreg$theta /
                               dndscv_output$genemuts$exp_syn_cv)
                         )
      ) /
        num_syn /
        number_of_tumors_in_this_subset

    } else{
      mutrates_vec <- rep(NA,length(dndscv_output$genemuts$exp_syn_cv))
      syn_sites <- num_syn
      for(i in 1:length(mutrates_vec)){
        mutrates_vec[i] <-  max(dndscv_output$genemuts$exp_syn_cv[i],  ((dndscv_output$genemuts$n_syn[i] +
                                                                           dndscv_output$nbreg$theta -
                                                                           1) /
                                                                          (1 +
                                                                             (dndscv_output$nbreg$theta /
                                                                                dndscv_output$genemuts$exp_syn_cv[i])
                                                                          ))) /
          num_syn[i] /
          number_of_tumors_in_this_subset
      }
    }
    mutrates_list[[this_subset]] = mutrates_vec
  }
  mutrates_dt = as.data.table(mutrates_list)
  mutrates_dt[, gene := dndscv_out_list[[1]]$genemuts$gene_name] # all runs of dNdScv have same the genes in the same order
  setcolorder(mutrates_dt, "gene")
  
  # Genes in gr_genes that are not present in the covariates data (a few are missing, usually),
  # won't get rates calcualted by dNdScv. Here, we assign them the rate of the nearest gene,
  # as measured by center-to-center distance.
  missing_genes = setdiff(GenomicRanges::mcols(gr_genes)["names"][,1], mutrates_dt$gene)
  
  tmp = as.data.table(gr_genes)
  gene_intervals = tmp[, .(chr = as.character(seqnames[1]), start = min(start), end = max(end)), by = "names"]
  setnames(gene_intervals, "names", "gene")
  gene_intervals[, center := start + trunc((end - start)/2)]
  gene_intervals = gene_intervals[order(chr, center)]
  present_genes = gene_intervals[! gene %in% missing_genes]
  not_present = gene_intervals[gene %in% missing_genes]
  
  # find nearest genes by chr/center to the missing genes from the present genes,
  # and eliminate the rare tie by taking the first match for each gene
  nearest_genes = present_genes[not_present, on = c("chr", "center"), roll = "nearest"]
  nearest_genes = nearest_genes[! duplicated(i.gene)]
  
  setkey(mutrates_dt, "gene")
  nearest_rates = mutrates_dt[nearest_genes$gene, -"gene"]
  missing_rates = cbind(gene = nearest_genes$i.gene, nearest_rates)
  mutrates_dt = rbind(mutrates_dt, missing_rates)

  # keep just the main gene-level selection output from dNdScv, unless user wanted everything
  if(! save_all_dndscv_output) {
    for (i in 1:length(dndscv_out_list)) {
      dndscv_out_list[[i]] = dndscv_out_list[[i]]$sel_cv
    }
  }
  cesa@mutrates = mutrates_dt
  cesa@dndscv_out_list = dndscv_out_list
  return(cesa)
}

