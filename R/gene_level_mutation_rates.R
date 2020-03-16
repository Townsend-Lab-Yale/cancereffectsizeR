#' Gene-level mutation rates
#'
#' This function calculates gene-level neutral mutation rates based on counts
#' of nonsynonymous and synonymous mutations per gene under the dNdScv package's model, 
#' as described in Martincorena et al. (https://doi.org/10.1016/j.cell.2017.09.042).
#'
#' @param cesa CESAnalysis object
#' @param covariate_file "default" uses dNdScv hg19 default covariates, NULL runs without covariates, or give a custom prcomp RData file, or specify 
#' one of these for hg19 tissue-specific covariate data: "bladder_pca"  "breast_pca" "cesc_pca" "colon_pca" "esca_pca" "gbm_pca" "hnsc_pca" "kidney_pca" "lihc_pca" 
#' "lung_pca" "ov_pca" "pancreas_pca" "prostate_pca" "rectum_pca" "skcm_pca"  "stomach_pca"  "thca_pca" "ucec_pca"
#' @param save_all_dndscv_output default false; when true, saves all dndscv output, not just what's needed by CES (will make object size very large)
#' @return CESAnalysis object with gene-level mutation rates calculated
#' @export
# don't change this function at all without being sure you're not messing up tests
gene_level_mutation_rates <- function(cesa, covariate_file = 'default', save_all_dndscv_output = FALSE){
  RefCDS = get_genome_data(cesa, "RefCDS")
  # hacky way of forcing an object of name gr_genes into the dndscv::dndscv function environment,
  # since the object is required by dndscv but there's no argument to supply your own copy of it
  our_env = new.env(parent = environment(dndscv::dndscv))
  our_env$gr_genes = get_genome_data(cesa, "gr_genes")
  f = dndscv::dndscv
  environment(f) = our_env
  dndscv_input = dndscv_preprocess(cesa = cesa, covariate_file = covariate_file)
  dndscv_input = lapply(dndscv_input, function(x) { x[["refdb"]] = RefCDS; return(x)})
  message("Running dNdScv...")
  withCallingHandlers(
    {
      dndscv_raw_output = lapply(dndscv_input, function(x) do.call(f, x))
    }, error = function(e) {
      if (startsWith(conditionMessage(e), "bad 'file' argument"))  {
        stop("You need to update dNdScv. Try running \"devtools::update_packages(packages = \"dndscv\")\".")
      }
    }
  )
  cesa = dndscv_postprocess(cesa = cesa, dndscv_raw_output = dndscv_raw_output, save_all_dndscv_output = save_all_dndscv_output)
  cesa@status[["gene mutation rates"]] = paste0("Calculated with dNdScv using \"", covariate_file, "\" covariates")
  return(cesa)
}


#' Internal function to prepare for running dNdScv
dndscv_preprocess = function(cesa, covariate_file = "default") {
  if(is.null(covariate_file)){
    message("Warning: Calculating gene mutation rates with no covariate data (supply covariates if available).")
    warning("Calculating gene mutation rates with no covariate data (covariates should be supplied if available)")
  } else if (is(covariate_file, "character") && covariate_file[1] == "default") {
    if(names(cesa@genome_data_dir) == "hg19") {
      message("Loading dNdScv default covariates for hg19...")
      data("covariates_hg19",package = "dndscv", envir = environment())
      genes_in_pca <- rownames(covs)
      cv = "hg19"
    } else {
      stop("There is no default covariates data for this genome build, so you'll need to supply your own or run without.")
    }
  } else {
    this_cov_pca <- get(covariate_file) # To-do: clean this up
    cv = this_cov_pca$rotation
    genes_in_pca = rownames(cv)
  }

  # select MAF records to use for gene mutation rate calculation (via dNdScv)
  # for now, records from TGS samples are kept out; in the future, we could check out dNdScv's panel sequencing features
  dndscv_samples = cesa@samples[coverage == "exome" | coverage == "genome", Unique_Patient_Identifier]
  
  dndscv_maf = cesa@maf[Unique_Patient_Identifier %in% dndscv_samples,]

  dndscv_input = list()
  for (progression in cesa@progressions) {
    current_subset_tumors = cesa@samples[progression_name == progression, Unique_Patient_Identifier]
    mutations = dndscv_maf[dndscv_maf$Unique_Patient_Identifier %in% current_subset_tumors,]
    dndscv_input[[progression]] = list(mutations = mutations, gene_list = genes_in_pca, cv = cv)
  }
  return(dndscv_input)
}

#' Internal function to calculate gene-level mutation rates from dNdScv output
dndscv_postprocess = function(cesa, dndscv_raw_output, save_all_dndscv_output = FALSE) {
  # load RefCDS data if it's not already in the environment
  if(! "RefCDS" %in% ls()) {
    RefCDS = get_genome_data(cesa, "RefCDS")
  }
  dndscv_out_list = dndscv_raw_output
  names(dndscv_out_list) = cesa@progressions
   # Get RefCDS data on number of synonymous mutations possible at each site
  # Per dNdScv docs, L matrices list "number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context"
  num_syn = sapply(RefCDS, function(x) colSums(x$L)[1])
  names(num_syn) = sapply(RefCDS, function(x) x$gene_name)

  dndscv_genes = dndscv_out_list[[1]]$genemuts$gene_name # dndscv uses same set of genes for each stage
  num_syn = num_syn[names(num_syn) %in% dndscv_genes]

  mutrates_list <- vector(mode = "list",length = length(cesa@progressions))
  names(mutrates_list) <- cesa@progressions


  message("Using dNdScv output to calculate gene-level mutation rates...")
  for(this_subset in 1:length(mutrates_list)){
    dndscv_output = dndscv_out_list[[this_subset]] 
    number_of_tumors_in_this_subset <- length(unique(dndscv_output$annotmuts$sampleID))
    if(dndscv_output$nbreg$theta>1){

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
    names(mutrates_vec) <- dndscv_output$genemuts$gene_name
    mutrates_list[[this_subset]] = mutrates_vec
  }

  # keep just the main gene-level selection output from dNdScv, unless user wanted everything
  # currently also need annotmuts for annotate_gene_maf
  if(! save_all_dndscv_output) {
    for (i in 1:length(dndscv_out_list)) {
      # filter out genes with 0 mutations (to keep object size small, mainly for dev purposes)
      sel_cv = dndscv_out_list[[i]]$sel_cv
      #no_dndscv_mutations = (sel_cv$n_syn == 0 & sel_cv$n_mis == 0 & sel_cv$n_non == 0 & sel_cv$n_spl == 0)
      dndscv_out_list[[i]] = list(sel_cv = sel_cv, annotmuts = dndscv_out_list[[i]]$annotmuts)
    }
  }
  cesa@mutrates_list = mutrates_list
  cesa@dndscv_out_list = dndscv_out_list
  return(cesa)
}

