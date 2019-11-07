#' Gene-level mutation rate function
#'
#' This function calculates the maximum likelihood estimate for the expected
#' number of synonymous mutations in a gene under the dNdScv model, as described within
#' Martincorena, I., Raine, K. M., Gerstung, M., Dawson, K. J., Haase, K., Van Loo, P., … Campbell, P. J.
#' (2017). Universal Patterns of Selection in Cancer and Somatic Tissues. Cell. https://doi.org/10.1016/j.cell.2017.09.042
#' and the associated R package dndscv .
#'
#' @param cesa CESAnalysis object
#' @param covariate_file Either NULL and usefs dNdScv default covariates, or one of these:  "bladder_pca"  "breast_pca"
#' "cesc_pca" "colon_pca" "esca_pca" "gbm_pca" "hnsc_pca" "kidney_pca" "lihc_pca" "lung_pca" "ov_pca" "pancreas_pca" "prostate_pca" "rectum_pca" "skin_pca"  "stomach_pca"  "thca_pca" "ucec_pca"
#'
#' @return
#' @export
#'
#'
#'

# Runs dNdScv and calculates gene-level mutation rates
gene_level_mutation_rates <- function(cesa, covariate_file = NULL){
  

  if(is.null(covariate_file)){
    message("Loading dNdScv default covariates for hg19...")
    data("covariates_hg19",package = "dndscv")
    genes_in_pca <- rownames(covs)
  }else{
    message("Loading tissue covariates file...")
    this_cov_pca <- get(data(list=covariate_file, package="cancereffectsizeR"))
    genes_in_pca <- rownames(this_cov_pca$rotation)
  }

  path_to_library <- dir(.libPaths(),full.names=T)[grep(dir(.libPaths(),full.names=T),pattern="cancereffectsizeR")][1] # find the path to this package
  data("RefCDS_TP53splice",package = "cancereffectsizeR")

  # store dndscv data in a list, split by tumor progression subsets

  dndscv_out_list <- vector(mode = "list",length = length(cesa@progressions@order))
  names(dndscv_out_list) <- names(cesa@progressions@order)

  # select MAF records to use for gene mutation rate calculation (via dNdScv)
  # for now, records from TGS samples are kept out; in the future, we could check out dNdScv's panel sequencing features
  exome_samples = cesa@coverage$samples_by_coverage[["exome"]]
  dndscv_maf = cesa@maf[cesa@maf$Unique_Patient_Identifier %in% exome_samples,]
  
  # dndscv output for each subset
  message("Running dNdScv...")
  for(this_subset in 1:length(dndscv_out_list)){
    current_subset_tumors = get_progression_tumors(cesa@progressions, this_subset)
    dndscv_out_list[[this_subset]] <- dndscv::dndscv(
      mutations = dndscv_maf[dndscv_maf$Unique_Patient_Identifier %in% current_subset_tumors,],
      gene_list = genes_in_pca,
      cv = if(is.null(covariate_file)){ "hg19"}else{ this_cov_pca$rotation},
      refdb = paste(path_to_library,"/data/RefCDS_TP53splice.RData",sep=""))
  }



  RefCDS_our_genes <- RefCDS[which(sapply(RefCDS, function(x) x$gene_name) %in% dndscv_out_list[[1]]$genemuts$gene_name)]
  mutrates_list <- vector(mode = "list",length = length(cesa@progressions@order))
  names(mutrates_list) <- names(cesa@progressions@order)


  message("Calculating gene-level mutation rates...")
  for(this_subset in 1:length(mutrates_list)){
    dndscv_output = dndscv_out_list[[this_subset]]
    RefCDS_gene_names <- sapply(RefCDS, function(x) x$gene_name)
    names(RefCDS) <- RefCDS_gene_names
    this_RefCDS <- RefCDS[as.character(dndscv_output$genemuts$gene_name)]
    number_of_syn_mutations <- sapply(this_RefCDS, function(x) colSums(x$L)[1])


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
        number_of_syn_mutations /
        number_of_tumors_in_this_subset

    } else{

      mutrates_vec <- rep(NA,length(dndscv_output$genemuts$exp_syn_cv))

      syn_sites <- number_of_syn_mutations
      for(i in 1:length(mutrates_vec)){
        mutrates_vec[i] <-  max(dndscv_output$genemuts$exp_syn_cv[i],  ((dndscv_output$genemuts$n_syn[i] +
                                                                           dndscv_output$nbreg$theta -
                                                                           1) /
                                                                          (1 +
                                                                             (dndscv_output$nbreg$theta /
                                                                                dndscv_output$genemuts$exp_syn_cv[i])
                                                                          ))) /
          number_of_syn_mutations[i] /
          number_of_tumors_in_this_subset
      }
    }
    names(mutrates_vec) <- dndscv_output$genemuts$gene_name
    mutrates_list[[this_subset]] = mutrates_vec
  }
  cesa@mutrates_list = mutrates_list
  cesa@dndscv_out_list = dndscv_out_list
  cesa@refcds_data = RefCDS_our_genes 
  return(cesa)
}
