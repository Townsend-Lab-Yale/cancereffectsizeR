#' Gene-level mutation rate function
#'
#' This function calculates the maximum likelihood estimate for the expected
#' number of synonymous mutations in a gene under the dNdScv model, as described within
#' Martincorena, I., Raine, K. M., Gerstung, M., Dawson, K. J., Haase, K., Van Loo, P., … Campbell, P. J.
#' (2017). Universal Patterns of Selection in Cancer and Somatic Tissues. Cell. https://doi.org/10.1016/j.cell.2017.09.042
#' and the associated R package dndscv .
#'
#'
#' @param dndscv_output output from dndscv
#' @param RefCDS_object object in the style of 'RefCDS' within the dndscv package
#'
#' @return
#' @export
#'
#' @examples
#'
#'
gene_level_mutation_rates <- function(dndscv_output, RefCDS_object){


  RefCDS_gene_names <- sapply(RefCDS_object, function(x) x$gene_name)
  names(RefCDS_object) <- RefCDS_gene_names
  this_RefCDS <- RefCDS_object[as.character(dndscv_output$genemuts$gene_name)]
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

  }else{

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




  return(mutrates_vec)


}
