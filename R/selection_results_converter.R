#' Selection results converter
#'
#' @importFrom magrittr "%>%"
#'
#' @param cesa CESAnalysis object with selection intensity results
#' @param min_recurrence   minimum number of samples in which a mutation must appear for the mutation to be included in results
#'
#' @return
#' @export
#'
#' @examples
#'
selection_results_converter <- function(cesa, min_recurrence = 2){

  gene_adder <- function(x){
    holder <- x$selection_results
    holder$gene <- x$gene_name
    return(holder)
  }

  selection_results = cesa@selection_results
  progression_names <- names(cesa@progressions@order)



  selection_data <- selection_results %>%
    lapply(., gene_adder) %>%
    dplyr::bind_rows() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(variant =
             dplyr::case_when(sum(charToRaw(variant) == charToRaw(" "))==0 ~ paste(gene, variant),
                       sum(charToRaw(variant) == charToRaw(" "))>0 ~ variant)) %>%
    dplyr::ungroup() %>%
    tidyr::unnest()

    selection_data$subset=rep(progression_names,nrow(selection_data)/length(progression_names))



  selection_data_df <- as.data.frame(selection_data,stringsAsFactors =F)
  
  
  
  # returns vector with number of tumors with coverage per progression stage
  # get_num_tumors_with_coverage(variant_aa) {
  #     pass
  # }
  
  maf = cesa@annotated.snv.maf
  
  tumors_by_stage = lapply(progression_names, function(x) get_progression_tumors(cesa@progressions, x))
  names(tumors_by_stage) = progression_names
  num_tumors_by_stage = lapply(tumors_by_stage, length)
  names(num_tumors_by_stage) = progression_names
  table_by_stage = lapply(tumors_by_stage, function(x) table(maf$unique_variant_ID[maf$Unique_Patient_Identifier %in% x] ))
  
  variant_freq = numeric(nrow(selection_data_df))
  population_proportion = variant_freq
  
  for (i in 1:nrow(selection_data_df)) {
    stage = selection_data_df[i, "subset"]
    variant_freq[i] = as.numeric(table_by_stage[[stage]][selection_data_df$unique_variant_ID[i]])
    population_proportion[i] = variant_freq[i]/num_tumors_by_stage[[stage]]
  }

  selection_data_df$variant_freq = variant_freq
  selection_data_df$population_proportion = population_proportion
  
  selection_data_df <- selection_data_df[! is.na(selection_data_df$variant_freq) & selection_data_df$variant_freq>=min_recurrence,]
  selection_data_df$selection_intensity <- as.numeric(selection_data_df$selection_intensity)

  dndscv_results <- vector(mode = "list",length = length(progression_names))
  names(dndscv_results) <- names(cesa@dndscv_out_list)

  for(subset_index in 1:length(dndscv_results)){
    dndscv_results[[subset_index]] <- cesa@dndscv_out_list[[subset_index]]$sel_cv
    rownames(dndscv_results[[subset_index]]) <- dndscv_results[[subset_index]]$gene_name

    if(length(which(selection_data_df$subset == progression_names[subset_index]))>0){
    selection_data_df[which(selection_data_df$subset ==
                              progression_names[subset_index]),
                      "dndscv_q"] <-
      dndscv_results[[subset_index]][
        selection_data_df[which(selection_data_df$subset ==
                                  progression_names[subset_index]),"gene"],
        "qallsubs_cv"]
    }
  }

  results_output <- selection_data_df[order(selection_data_df$selection_intensity,decreasing = T),]


  return(results_output)
}
