#' Selection results converter
#'
#' @importFrom magrittr "%>%"
#'
#' @param cesa CESAnalysis object with selection intensity results
#' @param subset_greater_than_freq subsetting the data for only substitutions in more than this number of tumors
#'
#' @return
#' @export
#'
#' @examples
#'
selection_results_converter <- function(cesa, subset_greater_than_freq=1){

  gene_adder <- function(x){
    holder <- x$selection_results
    holder$gene <- x$gene_name
    return(holder)
  }

  selection_results = cesa@selection_results
  levels_in_selection_analysis <- names(cesa@progressions@order)

  subset_col = factor(levels_in_selection_analysis, levels = levels_in_selection_analysis, ordered = T)
  unique_tumors = sapply(cesa@progressions@order, function(x) length(get_progression_tumors(cesa@progressions, x)))

  subsets = dplyr::tibble(subset_col = subset_col, unique_tumors = unique_tumors)

  selection_data <- selection_results %>%
    lapply(., gene_adder) %>%
    dplyr::bind_rows() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(variant =
             dplyr::case_when(sum(charToRaw(variant) == charToRaw(" "))==0 ~ paste(gene, variant),
                       sum(charToRaw(variant) == charToRaw(" "))>0 ~ variant)) %>%
    dplyr::ungroup() %>%
    tidyr::unnest()

    selection_data$subset=rep(levels_in_selection_analysis,nrow(selection_data)/length(levels_in_selection_analysis))

  selection_data$population_proportion <- NA

  for(subset_index in 1:nrow(subsets)){
    selection_data[which(selection_data[,"subset"] == as.character(as.matrix(subsets[subset_index,1]))),"population_proportion"] <- selection_data[which(selection_data[,"subset"] == as.character(as.matrix(subsets[subset_index,1]))),"variant_freq"] / as.numeric(subsets[subset_index,"unique_tumors"])
  }

  selection_data <- selection_data %>%
    dplyr::mutate(population_percent =
             paste(round(population_proportion*100,2), "%", sep=""))


  selection_data_df <- as.data.frame(selection_data,stringsAsFactors =F)
  selection_data_df$variant_freq <- as.numeric(selection_data_df$variant_freq)
  selection_data_df <- selection_data_df[selection_data_df$variant_freq>subset_greater_than_freq,]
  selection_data_df$selection_intensity <- as.numeric(selection_data_df$selection_intensity)

  dndscv_results <- vector(mode = "list",length = length(levels_in_selection_analysis))
  names(dndscv_results) <- names(cesa@dndscv_out_list)

  for(subset_index in 1:length(dndscv_results)){
    dndscv_results[[subset_index]] <- cesa@dndscv_out_list[[subset_index]]$sel_cv
    rownames(dndscv_results[[subset_index]]) <- dndscv_results[[subset_index]]$gene_name

    if(length(which(selection_data_df$subset == levels_in_selection_analysis[subset_index]))>0){
    selection_data_df[which(selection_data_df$subset ==
                              levels_in_selection_analysis[subset_index]),
                      "dndscv_q"] <-
      dndscv_results[[subset_index]][
        selection_data_df[which(selection_data_df$subset ==
                                  levels_in_selection_analysis[subset_index]),"gene"],
        "qallsubs_cv"]
    }
  }

  # selection_data_df$dndscv_q  <- dndscv_results[selection_data_df[,"gene"],"qallsubs_cv"]

  results_output <- selection_data_df[order(selection_data_df$selection_intensity,decreasing = T),]


  return(results_output)
}
