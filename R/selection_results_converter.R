#' Selection results converter
#'
#' @import dplyr
#' @import magrittr
#'
#' @param results_input List of the selection results from get_gene_results within effect_size_SNV
#' @param subset_greater_than_freq subsetting the data for only substitutions in more than this number of tumors
#'
#' @return
#' @export
#'
#' @examples
#'
selection_results_converter <- function(results_input, subset_greater_than_freq=1){

  gene_adder <- function(x){
    holder <- x$selection_results
    holder$gene <- x$gene_name
    return(holder)
  }

  selection_results <- results_input$selection_output

  subsets <- results_input$MAF %>%
    dplyr::group_by(subset_col) %>%
    summarize(unique_tumors = n_distinct(Unique_patient_identifier))

  levels_in_selection_analysis <- levels(results_input$MAF[,"subset_col"])

  # gene_names <- sapply(selection_results, function(x) x$gene_name)

  # selection_data <- sapply(selection_results, function(x) x$selection_results)

  selection_data <- selection_results %>%
    map(gene_adder) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(variant =
             case_when(sum(charToRaw(variant) == charToRaw(" "))==0 ~ paste(gene, variant),
                       sum(charToRaw(variant) == charToRaw(" "))>0 ~ variant)) %>%
    ungroup() %>%
    unnest() %>%
    add_column(subset=rep(levels_in_selection_analysis,nrow(.)/length(levels_in_selection_analysis)))

  selection_data$population_proportion <- NA

  for(subset_index in 1:nrow(subsets)){
    selection_data[which(selection_data[,"subset"] == as.character(as.matrix(subsets[subset_index,1]))),"population_proportion"] <- selection_data[which(selection_data[,"subset"] == as.character(as.matrix(subsets[subset_index,1]))),"variant_freq"] / as.numeric(subsets[subset_index,"unique_tumors"])
  }

  selection_data <- selection_data %>%
    mutate(population_percent =
             paste(round(population_proportion*100,2), "%", sep=""))


  selection_data_df <- as.data.frame(selection_data,stringsAsFactors =F)
  selection_data_df$variant_freq <- as.numeric(selection_data_df$variant_freq)
  selection_data_df <- selection_data_df[selection_data_df$variant_freq>subset_greater_than_freq,]
  selection_data_df$selection_intensity <- as.numeric(selection_data_df$selection_intensity)

  dndscv_results <- vector(mode = "list",length = length(levels_in_selection_analysis))
  names(dndscv_results) <- names(results_input$dndscvout)

  for(subset_index in 1:length(dndscv_results)){
    dndscv_results[[subset_index]] <- results_input$dndscvout[[subset_index]]$sel_cv
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
