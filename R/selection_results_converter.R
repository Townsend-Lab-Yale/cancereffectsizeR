#' Selection results converter
#'
#' @param results_input List of the selection results from get_gene_results
#' @param subset_greater_than_freq subsetting the data for only substitutions in more than this number of tumors
#'
#' @return
#' @export
#'
#' @examples
#'
selection_results_converter <- function(results_input, subset_greater_than_freq=1){

  selection_results <- results_input$selection_output

  gene_names <- sapply(selection_results, function(x) x$gene_name)

  selection_data <- sapply(selection_results, function(x) x$selection_results)

  # adding gene names to the variant names
  for(i in 1:length(selection_data)){
    selection_data[[i]][which(sapply(strsplit(selection_data[[i]][,1],split = " "), function(x) length(x))==1),1] <- paste(gene_names[i],selection_data[[i]][which(sapply(strsplit(selection_data[[i]][,1],split = " "), function(x) length(x))==1),1],sep=" ")
    rownames(selection_data[[i]]) <- as.character(selection_data[[i]][,1])
    selection_data[[i]] <- cbind(selection_data[[i]], gene_names[i])

  }
  selection_data <- do.call(rbind,selection_data)

  selection_data_df <- as.data.frame(selection_data,stringsAsFactors =F)
  selection_data_df$variant_freq <- as.numeric(selection_data_df$variant_freq)
  selection_data_df <- selection_data_df[selection_data_df$variant_freq>subset_greater_than_freq,]
  selection_data_df$selection_intensity <- as.numeric(selection_data_df$selection_intensity)

  dndscv_results <- selection_output$dndscvout$sel_cv
  rownames(dndscv_results) <- dndscv_results$gene_name

  selection_data_df$dndscv_q  <- dndscv_results[selection_data_df$V6,"qallsubs_cv"]

  results_output <- selection_data_df[order(selection_data_df$selection_intensity,decreasing = T),]


  colnames(results_output)[6] <- "Gene"

  return(results_output)
}
