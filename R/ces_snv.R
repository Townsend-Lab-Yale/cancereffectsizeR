#' Calculate SNV selection intensity
#' @param cesa CESAnalysis object
#' @param gene which genes to calculate effect sizes within; defaults to all genes with recurrent mutations in data set
#' @param cores number of cores to use
#' @param include_genes_without_recurrent_mutations default false; will increase runtime and won't find anything interesting
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export

ces_snv <- function(cesa = NULL,
                            genes = "all",
                            cores = 1,
                            include_genes_without_recurrent_mutations = F,
                            find_CI=T) 
{

  if (! "character" %in% class(genes)) {
    stop("Expected argument \"genes\" to take a character vector.")
  }

  # using the "SNV" genes
  snv.maf = cesa@annotated.snv.maf
  genes_in_dataset = unique(snv.maf$Gene_name)
  if(length(genes_in_dataset) == 0) {
    stop("The SNV mutation data set is empty!")
  }


  # Possibly make this only load when needed (but it's usually needed)
  data("RefCDS_TP53splice", package = "cancereffectsizeR", envir = environment())
  names(RefCDS) = sapply(RefCDS, function(x) x$gene_name)


  if(genes[1] =="all") {
    if (include_genes_without_recurrent_mutations) {
      genes_to_analyze <- genes_in_dataset
    } else {
      tmp = table(snv.maf$unique_variant_ID)
      recurrent_variants = names(tmp[tmp > 1])
      has_recurrent = snv.maf$unique_variant_ID %in% recurrent_variants
      genes_to_analyze = unique(snv.maf[has_recurrent, "Gene_name"])
    }
    message(paste(length(genes_in_dataset) - length(genes_to_analyze), "genes in the data set have no recurrent SNV mutations."))
    message(paste("Calculating selection intensity for recurrent SNV mutations across", length(genes_to_analyze), "genes."))
  } else{
    genes = unique(genes)
    genes_to_analyze <- genes[genes %in% genes_in_dataset]
    missing_genes = genes[! genes %in% genes_in_dataset]
    invalid_genes = missing_genes[! missing_genes %in% names(RefCDS)]
    num_invalid = length(invalid_genes)
    if(num_invalid > 0) {
      additional_msg = ""
      if(num_invalid > 50) {
        invalid_genes = invalid_genes[1:40]
        additional_msg = paste0(" (and ", num_invalid - 40, " more)")
      }
      list_of_invalid = paste(invalid_genes, collapse = ", ")
      stop(paste0("Note: The following requested genes have no reference data (for genome build ", cesa@genome_build, "):\n\t",
                  list_of_invalid, additional_msg, "\n"))
    }
    if (length(genes_to_analyze) == 0) {
      stop("None of the requested genes have mutations in the SNV data set.")
    }

    num_missing = length(missing_genes)
    if(num_missing > 0) {
      additional_msg = ""
      if(num_missing > 50) {
        missing_genes = missing_genes[1:40]
        additional_msg = paste0(" (and ", num_missing - 40, " more)")
      }
      list_of_missing = paste(missing_genes, collapse = ", ")

      message(paste0("The following requested genes have no mutations in the SNV data set, so they won't be analyzed:\n\t",
        list_of_missing, additional_msg))
    }
  }




  selection_results <- vector("list",length = length(genes_to_analyze))
  names(selection_results) <- genes_to_analyze
  
  

  all_tumors = unique(snv.maf$Unique_Patient_Identifier)

  selection_results <- pbapply::pblapply(genes_to_analyze, get_gene_results, cesa = cesa,
                                            all_tumors = all_tumors,find_CI=find_CI, RefCDS = RefCDS, cl = cores)

  cesa@selection_results = selection_results
  return(cesa)
}




#' Single-stage SNV effect size analysis (gets called by ces_snv)
get_gene_results <- function(gene_to_analyze, cesa, all_tumors, find_CI, RefCDS) {
  snv.maf = cesa@annotated.snv.maf
  progressions = cesa@progressions
  current_gene_maf = snv.maf[snv.maf$Gene_name == gene_to_analyze,]
  these_mutation_rates <-
    cancereffectsizeR::mutation_rate_calc(
      this_MAF = current_gene_maf,
      gene = gene_to_analyze,
      gene_mut_rate = cesa@mutrates_list,
      trinuc_proportion_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix,
      gene_refcds = RefCDS[[gene_to_analyze]],
      all_tumors = all_tumors,
      progressions = progressions)


  # begin populating these_selection_results, which are selection results for all variants in given gene
  ## rows = all variants found in mutatioin rate matrix columns
  these_selection_results <- dplyr::tibble(variant = colnames(these_mutation_rates$mutation_rate_matrix),
                                           selection_intensity = vector(mode = "list",
                                                                        length = ncol(these_mutation_rates$mutation_rate_matrix)),
                                           unsure_gene_name=NA,
                                           unique_variant_ID=NA,
                                           loglikelihood=NA)

  if(length(progressions@order) == 1 & find_CI) {
    # add columns corresponding to .999% and .95% confidence intervals
    these_selection_results <- cbind(these_selection_results,
                                     dplyr::tibble(ci_low_999=rep(NA,ncol(these_mutation_rates$mutation_rate_matrix)),ci_high_999=NA,ci_low_95=NA,ci_high_95=NA))

  }

  num_selection_results = nrow(these_selection_results)
  for(j in 1:num_selection_results){
    variant = colnames(these_mutation_rates$mutation_rate_matrix)[j]
    current_locus = GenomicRanges::makeGRangesFromDataFrame(snv.maf[snv.maf$unique_variant_ID_AA == variant, ][1,], seqnames.field = "Chromosome",
                                             start.field = "Start_Position", end.field = "Start_Position")
    eligible_tumors = cancereffectsizeR:::get_tumors_with_coverage(coverage = cesa@coverage, locus = current_locus)
    
    # toss out tumors that do not have any SNVs to analyze
    eligible_tumors = eligible_tumors[eligible_tumors %in% all_tumors]
    
    optimization_output <- cancereffectsizeR::optimize_gamma(
      MAF_input = current_gene_maf,
      eligible_tumors = eligible_tumors,
      progressions = progressions,
      gene=gene_to_analyze,
      variant=variant,
      specific_mut_rates=these_mutation_rates$mutation_rate_matrix)

    these_selection_results[j,c("selection_intensity")][[1]] <-
      list(optimization_output$par)


    names(these_selection_results[j,c("selection_intensity")][[1]][[1]]) <- names(progressions@order)

    these_selection_results[j,"loglikelihood"] <- optimization_output$value


    if(length(progressions@order) == 1 & find_CI){
      # find CI function
      CI_results <- cancereffectsizeR::CI_finder(gamma_max = optimization_output$par,
                                                 MAF_input= current_gene_maf,
                                                 eligible_tumors = eligible_tumors,
                                                 progressions = progressions,
                                                 gene=gene_to_analyze,
                                                 variant=colnames(these_mutation_rates$mutation_rate_matrix)[j],
                                                 specific_mut_rates=these_mutation_rates$mutation_rate_matrix)


      these_selection_results[j,"ci_low_999"] <- CI_results$lower_CI
      these_selection_results[j,"ci_high_999"] <- CI_results$upper_CI

      CI_results <- cancereffectsizeR::CI_finder(gamma_max = optimization_output$par,
                                                 MAF_input= current_gene_maf,
                                                 eligible_tumors = eligible_tumors,
                                                 progressions = progressions,
                                                 gene=gene_to_analyze,
                                                 variant=colnames(these_mutation_rates$mutation_rate_matrix)[j],
                                                 specific_mut_rates=these_mutation_rates$mutation_rate_matrix,
                                                 log_units_down = 1.92 # 95% confidence interval
                                                 )


      these_selection_results[j,"ci_low_95"] <- CI_results$lower_CI
      these_selection_results[j,"ci_high_95"] <- CI_results$upper_CI


    }


  }

  these_selection_results[,"unsure_gene_name"] <- these_mutation_rates$unsure_genes_vec
  these_selection_results[,"unique_variant_ID"] <- these_mutation_rates$unique_variant_ID

  return(list(gene_name=gene_to_analyze, selection_results=these_selection_results))
}





