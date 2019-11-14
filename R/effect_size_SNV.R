#' Calculate SNV selection intensity
#' @param cesa CESAnalysis object
#' @param gene which genes to calculate effect sizes within; defaults to all genes with recurrent mutations in data set
#' @param include_genes_without_recurrent_mutations default false; will greatly slow runtime and won't find anything interesting
#' @param analysis choose SNV, gene-level-epistasis, or top-epistasis
#' @param cores number of cores to use
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export

effect_size_SNV <- function(
  cesa = NULL,
  genes = "all",
  analysis = c("SNV", "gene-level-epistasis", "top-epistasis"),
  epistasis_top_prev_number = NULL,
  cores = 1,
  include_genes_without_recurrent_mutations = F,
  full_gene_epistasis_lower_optim = 1e-3,
  full_gene_epistasis_upper_optim=1e9,
  full_gene_epistasis_fnscale=-1e-16,
  q_threshold_for_gene_level=0.1,
  find_CI=T) {

  # validate analysis choice
  analysis = match.arg(analysis)

  if (analysis != "top-epistasis" && ! is.null(epistasis_top_prev_number)) {
    stop("Error: epistasis_top_prev_number is only used for \"top-epistasis\" analysis")
  } else if(! is.null(epistasis_top_prev_number)) {
    stop("Error: epistasis_top_prev_number must be defined for top-epistasis analysis")
  }

  # can't combine epistasis analysis with multi-stage yet
  if (analysis != "SNV" && length(cesa@progressions@order) > 1) {
    stop("Epistasis analysis is not compatible yet with multi-stage analyses. You'll have to re-run from the beginning.")
  }

  if (! "character" %in% class(genes)) {
    stop("Expected -genes to take a character vector.")
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




  data("gene_trinuc_comp", package = "cancereffectsizeR", envir = environment())
  names(gene_trinuc_comp) <- sapply(gene_trinuc_comp, function(x) x$gene_name) # formally used names of RefCDS, but they're the same


  data("AA_mutation_list", package="cancereffectsizeR")
  data("AA_translations", package="cancereffectsizeR")

  selection_results <- vector("list",length = length(genes_to_analyze))
  names(selection_results) <- genes_to_analyze
  
  
  # divide MAF data by gene
  MAF = snv.maf # temp
  mafs = new.env(hash=TRUE)
  for (gene in genes_to_analyze) {
    gene_maf = MAF[MAF[,"Gene_name"] == gene,]
    mafs[[gene]][["maf"]] = gene_maf
  }

  all_tumors = unique(snv.maf$Unique_Patient_Identifier)

  if (analysis == "SNV") {
    selection_results <- pbapply::pblapply(genes_to_analyze, get_gene_results, cesa = cesa,
                                            gene_mafs = mafs, gene_trinuc_comp = gene_trinuc_comp,
                                            all_tumors = all_tumors,find_CI=find_CI, RefCDS = RefCDS, cl = cores)
  } else if(analysis == "gene-level-epistasis") {


    cesa_subset <- cesa@annotated.snv.maf[cesa@annotated.snv.maf$Gene_name %in% genes_to_analyze,]

    cesa_subset$identifier <- cesa_subset$unique_variant_ID_AA
    cesa_subset$identifier[
      which(sapply(strsplit(cesa_subset$identifier,split = " "),
                   function(x) length(x))==1)
      ] <- paste(cesa_subset$Gene_name[
        which(sapply(strsplit(cesa_subset$identifier,split = " "),
                     function(x) length(x))==1)],
        cesa_subset$identifier[
          which(sapply(strsplit(cesa_subset$identifier,split = " "),
                       function(x) length(x))==1)
          ],sep=" ")

    cesa_subset_table <- table(cesa_subset$identifier)
    cesa_subset$recurrent_val <- cesa_subset_table[cesa_subset$identifier]

    recurrently_subed_genes <- unique(cesa_subset[cesa_subset$recurrent_val>1,"Gene_name"])

    if(length(recurrently_subed_genes) < 2){
     stop("Less than 2 of the 'genes' you specified have recurrent variants.
          Choose 2 genes or more with recurrent variants to measure epistasis
          among those recurrently substituted variants. ")
    }

    genes_to_analyze <- recurrently_subed_genes

    selection_epistasis_results <- t(utils::combn(genes_to_analyze,2))
    selection_epistasis_results <- data.frame(t(selection_epistasis_results),stringsAsFactors=F)
    rownames(selection_epistasis_results) <- c("Variant_1","Variant_2")
    selection_epistasis_results_list <- as.list(selection_epistasis_results)


    selection_results = pbapply::pblapply(X = selection_epistasis_results_list,
                                           FUN = epistasis_gene_level,
                                           MAF=cesa@annotated.snv.maf,
                                           trinuc_proportion_matrix=cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix,
                                           cesa=cesa,
                                           gene_trinuc_comp = gene_trinuc_comp,
                                           all_tumors = all_tumors,
                                           RefCDS = RefCDS,
                                           cl = cores)
  } else if(analysis == "top-epistasis") {
    selection_results = top_epistasis()
  }
  cesa@selection_results = selection_results
  return(cesa)
}




#' Single-stage SNV effect size analysis (gets called by effect_size_SNV)
get_gene_results <- function(gene_to_analyze, cesa, gene_mafs, gene_trinuc_comp, all_tumors, find_CI, RefCDS) {
  mutrates_list = cesa@mutrates_list
  trinuc_proportion_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  progressions = cesa@progressions
  these_mutation_rates <-
    cancereffectsizeR::mutation_rate_calc(
      this_MAF = gene_mafs[[gene_to_analyze]][["maf"]],
      gene = gene_to_analyze,
      gene_mut_rate = mutrates_list,
      trinuc_proportion_matrix = trinuc_proportion_matrix,
      gene_trinuc_comp = gene_trinuc_comp,
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
    MAF_input = gene_mafs[[gene_to_analyze]][["maf"]]
    current_locus = GenomicRanges::makeGRangesFromDataFrame(MAF_input[MAF_input$unique_variant_ID_AA == variant, ][1,], seqnames.field = "Chromosome",
                                             start.field = "Start_Position", end.field = "Start_Position")
    eligible_tumors = cancereffectsizeR:::get_tumors_with_coverage(coverage = cesa@coverage, locus = current_locus)
    
    # toss out tumors that do not have any SNVs to analyze
    eligible_tumors = eligible_tumors[eligible_tumors %in% all_tumors]
    
    optimization_output <- cancereffectsizeR::optimize_gamma(
      MAF_input = MAF_input,
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
                                                 MAF_input= gene_mafs[[gene_to_analyze]][["maf"]],
                                                 eligible_tumors = eligible_tumors,
                                                 progressions = progressions,
                                                 gene=gene_to_analyze,
                                                 variant=colnames(these_mutation_rates$mutation_rate_matrix)[j],
                                                 specific_mut_rates=these_mutation_rates$mutation_rate_matrix)


      these_selection_results[j,"ci_low_999"] <- CI_results$lower_CI
      these_selection_results[j,"ci_high_999"] <- CI_results$upper_CI

      CI_results <- cancereffectsizeR::CI_finder(gamma_max = optimization_output$par,
                                                 MAF_input= gene_mafs[[gene_to_analyze]][["maf"]],
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






## Gets called by effect_size_SNV for gene-level epistasis
epistasis_gene_level = function(genes_to_analyze,
                                MAF,
                                trinuc_proportion_matrix,
                                cesa,
                                gene_trinuc_comp,
                                all_tumors,
                                RefCDS) {


  mutrates_list = cesa@mutrates_list
  trinuc_proportion_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  progressions = cesa@progressions


  get_gene_results_epistasis_bygene <- function(variant_combo_list) {

    # print(variant_combo_list)

    variant1 <- variant_combo_list[1]
    variant2 <- variant_combo_list[2]

    MAF_input1=MAF[MAF[,"Gene_name"] == variant1 &
                     MAF[,"Reference_Allele"] %in% c("A","T","G","C") &
                     MAF[,"Tumor_Allele"] %in% c("A","T","G","C"),]

    MAF_input2=MAF[MAF[,"Gene_name"] == variant2 &
                     MAF[,"Reference_Allele"] %in% c("A","T","G","C") &
                     MAF[,"Tumor_Allele"] %in% c("A","T","G","C"),]
    
    # when including target gene sequencing data, need to throw out any samples that don't have coverage at ALL variant sites in these genes
    # this could be a problem if some of the "exome" samples are actually whole-genome if they haven't been trimmed strictly enough
    
    # this is clunky but okay for now and will change with next refactor
    eligible_tumors = all_tumors # samples without coverage at all sites will get intersected out
    for (maf in list(MAF_input1, MAF_input2)) {
      for (i in 1:nrow(maf)) {
        current_locus = GenomicRanges::makeGRangesFromDataFrame(maf[i,], seqnames.field = "Chromosome",
                                                                start.field = "Start_Position", end.field = "Start_Position")
        tumors_covering_locus = cancereffectsizeR:::get_tumors_with_coverage(coverage = cesa@coverage, locus = current_locus)
        eligible_tumors = intersect(eligible_tumors, tumors_covering_locus) 
      }
    }
    
    # restrict MAFs to eligible tumors
    MAF_input1 = MAF_input1[MAF_input1$Unique_Patient_Identifier %in% eligible_tumors,]
    MAF_input2 = MAF_input2[MAF_input2$Unique_Patient_Identifier %in% eligible_tumors,]
    
    variant_freq_1 <- table(MAF_input1$unique_variant_ID_AA)
    variant_freq_2 <- table(MAF_input2$unique_variant_ID_AA)

    # only run the selection algorithm if there are 2 or more tumors with
    # recurrent variants of each gene present.
    these_mutation_rates1 <-
      cancereffectsizeR::mutation_rate_calc(
        this_MAF = MAF_input1,
        gene = variant1,
        gene_mut_rate = mutrates_list,
        trinuc_proportion_matrix = trinuc_proportion_matrix,
        gene_trinuc_comp = gene_trinuc_comp,
        gene_refcds = RefCDS[[variant1]],
        all_tumors = eligible_tumors,
        progressions = progressions)

    these_mutation_rates2 <-
      cancereffectsizeR::mutation_rate_calc(
        this_MAF = MAF_input2,
        gene = variant2,
        gene_mut_rate = mutrates_list,
        trinuc_proportion_matrix = trinuc_proportion_matrix,
        gene_trinuc_comp = gene_trinuc_comp,
        gene_refcds = RefCDS[[variant2]],
        all_tumors = eligible_tumors,
        progressions = progressions)




    # since we are looking at selection at the gene level,
    # we only consider selection at sites that are recurrently
    # substituted.

    these_mutation_rates1$mutation_rate_matrix <- as.matrix(these_mutation_rates1$mutation_rate_matrix[,colnames(these_mutation_rates1$mutation_rate_matrix) %in% names(variant_freq_1[variant_freq_1>1])])
    these_mutation_rates2$mutation_rate_matrix <- as.matrix(these_mutation_rates2$mutation_rate_matrix[,colnames(these_mutation_rates2$mutation_rate_matrix) %in% names(variant_freq_2[variant_freq_2>1])])

    these_selection_results <- vector(mode = "list",length = 10)
    names(these_selection_results) <- c("Variant_1",
                                        "Variant_2",
                                        "Gamma_1",
                                        "Gamma_2",
                                        "Gamma_1_2background",
                                        "Gamma_2_1background",
                                        "tumors_with_ONLY_variant1_substituted",
                                        "tumors_with_ONLY_variant2_substituted",
                                        "tumors_with_both_substituted",
                                        "tumors_with_neither_substituted")

    these_selection_results[1:2] <- c(variant1,variant2)

    these_selection_results[3:6] <-
      cancereffectsizeR::optimize_gamma_epistasis_full_gene(
        MAF_input1=MAF_input1,
        MAF_input2=MAF_input2,
        all_tumors=eligible_tumors,
        gene1=variant1,
        gene2=variant2,
        specific_mut_rates1=these_mutation_rates1$mutation_rate_matrix,
        specific_mut_rates2=these_mutation_rates2$mutation_rate_matrix,
        variant_freq_1 = variant_freq_1,
        variant_freq_2 = variant_freq_2,
        full_gene_epistasis_lower_optim=full_gene_epistasis_lower_optim,
        full_gene_epistasis_upper_optim=full_gene_epistasis_upper_optim,
        full_gene_epistasis_fnscale=full_gene_epistasis_fnscale)


    # only the tumors containing a recurrent variant factor into the selection analysis
    tumors_with_variant1_mutated <- MAF_input1[which(MAF_input1$unique_variant_ID_AA %in% names(which(variant_freq_1>1))),"Unique_Patient_Identifier"]
    tumors_with_variant2_mutated <- MAF_input2[which(MAF_input2$unique_variant_ID_AA %in% names(which(variant_freq_2>1))),"Unique_Patient_Identifier"]

    tumors_with_both_mutated <- base::intersect(tumors_with_variant1_mutated,tumors_with_variant2_mutated)

    tumors_with_ONLY_variant1_mutated <- tumors_with_variant1_mutated[which(!tumors_with_variant1_mutated %in% tumors_with_variant2_mutated)]

    tumors_with_ONLY_variant2_mutated <- tumors_with_variant2_mutated[which(!tumors_with_variant2_mutated %in% tumors_with_variant1_mutated)]

    # not all tumors with no mutations in the genes, just all that had coverage at all variant sites in both genes and no specific substitutions recurrent in other tumors
    tumors_with_neither_mutated <- setdiff(eligible_tumors, c(tumors_with_both_mutated,tumors_with_ONLY_variant1_mutated,tumors_with_ONLY_variant2_mutated))


    these_selection_results$tumors_with_ONLY_variant1_substituted <- tumors_with_ONLY_variant1_mutated
    these_selection_results$tumors_with_ONLY_variant2_substituted <- tumors_with_ONLY_variant2_mutated
    these_selection_results$tumors_with_both_substituted <- tumors_with_both_mutated
    these_selection_results$tumors_with_neither_substituted <- tumors_with_neither_mutated
    return(these_selection_results)
  }

  selection_results <- get_gene_results_epistasis_bygene(genes_to_analyze)
  return(selection_results)
}










top_epistasis = function(cesa, genes_to_analyze) {
  MAF$identifier <- MAF$unique_variant_ID_AA
  MAF$identifier[
    which(sapply(strsplit(MAF$identifier,split = " "),
                 function(x) length(x))==1)
    ] <- paste(MAF$Gene_name[
      which(sapply(strsplit(MAF$identifier,split = " "),
                   function(x) length(x))==1)],
      MAF$identifier[
        which(sapply(strsplit(MAF$identifier,split = " "),
                     function(x) length(x))==1)
        ],sep=" ")

  # Sort by top prevalences
  ID_prevalence <- table(MAF$identifier)[order(table(MAF$identifier),decreasing = T)]

  # Find the variants that match the numerical frequencies of the top variants
  # (in case there are a few sharing the nth value)
  ID_prevalence_top <- ID_prevalence[which(ID_prevalence %in% ID_prevalence[1:epistasis_top_prev_number])]
  selection_epistasis_results <- t(utils::combn(names(ID_prevalence_top),2))
  selection_epistasis_results <- data.frame(t(selection_epistasis_results),stringsAsFactors=F)
  rownames(selection_epistasis_results) <- c("Variant_1","Variant_2")
  selection_epistasis_results_list <- as.list(selection_epistasis_results)

  get_gene_results_epistasis <- function(variant_combo_list) {
    variant1 <- variant_combo_list[1]
    variant2 <- variant_combo_list[2]

    variant1_MAFindex <- which(MAF$identifier==variant1)[1]
    variant2_MAFindex <- which(MAF$identifier==variant2)[1]


    these_mutation_rates1 <-
      cancereffectsizeR::mutation_rate_calc(
        this_MAF = MAF[MAF[,"Gene_name"] == MAF[variant1_MAFindex,"Gene_name"] &
                         MAF[,ref_column] %in% c("A","T","G","C") &
                         MAF[,alt_column] %in% c("A","T","G","C"),],
        gene = MAF[variant1_MAFindex,"Gene_name"],
        gene_mut_rate = mutrates_list,
        trinuc_proportion_matrix = trinuc_proportion_matrix,
        gene_trinuc_comp = gene_trinuc_comp,
        gene_refcds = RefCDS[[gene]],
        tumor_subsets = tumors,subset_col=subset_col)

    these_mutation_rates2 <-
      cancereffectsizeR::mutation_rate_calc(
        this_MAF = MAF[MAF[,"Gene_name"] == MAF[variant2_MAFindex,"Gene_name"] &
                         MAF[,ref_column] %in% c("A","T","G","C") &
                         MAF[,alt_column] %in% c("A","T","G","C"),],
        gene = MAF[variant2_MAFindex,"Gene_name"],
        gene_mut_rate = mutrates_list,
        trinuc_proportion_matrix = trinuc_proportion_matrix,
        gene_trinuc_comp = gene_trinuc_comp,
        gene_refcds = RefCDS[[gene]],
        tumor_subsets = tumors,subset_col=subset_col)




    these_selection_results <- c(variant1,variant2,NA,NA,NA,NA)
    names(these_selection_results) <- c("Variant_1", "Variant_2", "Gamma_1", "Gamma_2",
                                        "Gamma_1_2background", "Gamma_2_1background")
    these_selection_results[3:6] <- cancereffectsizeR::optimize_gamma_epistasis(
      MAF_input1=MAF[MAF[,"Gene_name"] == MAF[variant1_MAFindex,"Gene_name"] &
                       MAF[,ref_column] %in% c("A","T","G","C") &
                       MAF[,alt_column] %in% c("A","T","G","C"),],
      MAF_input2=MAF[MAF[,"Gene_name"] == MAF[variant2_MAFindex,"Gene_name"] &
                       MAF[,ref_column] %in% c("A","T","G","C") &
                       MAF[,alt_column] %in% c("A","T","G","C"),],
      all_tumors=tumors,
      gene1=MAF[variant1_MAFindex,"Gene_name"],
      gene2=MAF[variant2_MAFindex,"Gene_name"],
      variant1= MAF[variant1_MAFindex,"unique_variant_ID_AA"],
      variant2= MAF[variant2_MAFindex,"unique_variant_ID_AA"],
      specific_mut_rates1=these_mutation_rates1$mutation_rate_matrix,
      specific_mut_rates2=these_mutation_rates2$mutation_rate_matrix)
    return(these_selection_results)
  }
  selection_results <- pbapply::pblapply(selection_epistasis_results_list, get_gene_results_epistasis, cl = cores)
  return(selection_results)
}

