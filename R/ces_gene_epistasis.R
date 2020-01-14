#' Calculate gene-level selection intensity under an assumption of pairwise gene epistasis
#' @param cesa CESAnalysis object
#' @param genes which genes to calculate effect sizes within; defaults to all genes with recurrent mutations in data set
#' @param cores number of cores to use
#' @param full_gene_epistasis_lower_optim  optimx parameter
#' @param full_gene_epistasis_upper_optim optimx parameter
#' @param full_gene_epistasis_fnscale optimx parameter
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export


ces_gene_epistasis = function(cesa = NULL, genes, cores = 1, full_gene_epistasis_lower_optim = 1e-3, 
						full_gene_epistasis_upper_optim = 1e9, full_gene_epistasis_fnscale = -1e-16)
{
  # can't combine epistasis analysis with multi-stage yet
	if (length(cesa@progressions@order) > 1) {
	  stop("Epistasis analysis is not compatible yet with multi-stage analyses. You'll have to create a new single-stage analysis.")
	}

  if (length(genes) < 2) {
    stop("Supply at least two genes to analyze.")
  }

  # temporarily hard-coded RefCDS data
  load(system.file("genomes/hg19/ces_hg19_tp53_splice_refcds_gr_genes.rda", package = "cancereffectsizeR"))


	genes = unique(genes)
	genes_in_dataset = unique(cesa@annotated.snv.maf$Gene_name)
	genes_to_analyze = genes[genes %in% genes_in_dataset]


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

	recurrently_subed_genes <- unique(cesa_subset$Gene_name[cesa_subset$recurrent_val > 1])

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
                                           all_tumors = unique(cesa@annotated.snv.maf$Unique_Patient_Identifier),
                                           RefCDS = RefCDS,
                                           cl = cores)
    cesa@selection_results = selection_results
    return(cesa)
}

epistasis_gene_level = function(genes_to_analyze,
                                MAF,
                                trinuc_proportion_matrix,
                                cesa,
                                all_tumors,
                                RefCDS) {
  mutrates_list = cesa@mutrates_list
  trinuc_proportion_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
  progressions = cesa@progressions


  get_gene_results_epistasis_bygene <- function(variant_combo_list) {

    # print(variant_combo_list)

    variant1 <- variant_combo_list[1]
    variant2 <- variant_combo_list[2]

    bases = c("A","T","G","C") 
    MAF_input1=data.frame(MAF[Gene_name == variant1 &
                     Reference_Allele %in% bases &
                     Tumor_Allele %in% bases])

    MAF_input2=data.frame(MAF[Gene_name == variant2 &
                     Reference_Allele %in% bases &
                     Tumor_Allele %in% bases])
    
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
        gene_refcds = RefCDS[[variant1]],
        all_tumors = eligible_tumors,
        progressions = progressions)

    these_mutation_rates2 <-
      cancereffectsizeR::mutation_rate_calc(
        this_MAF = MAF_input2,
        gene = variant2,
        gene_mut_rate = mutrates_list,
        trinuc_proportion_matrix = trinuc_proportion_matrix,
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






## currently broken (and will likely be replaced with other functions)
# top_epistasis = function(cesa, genes_to_analyze) {
#   MAF$identifier <- MAF$unique_variant_ID_AA
#   MAF$identifier[
#     which(sapply(strsplit(MAF$identifier,split = " "),
#                  function(x) length(x))==1)
#     ] <- paste(MAF$Gene_name[
#       which(sapply(strsplit(MAF$identifier,split = " "),
#                    function(x) length(x))==1)],
#       MAF$identifier[
#         which(sapply(strsplit(MAF$identifier,split = " "),
#                      function(x) length(x))==1)
#         ],sep=" ")

#   # Sort by top prevalences
#   ID_prevalence <- table(MAF$identifier)[order(table(MAF$identifier),decreasing = T)]

#   # Find the variants that match the numerical frequencies of the top variants
#   # (in case there are a few sharing the nth value)
#   ID_prevalence_top <- ID_prevalence[which(ID_prevalence %in% ID_prevalence[1:epistasis_top_prev_number])]
#   selection_epistasis_results <- t(utils::combn(names(ID_prevalence_top),2))
#   selection_epistasis_results <- data.frame(t(selection_epistasis_results),stringsAsFactors=F)
#   rownames(selection_epistasis_results) <- c("Variant_1","Variant_2")
#   selection_epistasis_results_list <- as.list(selection_epistasis_results)

#   get_gene_results_epistasis <- function(variant_combo_list) {
#     variant1 <- variant_combo_list[1]
#     variant2 <- variant_combo_list[2]

#     variant1_MAFindex <- which(MAF$identifier==variant1)[1]
#     variant2_MAFindex <- which(MAF$identifier==variant2)[1]


#     these_mutation_rates1 <-
#       cancereffectsizeR::mutation_rate_calc(
#         this_MAF = MAF[MAF[,"Gene_name"] == MAF[variant1_MAFindex,"Gene_name"] &
#                          MAF[,ref_column] %in% c("A","T","G","C") &
#                          MAF[,alt_column] %in% c("A","T","G","C"),],
#         gene = MAF[variant1_MAFindex,"Gene_name"],
#         gene_mut_rate = mutrates_list,
#         trinuc_proportion_matrix = trinuc_proportion_matrix,
#         gene_trinuc_comp = gene_trinuc_comp,
#         gene_refcds = RefCDS[[gene]],
#         tumor_subsets = tumors,subset_col=subset_col)

#     these_mutation_rates2 <-
#       cancereffectsizeR::mutation_rate_calc(
#         this_MAF = MAF[MAF[,"Gene_name"] == MAF[variant2_MAFindex,"Gene_name"] &
#                          MAF[,ref_column] %in% c("A","T","G","C") &
#                          MAF[,alt_column] %in% c("A","T","G","C"),],
#         gene = MAF[variant2_MAFindex,"Gene_name"],
#         gene_mut_rate = mutrates_list,
#         trinuc_proportion_matrix = trinuc_proportion_matrix,
#         gene_trinuc_comp = gene_trinuc_comp,
#         gene_refcds = RefCDS[[gene]],
#         tumor_subsets = tumors,subset_col=subset_col)




#     these_selection_results <- c(variant1,variant2,NA,NA,NA,NA)
#     names(these_selection_results) <- c("Variant_1", "Variant_2", "Gamma_1", "Gamma_2",
#                                         "Gamma_1_2background", "Gamma_2_1background")
#     these_selection_results[3:6] <- cancereffectsizeR::optimize_gamma_epistasis(
#       MAF_input1=MAF[MAF[,"Gene_name"] == MAF[variant1_MAFindex,"Gene_name"] &
#                        MAF[,ref_column] %in% c("A","T","G","C") &
#                        MAF[,alt_column] %in% c("A","T","G","C"),],
#       MAF_input2=MAF[MAF[,"Gene_name"] == MAF[variant2_MAFindex,"Gene_name"] &
#                        MAF[,ref_column] %in% c("A","T","G","C") &
#                        MAF[,alt_column] %in% c("A","T","G","C"),],
#       all_tumors=tumors,
#       gene1=MAF[variant1_MAFindex,"Gene_name"],
#       gene2=MAF[variant2_MAFindex,"Gene_name"],
#       variant1= MAF[variant1_MAFindex,"unique_variant_ID_AA"],
#       variant2= MAF[variant2_MAFindex,"unique_variant_ID_AA"],
#       specific_mut_rates1=these_mutation_rates1$mutation_rate_matrix,
#       specific_mut_rates2=these_mutation_rates2$mutation_rate_matrix)
#     return(these_selection_results)
#   }
#   selection_results <- pbapply::pblapply(selection_epistasis_results_list, get_gene_results_epistasis, cl = cores)
#   return(selection_results)
# }