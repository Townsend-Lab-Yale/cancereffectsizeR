#' SNV effect size calculation function
#'
#' Function that determines mutation rate at the gene- and trinucleotide- level
#' and then calculates cancer effect size.
#' Works best on datasets with a large number of tumors, each with a large
#' number of SNV (The \code{deconstructSigs} package, used here to calculate
#' trinucleotide signatues, suggests at least 50 SNV per tumor used in the
#' signature calculation). Requres data to be in hg coordinates, you can use
#' the \code{hg_converter()} function provided with this package to convert
#' from one coordinate system to another.
#'
#'
#' @param cesa CESAnalysis object
#' @param genes_for_effect_size_analysis genes to calculate effect sizes within. If unspecified, defaults to "all".
#' @param cores number of cores to use
#' @param tumor_specific_rate_choice weights tumor-specific rates by their relative proportional substitution count (not recommended)
#' @param trinuc_all_tumors Calculates trinucleotide signatures within all tumors (even those with < 50 variants)
#'
#' @export





effect_size_SNV <- function(
                            cesa = NULL,
                            covariate_file=NULL,
                            genes_for_effect_size_analysis = "all",
                            cores = 1,
                            tumor_specific_rate_choice = F,
                            trinuc_all_tumors = T,
                            epistasis_top_prev_number = NULL,
                            epistasis_gene_level = F,
                            full_gene_epistasis_lower_optim = 1e-3,
                            full_gene_epistasis_upper_optim=1e9,
                            full_gene_epistasis_fnscale=-1e-16,
                            q_threshold_for_gene_level=0.1,
                            trinuc_algorithm_choice="weighted"){



  MAF = cesa@maf
  if(!is.null(subset_col) & !is.null(epistasis_top_prev_number)){
    stop("You can either measure selection of a linear evolving system
        (i.e., you provide `subset_col` with a value), or
        you can measure selection under epistasis
        (i.e., you provide `epistasis_top_prev_number` with a value...
        but, you cannot do both, yet. Please rerun using valid inputs.")
  }



  version_message <- "PLEASE DOWNLOAD v0.1.0, USING INSTRUCTIONS HERE: \nhttps://github.com/Townsend-Lab-Yale/cancereffectsizeR/blob/master/user_guide/cancereffectsizeR_user_guide.md \nOR devtools::install_github(`Townsend-Lab-Yale/cancereffectsizeR@0.1.0`) \nOR https://github.com/Townsend-Lab-Yale/cancereffectsizeR/releases/tag/0.1.0
\nEverything after v0.1.0 is experimental. v0.1.0 is associated with this publication: \nhttps://doi.org/10.1093/jnci/djy168"

  warning(version_message)
  message(version_message)

  message("Calculating trinucleotide composition and signatures...")
  cesa = cancereffectsizeR::trinucleotide_mutation_weights(cesa, algorithm_choice = trinuc_algorithm_choice)

  # Calculate gene-level mutation rates using dNdScv
  cesa = cancereffectsizeR::gene_level_mutation_rates(cesa, covariate_file)

  # Assign genes to MAF, keeping assignments consistent with dndscv when possible
  cesa = cancereffectsizeR::annotate_gene_maf(cesa)

  # for each SNV, calculate the gene- and tumor- and mutation-specific mutation rate
  cesa = cancereffectsizeR::snv_mutation_rates(cesa)

  message("Calculating selection intensity...")

  genes_in_dataset = unique(cesa@annotated.maf$Gene_name)
  if(all(genes_for_effect_size_analysis=="all")) {
    genes_to_analyze <- genes_in_dataset
  } else{
    genes_to_analyze <- genes_for_effect_size_analysis[genes_for_effect_size_analysis %in% genes_in_dataset]
  }

  if(is.null(epistasis_top_prev_number)){


  } else {


    if(epistasis_gene_level){ # if epistasis calc is on gene level ----

      CDS_sizes <- sapply(RefCDS_our_genes,function(x) x$CDS_length)
      names(CDS_sizes) <- names(RefCDS_our_genes)


      # want to prioritize genes with recurrent variants
      MAF$gene_AA_unique <- paste(MAF$Gene_name, MAF$unique_variant_ID_AA)
      gene_AA_table <- table(MAF$gene_AA_unique)
      MAF$gene_AA_tally <- gene_AA_table[MAF$gene_AA_unique]
      MAF$gene_has_recurrent_variant <- F
      if(length(which(MAF$gene_AA_tally>1))>0){
        genes_with_recurrent_variants <- MAF$Gene_name[which(MAF$gene_AA_tally>1)]
        MAF$gene_has_recurrent_variant[which(MAF$Gene_name %in% genes_with_recurrent_variants)] <- T
      }else{
        stop("There are no recurrent variants in the dataset, so you cannot
            run the gene-by-gene epistasis analysis")
      }

      MAF_all_sure <- MAF
      MAF_all_sure <- MAF_all_sure[MAF_all_sure$gene_has_recurrent_variant,]

      # list of tumors that contain recurrent variants per gene
      tumors_per_gene <- vector(mode = "list",length = length(unique(MAF_all_sure$Gene_name)))
      recurrent_gene_list <- unique(MAF_all_sure$Gene_name)
      names(tumors_per_gene) <- recurrent_gene_list
      for(gene in 1:length(tumors_per_gene)){
        tumors_per_gene[[gene]] <- MAF_all_sure[which(MAF_all_sure$Gene_name == recurrent_gene_list[gene] & MAF_all_sure$gene_AA_tally>1),sample_ID_column]
      }

      # see which genes have recurrent variants >1 tumors

      selection_epistasis_results_list <- list()

      for(gene in 1:(length(tumors_per_gene)-1)){

        for(other_genes in (gene+1):length(tumors_per_gene)){

          # if there are more than one tumors with shared recurrent variants among the genes
          if(length(which(tumors_per_gene[[gene]] %in% tumors_per_gene[[other_genes]]))>1){
            selection_epistasis_results_list[[(length(selection_epistasis_results_list)+1)]] <- c(recurrent_gene_list[gene],recurrent_gene_list[other_genes])
          }

        }

      }

      get_gene_results_epistasis_bygene <- function(variant_combo_list) {

        # print(variant_combo_list)

        variant1 <- variant_combo_list[1]
        variant2 <- variant_combo_list[2]

        variant1_MAFindex <- which(MAF$identifier==variant1)[1]
        variant2_MAFindex <- which(MAF$identifier==variant2)[1]



        MAF_input1=MAF[MAF[,"Gene_name"] == variant1 &
                         MAF[,ref_column] %in% c("A","T","G","C") &
                         MAF[,alt_column] %in% c("A","T","G","C"),]

        MAF_input2=MAF[MAF[,"Gene_name"] == variant2 &
                         MAF[,ref_column] %in% c("A","T","G","C") &
                         MAF[,alt_column] %in% c("A","T","G","C"),]


        variant_freq_1 <- table(MAF_input1$unique_variant_ID_AA)
        variant_freq_2 <- table(MAF_input2$unique_variant_ID_AA)

        # only run the selection algorithm if there are 2 or more tumors with
        # recurrent variants of each gene present.
        # if(length(tumors_with_both_mutated) > 1) {


        these_mutation_rates1 <-
          cancereffectsizeR::mutation_rate_calc(
            this_MAF = MAF[MAF[,"Gene_name"] == variant1 &
                             MAF[,ref_column] %in% c("A","T","G","C") &
                             MAF[,alt_column] %in% c("A","T","G","C"),],
            gene = variant1,
            gene_mut_rate = mutrates_list,
            trinuc_proportion_matrix = trinucleotide_weights$trinuc_proportion_matrix,
            gene_trinuc_comp = gene_trinuc_comp,
            RefCDS = RefCDS_our_genes,
            relative_substitution_rate=relative_substitution_rate,
            tumor_specific_rate=tumor_specific_rate_choice,
            tumor_subsets = tumors,
            subset_col=subset_col)

        these_mutation_rates2 <-
          cancereffectsizeR::mutation_rate_calc(
            this_MAF = MAF[MAF[,"Gene_name"] == variant2 &
                             MAF[,ref_column] %in% c("A","T","G","C") &
                             MAF[,alt_column] %in% c("A","T","G","C"),],
            gene = variant2,
            gene_mut_rate = mutrates_list,
            trinuc_proportion_matrix = trinucleotide_weights$trinuc_proportion_matrix,
            gene_trinuc_comp = gene_trinuc_comp,
            RefCDS = RefCDS_our_genes,
            relative_substitution_rate=relative_substitution_rate,
            tumor_specific_rate=tumor_specific_rate_choice,
            tumor_subsets = tumors,subset_col=subset_col)




        # since we are looking at selection at the gene level,
        # we only consider selection at sites that are recurrently
        # substituted.


        these_mutation_rates1$mutation_rate_matrix <- as.matrix(these_mutation_rates1$mutation_rate_matrix[,colnames(these_mutation_rates1$mutation_rate_matrix) %in% names(these_mutation_rates1$variant_freq$no_subset[these_mutation_rates1$variant_freq$no_subset>1])])

        these_mutation_rates2$mutation_rate_matrix <- as.matrix(these_mutation_rates2$mutation_rate_matrix[,colnames(these_mutation_rates2$mutation_rate_matrix) %in% names(these_mutation_rates2$variant_freq$no_subset[these_mutation_rates2$variant_freq$no_subset>1])])

        these_selection_results <- c(variant1,variant2,NA,NA,NA,NA)
        names(these_selection_results) <- c("Variant_1",
                                            "Variant_2",
                                            "Gamma_1",
                                            "Gamma_2",
                                            "Gamma_1_2background",
                                            "Gamma_2_1background")






        these_selection_results[3:6] <-
          cancereffectsizeR::optimize_gamma_epistasis_full_gene(
            MAF_input1=MAF[MAF[,"Gene_name"] == variant1 &
                             MAF[,ref_column] %in% c("A","T","G","C") &
                             MAF[,alt_column] %in% c("A","T","G","C"),],
            MAF_input2=MAF[MAF[,"Gene_name"] == variant2 &
                             MAF[,ref_column] %in% c("A","T","G","C") &
                             MAF[,alt_column] %in% c("A","T","G","C"),],
            all_tumors=tumors,
            gene1=variant1,
            gene2=variant2,
            specific_mut_rates1=these_mutation_rates1$mutation_rate_matrix,
            specific_mut_rates2=these_mutation_rates2$mutation_rate_matrix,
            variant_freq_1 = these_mutation_rates1$variant_freq$no_subset,
            variant_freq_2 = these_mutation_rates2$variant_freq$no_subset,
            full_gene_epistasis_lower_optim=full_gene_epistasis_lower_optim,
            full_gene_epistasis_upper_optim=full_gene_epistasis_upper_optim,
            full_gene_epistasis_fnscale=full_gene_epistasis_fnscale)




        return(these_selection_results)
      }

      selection_results <- parallel::mclapply(selection_epistasis_results_list, get_gene_results_epistasis_bygene, mc.cores = cores)

    } else{

      # find the prevalences from the MAF data
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

      # function to calculate selection results

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
            trinuc_proportion_matrix = trinucleotide_weights$trinuc_proportion_matrix,
            gene_trinuc_comp = gene_trinuc_comp,
            RefCDS = RefCDS_our_genes,
            relative_substitution_rate=relative_substitution_rate,
            tumor_specific_rate=tumor_specific_rate_choice,
            tumor_subsets = tumors,subset_col=subset_col)

        these_mutation_rates2 <-
          cancereffectsizeR::mutation_rate_calc(
            this_MAF = MAF[MAF[,"Gene_name"] == MAF[variant2_MAFindex,"Gene_name"] &
                             MAF[,ref_column] %in% c("A","T","G","C") &
                             MAF[,alt_column] %in% c("A","T","G","C"),],
            gene = MAF[variant2_MAFindex,"Gene_name"],
            gene_mut_rate = mutrates_list,
            trinuc_proportion_matrix = trinucleotide_weights$trinuc_proportion_matrix,
            gene_trinuc_comp = gene_trinuc_comp,
            RefCDS = RefCDS_our_genes,
            relative_substitution_rate=relative_substitution_rate,
            tumor_specific_rate=tumor_specific_rate_choice,
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

      selection_results <- parallel::mclapply(selection_epistasis_results_list, get_gene_results_epistasis, mc.cores = cores)
    }
  }




  return(list(selection_output=selection_results,
              mutation_rates=mutrates_list,
              trinuc_data=list(trinuc_proportion_matrix=trinucleotide_weights$trinuc_proportion_matrix,
                               signatures_output=trinucleotide_weights$signatures_output_list),
              dndscvout=dndscv_out_list,
              MAF=MAF))

}



