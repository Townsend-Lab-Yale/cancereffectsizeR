#' Single-stage SNV effect size analysis (gets called by effect_size_SNV)
get_gene_results <- function(gene_to_analyze, cesa, gene_mafs, gene_trinuc_comp, 
							 tumor_specific_rate_choice, all_tumors) {
	mutrates_list = cesa@mutrates_list
	relative_substitution_rates = cesa@relative_substitution_rates
	trinuc_proportion_matrix = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix
	RefCDS = cesa@refcds_data
	progressions = cesa@progressions
	these_mutation_rates <-
			cancereffectsizeR::mutation_rate_calc(
			this_MAF = gene_mafs[[gene_to_analyze]][["maf"]],
			gene = gene_to_analyze,
			gene_mut_rate = mutrates_list,
			trinuc_proportion_matrix = trinuc_proportion_matrix,
			gene_trinuc_comp = gene_trinuc_comp,
			RefCDS = RefCDS,
			all_tumors = all_tumors,
			progressions = progressions,
			relative_substitution_rates=relative_substitution_rates,
			tumor_specific_rate=tumor_specific_rate_choice)




	these_selection_results <- dplyr::tibble(variant = colnames(these_mutation_rates$mutation_rate_matrix),
											selection_intensity = vector(mode = "list",
											length = ncol(these_mutation_rates$mutation_rate_matrix)),
											unsure_gene_name=NA,
											variant_freq=vector(mode = "list", length = ncol(these_mutation_rates$mutation_rate_matrix)),
											unique_variant_ID=NA,
											loglikelihood=NA)



	num_selection_results = nrow(these_selection_results)
	for(j in 1:num_selection_results){
		optimization_output <- cancereffectsizeR::optimize_gamma(
			MAF_input= gene_mafs[[gene_to_analyze]][["maf"]],
			all_tumors=all_tumors,
			progressions = progressions,
			gene=gene_to_analyze,
			variant=colnames(these_mutation_rates$mutation_rate_matrix)[j],
			specific_mut_rates=these_mutation_rates$mutation_rate_matrix)

		these_selection_results[j,c("selection_intensity")][[1]] <-
			list(optimization_output$par)


		names(these_selection_results[j,c("selection_intensity")][[1]][[1]]) <- names(progressions@order)

		these_selection_results[j,"loglikelihood"] <- optimization_output$value


		freq_vec <- NULL
		for(this_level in 1:length(these_mutation_rates$variant_freq)){
			freq_vec <- c(freq_vec,these_mutation_rates$variant_freq[[this_level]][as.character(these_selection_results[j,"variant"])])
		} 
		these_selection_results[j,"variant_freq"][[1]] <- list(freq_vec)
		names(these_selection_results[j,"variant_freq"][[1]][[1]]) <- progressions@order

	}

	these_selection_results[,"unsure_gene_name"] <- these_mutation_rates$unsure_genes_vec
	these_selection_results[,"unique_variant_ID"] <- these_mutation_rates$unique_variant_ID

	print(gene_to_analyze)
	return(list(gene_name=gene_to_analyze, RefCDS[[gene_to_analyze]]$gene_id,selection_results=these_selection_results))
}


## Gets called by effect_size_SNV for gene-level epistasis
epistasis_gene_level = function(genes_to_analyze, MAF, trinuc_proportion_matrix) {

	CDS_sizes <- sapply(cesa@refcds_data, function(x) x$CDS_length)
	names(CDS_sizes) <- names(cesa@refcds_data)


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
			trinuc_proportion_matrix = trinuc_proportion_matrix,
			gene_trinuc_comp = gene_trinuc_comp,
			RefCDS = cesa@refcds_data,
			relative_substitution_rates=relative_substitution_rates,
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
			trinuc_proportion_matrix = trinuc_proportion_matrix,
			gene_trinuc_comp = gene_trinuc_comp,
			RefCDS = cesa@refcds_data,
			relative_substitution_rates=relative_substitution_rates,
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
			all_tumors=all_tumors,
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
			RefCDS = RefCDS_our_genes,
			relative_substitution_rates=relative_substitution_rates,
			tumor_specific_rate=tumor_specific_rate_choice,
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
			RefCDS = RefCDS_our_genes,
			relative_substitution_rates=relative_substitution_rates,
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
	return(selection_results)
}


#' @name effect_size_SNV
#' @param cesa CESAnalysis object
#' @param gene which genes to calculate effect sizes within; defaults to all genes in data set
#' @param analysis choose SNV, gene-level-epistasis, or top-epistasis
#' @param cores number of cores to use
#' @param ignore_progression_stages don't calculate stage-specific selection intensities even if stage data is present
#' @param tumor_specific_rate_choice weights tumor-specific rates by their relative proportional substitution count (not recommended)
#' @return CESAnalysis object with selection results added for the chosen analysis
#' @export
effect_size_SNV <- function(
		cesa = NULL,
		genes = "all",
		analysis = c("SNV", "gene-level-epistasis", "top-epistasis"),
		epistasis_top_prev_number = NULL,
		cores = 1,
		ignore_progression_stages = F,
		tumor_specific_rate_choice = F,
		full_gene_epistasis_lower_optim = 1e-3,
		full_gene_epistasis_upper_optim=1e9,
		full_gene_epistasis_fnscale=-1e-16,
		q_threshold_for_gene_level=0.1) {

	# validate analysis choice
	analysis = match.arg(analysis)

	if (analysis != "top-epistasis" && ! is.null(epistasis_top_prev_number)) {
		stop("Error: epistasis_top_prev_number is only used for \"top-epistasis\" analysis")
	} else if(! is.null(epistasis_top_prev_number)) {
		stop("Error: epistasis_top_prev_number must be defined for top-epistasis analysis")
	}

	# can't combine epistasis analysis with multi-stage yet
	if (analysis != "SNV" && length(cesa@progressions@order) > 1) {
		message(paste("Warning: Tumor progression stage-specific analysis cannot be combined with an epistasis ",
					"analysis yet, so SNV selection intensities will be treated as constant across stages."))
	} else if(length(cesa@progressions@order) > 1 && ! ignore_progression_stages) {
		message(paste("Note: Tumors are annotated with progression stages, so SNV selection intensities will be estimated ",
			"for each progression. To ignore stage information, re-run with ignore_progression_stages=TRUE"))
	}

	# using the "SNV" genes
	snv.maf = cesa@annotated.snv.maf
	genes_in_dataset = unique(snv.maf$Gene_name) 
	if(genes[1] =="all") {
		genes_to_analyze <- genes_in_dataset
	} else{
		genes_to_analyze <- genes[genes %in% genes_in_dataset]
	}
	
	data("gene_trinuc_comp", package = "cancereffectsizeR")
	names(gene_trinuc_comp) <- sapply(gene_trinuc_comp, function(x) x$gene_name) # formally used names of RefCDS, but they're the same


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
		selection_results <- parallel::mclapply(genes_to_analyze, get_gene_results, mc.cores = cores, cesa = cesa,
										gene_mafs = mafs, gene_trinuc_comp = gene_trinuc_comp,
										tumor_specific_rate_choice = tumor_specific_rate_choice,
										all_tumors = all_tumors)
	} else if(analysis == "gene-level-epistasis") {
		selection_results = epistasis_gene_level()
	} else if(analysis == "top-epistasis") {
		selection_results = top_epistasis()
	}
	cesa@selection_results = selection_results
	return(cesa)
}


	
