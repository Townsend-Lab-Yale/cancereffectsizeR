

load("after_step_6.RData")

# script to test all genes in the LUAD dataset.
source("R/mutation_rate_calc.R")
source("R/selection_intensity_calc_function_ML.R")
selection_results <- vector("list",length = length(unique(MAF_input$Gene_name)))
names(selection_results) <- unique(MAF_input$Gene_name)

for(i in 1:length(unique(MAF_input$Gene_name))){


 these_mutation_rates <-  mutation_rate_calc(MAF = MAF_input,gene = unique(MAF_input$Gene_name)[i],gene_mut_rate = mutrates,tumor_trinucs = trinuc_proportion_matrix)

 these_selection_results <- matrix(nrow=ncol(these_mutation_rates), ncol=2,data=NA)
 rownames(these_selection_results) <- colnames(these_mutation_rates); colnames(these_selection_results) <- c("variant","selection_intensity")

 for(j in 1:nrow(these_selection_results)){
  these_selection_results[j,] <- c(colnames(these_mutation_rates)[j] , optimize_gamma(MAF_input=MAF_input, all_tumors=tumors, gene=unique(MAF_input$Gene_name)[i], variant=colnames(these_mutation_rates)[j], specific_mut_rates=these_mutation_rates))
 }

 selection_results[[i]] <- list(gene_name=unique(MAF_input$Gene_name)[i],selection_results=these_selection_results)



 print(round(i/length(unique(MAF_input$Gene_name)),2))
 print(unique(MAF_input$Gene_name)[i])

if(i %% 100 == 0){
  save(selection_results, file="selection_results_LUAD_initial_check.RData")
}

}

save(selection_results, file="selection_results_LUAD_initial_check.RData")


