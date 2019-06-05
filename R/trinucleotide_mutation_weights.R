#' Trinucleotide mutation weights
#'
#' @param MAF
#' @param sample_ID_column
#' @param chr_column
#' @param pos_column
#' @param ref_column
#' @param alt_column
#' @param algorithm_choice The choice of the algorithm used in determining trinucleotide weights.
#'   Defaults to `weighted`, where all tumors with >= 50 substitutions have mutation
#'   rates directly from deconstructSigs, and tumors with < 50 substitutions have that rate weighted
#'   by the average of all tumors >= 50 substitutions relative to the substitution count < 50 in those tumors.
#'   Other options currently include `all_average`, where all tumors have a mutation rate equal to the
#'   average of all tumors >= 50 substitutions, `nearest_neighbor` where all tumors <50 substitutions
#'   have the mutation rate of the tumor with >=50 substitutions closest to the mutational
#'   distribution profile of the low substitution tumor as determined by a distance matrix, and
#'   `all_calculated`, where all tumors have the mutation rate from the deconstructSigs output regardless of substitution number.
#' @param remove_recurrent If TRUE, removes recurrent variants prior to calculating trinucleotide signatures and weights. If FALSE, uses all variants.
#'
#' @return
#' @export
#' @import deconstructSigs
#'
#' @examples
#'
#'
#'
trinucleotide_mutation_weights <- function(MAF,
                                           sample_ID_column="Unique_patient_identifier",
                                           chr_column = "Chromosome",
                                           pos_column = "Start_Position",
                                           ref_column = "Reference_Allele",
                                           alt_column = "Tumor_allele",
                                           algorithm_choice = "weighted",
                                           remove_recurrent = TRUE){



  #' only want to find the signatures that reflect the underlying mutational
  #' processes within the tumor, i.e., those that reach detectable levels neutrally
  #' so we detect the signatues from non-recurrently substituted variants


  #'TODO: it is possible that the only mutation in a tumor is a recurrent variant,
  #' and if we remove this variant, then the tumor is lost from the mutation
  #' rate calculation. Should this be added back in if algorithm_choice
  #' is something involving averages?


  # pre-process ----

  if(remove_recurrent){
    MAF_input_deconstructSigs_preprocessed <-
      cancereffectsizeR::deconstructSigs_input_preprocess(MAF = MAF,
                                                          sample_ID_column = sample_ID_column,
                                                          chr_column = chr_column,
                                                          pos_column = pos_column,
                                                          ref_column = ref_column,
                                                          alt_column = alt_column)
  }else{
    MAF[,chr_column] <- paste("chr",trimws(MAF[,chr_column]),sep="")
    MAF_input_deconstructSigs_preprocessed <- MAF
  }


  tumors_with_a_mutation_rate <- unique(MAF_input_deconstructSigs_preprocessed[,sample_ID_column])



  substitution_counts <- table(MAF_input_deconstructSigs_preprocessed[,sample_ID_column])
  tumors_with_50_or_more <- names(which(substitution_counts>=50))
  tumors_with_less_than_50 <- setdiff(MAF_input_deconstructSigs_preprocessed[,sample_ID_column],tumors_with_50_or_more)

  # relative total substitution rate in all tumors
  median_substitutions <- median(as.numeric(substitution_counts))

  relative_substitution_rate <- substitution_counts/median_substitutions



  trinuc_breakdown_per_tumor <- deconstructSigs::mut.to.sigs.input(mut.ref =
                                                                     MAF_input_deconstructSigs_preprocessed,
                                                                   sample.id = sample_ID_column,
                                                                   chr = chr_column,
                                                                   pos = pos_column,
                                                                   ref = ref_column,
                                                                   alt = alt_column)


  deconstructSigs_trinuc_string <- colnames(trinuc_breakdown_per_tumor)

  # rm(MAF_input_deconstructSigs_preprocessed)

  # message("Building trinuc_proportion_matrix")


  trinuc_proportion_matrix <- matrix(data = NA,
                                     nrow = nrow(trinuc_breakdown_per_tumor),
                                     ncol = ncol(trinuc_breakdown_per_tumor))
  rownames(trinuc_proportion_matrix) <- rownames(trinuc_breakdown_per_tumor)
  colnames(trinuc_proportion_matrix) <- colnames(trinuc_breakdown_per_tumor)

  signatures_output_list <- vector(mode = "list",length = nrow(trinuc_breakdown_per_tumor))
  names(signatures_output_list) <- rownames(trinuc_breakdown_per_tumor)


  # algorithms ----


  if(algorithm_choice == "weighted"){

    for(tumor_name in 1:nrow(trinuc_breakdown_per_tumor)){

      signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                            signatures.ref = signatures.cosmic,
                                                            sample.id = rownames(trinuc_breakdown_per_tumor)[tumor_name],
                                                            contexts.needed = TRUE,
                                                            tri.counts.method = 'exome2genome')

      signatures_output_list[[tumor_name]] <- list(signatures_output = signatures_output,
                                                   substitution_count = length(which(MAF[,sample_ID_column] == rownames(trinuc_breakdown_per_tumor)[tumor_name])))


      trinuc_proportion_matrix[rownames(trinuc_breakdown_per_tumor)[tumor_name],] <- signatures_output$product/sum( signatures_output$product) #need it to sum to 1.

      # Not all trinuc weights in the cosmic dataset are nonzero for certain signatures
      # This leads to the rare occasion where a certain combination of signatues leads to ZERO
      # rate for particular trinucleotide contexts.
      # True rate is nonzero, as we do see those variants in those tumors, so renormalizing
      # the rates by adding the lowest nonzero rate to all the rates and renormalizing


      if(0 %in% trinuc_proportion_matrix[tumor_name,]){
        # finding the lowest nonzero rate
        lowest_rate <- min(trinuc_proportion_matrix[tumor_name,
                                                    -which(trinuc_proportion_matrix[tumor_name,]==0)])

        # adding it to the rates
        trinuc_proportion_matrix[tumor_name,] <- trinuc_proportion_matrix[tumor_name,] + lowest_rate

        # renormalizing to 1
        trinuc_proportion_matrix[tumor_name,] <-
          trinuc_proportion_matrix[tumor_name,] / sum(trinuc_proportion_matrix[tumor_name,])
      }

    }

    # finding the average of the >= 50 tumors
    # original_signatures <- signatures_output_list
    averaged_product <- 0
    # add up the products ...
    for(tumor_name_index in 1:length(tumors_with_50_or_more)){
      averaged_product <- averaged_product +
        trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name_index],]
    }

    # ... and divide by the length
    averaged_product <- averaged_product/length(tumors_with_50_or_more)



    # need to do the same for weights
    weight_matrix <- matrix(data = NA,
                            nrow = length(signatures_output_list),
                            ncol = ncol(signatures_output_list[[1]]$signatures_output$weights))
    rownames(weight_matrix) <- names(signatures_output_list)
    colnames(weight_matrix) <- colnames(signatures_output_list[[1]]$signatures_output$weights)

    for(weight_row in 1:nrow(weight_matrix)){
     weight_matrix[weight_row,] <-  as.numeric(signatures_output_list[[rownames(weight_matrix)[weight_row]]]$signatures_output$weights)
    }

    averaged_weight <- 0
    for(tumor_name_index in 1:length(tumors_with_50_or_more)){
      averaged_weight <- averaged_weight +
        weight_matrix[tumors_with_50_or_more[tumor_name_index],]
    }

    averaged_weight <- averaged_weight/length(tumors_with_50_or_more)
    #TODO: should this be forced to sum() to 1?



    # for tumors < 50 substitutions, weight their found signature product by the average
    # of the signatures over 50 relative to substitution load.
    if(length(tumors_with_less_than_50)>0){
      for(tumor_name_index in 1:length(tumors_with_less_than_50)){

        trinuc_proportion_matrix[tumors_with_less_than_50[tumor_name_index],]  <-
          ((substitution_counts[tumors_with_less_than_50[tumor_name_index]]/50) *
             trinuc_proportion_matrix[tumors_with_less_than_50[tumor_name_index],]) +
          (((50-substitution_counts[tumors_with_less_than_50[tumor_name_index]])/50) *
             averaged_product)

        signatures_output_list[[tumors_with_less_than_50[tumor_name_index]]]$signatures_output$weights <-
          ((substitution_counts[tumors_with_less_than_50[tumor_name_index]]/50) *
             signatures_output_list[[tumors_with_less_than_50[tumor_name_index]]]$signatures_output$weights) +
          (((50-substitution_counts[tumors_with_less_than_50[tumor_name_index]])/50) * averaged_weight)

        signatures_output_list[[tumors_with_less_than_50[tumor_name_index]]]$signatures_output$product <-
          trinuc_proportion_matrix[tumors_with_less_than_50[tumor_name_index],]



      }
    }


    # If a tumor had all variants removed in preprocessing
    # (meaning it only had recurrent variants)
    # then add back in the mutation rates as the average of the > 50
    if(length(signatures_output_list) < length(unique(MAF[,sample_ID_column]))){
      tumors_to_add <- unique(MAF[,sample_ID_column])[which(!unique(MAF[,sample_ID_column]) %in% names(signatures_output_list))]

      matrix_to_add <- matrix(data = NA, nrow = length(tumors_to_add), ncol = ncol(trinuc_breakdown_per_tumor))
      rownames(matrix_to_add) <- tumors_to_add
      colnames(matrix_to_add) <- colnames(trinuc_breakdown_per_tumor)

      for(matrix_row in 1:nrow(matrix_to_add)){
       matrix_to_add[matrix_row,] <- as.numeric(averaged_product)
       signatures_output_list[[rownames(matrix_to_add)]]$signatures_output$product <- matrix(data=averaged_product,nrow=1)
       rownames(signatures_output_list[[rownames(matrix_to_add)]]$signatures_output$product) <- rownames(matrix_to_add)
       colnames(signatures_output_list[[rownames(matrix_to_add)]]$signatures_output$product) <- names(averaged_product)
       signatures_output_list[[rownames(matrix_to_add)]]$signatures_output$weights <- matrix(data=averaged_weight, nrow=1)
       rownames(signatures_output_list[[rownames(matrix_to_add)]]$signatures_output$weights) <- rownames(matrix_to_add)
       colnames(signatures_output_list[[rownames(matrix_to_add)]]$signatures_output$weights) <- names(averaged_weight)

      }

      trinuc_proportion_matrix <- rbind(trinuc_proportion_matrix, matrix_to_add)
      tumors_with_a_mutation_rate <- rownames(trinuc_proportion_matrix)





    }



  }


  # takes the average of all profiles with >= 50 substitutions
  if(algorithm_choice == "all_average"){

    for(tumor_name in 1:nrow(trinuc_breakdown_per_tumor)){

      signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                            signatures.ref = signatures.cosmic,
                                                            sample.id = rownames(trinuc_breakdown_per_tumor)[tumor_name],
                                                            contexts.needed = TRUE,
                                                            tri.counts.method = 'exome2genome')

      signatures_output_list[[tumor_name]] <- list(signatures_output = signatures_output,
                                                   substitution_count = length(which(MAF[,sample_ID_column] == rownames(trinuc_breakdown_per_tumor)[tumor_name])))


      trinuc_proportion_matrix[rownames(trinuc_breakdown_per_tumor)[tumor_name],] <- signatures_output$product/sum( signatures_output$product) #need it to sum to 1.

      # Not all trinuc weights in the cosmic dataset are nonzero for certain signatures
      # This leads to the rare occasion where a certain combination of signatues leads to ZERO
      # rate for particular trinucleotide contexts.
      # True rate is nonzero, as we do see those variants in those tumors, so renormalizing
      # the rates by adding the lowest nonzero rate to all the rates and renormalizing


      if(0 %in% trinuc_proportion_matrix[tumor_name,]){
        # finding the lowest nonzero rate
        lowest_rate <- min(trinuc_proportion_matrix[tumor_name,
                                                    -which(trinuc_proportion_matrix[tumor_name,]==0)])

        # adding it to the rates
        trinuc_proportion_matrix[tumor_name,] <- trinuc_proportion_matrix[tumor_name,] + lowest_rate

        # renormalizing to 1
        trinuc_proportion_matrix[tumor_name,] <-
          trinuc_proportion_matrix[tumor_name,] / sum(trinuc_proportion_matrix[tumor_name,])
      }

    }


    # finding the average of the >= 50 tumors

    averaged_product <- 0
    # add up the weights ...
    for(tumor_name_index in 1:length(tumors_with_50_or_more)){
      averaged_product <- averaged_product + trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name_index],]
    }

    # ... and divide by the length
    averaged_product <- averaged_product/length(tumors_with_50_or_more)

    for(this_row in 1:nrow(trinuc_proportion_matrix)){
      trinuc_proportion_matrix[this_row,] <- averaged_product
    }

  }



  if(algorithm_choice == "nearest_neighbor"){
    for(tumor_name in 1:length(tumors_with_50_or_more)){
      signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                            signatures.ref = signatures.cosmic,
                                                            sample.id = tumors_with_50_or_more[tumor_name],
                                                            contexts.needed = TRUE,
                                                            tri.counts.method = 'exome2genome')

      signatures_output_list[[tumors_with_50_or_more[tumor_name]]] <-
        list(signatures_output = signatures_output,
             substitution_count = length(which(MAF[,sample_ID_column] == tumors_with_50_or_more[tumor_name])))


      trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] <- signatures_output$product/sum( signatures_output$product) #need it to sum to 1.

      # Not all trinuc weights in the cosmic dataset are nonzero for certain signatures
      # This leads to the rare occasion where a certain combination of signatues leads to ZERO
      # rate for particular trinucleotide contexts.
      # True rate is nonzero, as we do see those variants in those tumors, so renormalizing
      # the rates by adding the lowest nonzero rate to all the rates and renormalizing


      if(0 %in% trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],]){
        # finding the lowest nonzero rate
        lowest_rate <- min(trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],
                                                    -which(trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],]==0)])

        # adding it to the rates
        trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] <- trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] + lowest_rate

        # renormalizing to 1
        trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] <-
          trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] / sum(trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],])
      }

    }

    # 2. Find nearest neighbor to tumors with < 50 mutations, assign identical weights as neighbor ----
    # message(head(trinuc_proportion_matrix))

    # message("should have printed")

    distance_matrix <- as.matrix(dist(trinuc_breakdown_per_tumor))

    for(tumor_name in 1:length(tumors_with_less_than_50)){
      #find closest tumor that have over 50 mutations
      closest_tumor <- names(sort(distance_matrix[tumors_with_less_than_50[tumor_name],tumors_with_50_or_more]))[1]

      trinuc_proportion_matrix[tumors_with_less_than_50[tumor_name],] <- trinuc_proportion_matrix[closest_tumor,]

      signatures_output_list[[tumors_with_less_than_50[tumor_name]]] <- list(signatures_output = signatures_output_list[[closest_tumor]]$signatures_output, substitution_count = length(which(MAF[,sample_ID_column] == tumors_with_less_than_50[tumor_name])), tumor_signatures_used = closest_tumor)


    }

    # rm(distance_matrix)
    # now, trinuc_proportion_matrix has the proportion of all trinucs in every tumor.

    # collect the garbage
    # gc()

  }



  if(algorithm_choice == "all_calculated"){

    for(tumor_name in 1:nrow(trinuc_breakdown_per_tumor)){

      signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                            signatures.ref = signatures.cosmic,
                                                            sample.id = rownames(trinuc_breakdown_per_tumor)[tumor_name],
                                                            contexts.needed = TRUE,
                                                            tri.counts.method = 'exome2genome')

      signatures_output_list[[tumor_name]] <- list(signatures_output = signatures_output,
                                                   substitution_count = length(which(MAF[,sample_ID_column] == rownames(trinuc_breakdown_per_tumor)[tumor_name])))


      trinuc_proportion_matrix[rownames(trinuc_breakdown_per_tumor)[tumor_name],] <- signatures_output$product/sum( signatures_output$product) #need it to sum to 1.

      # Not all trinuc weights in the cosmic dataset are nonzero for certain signatures
      # This leads to the rare occasion where a certain combination of signatues leads to ZERO
      # rate for particular trinucleotide contexts.
      # True rate is nonzero, as we do see those variants in those tumors, so renormalizing
      # the rates by adding the lowest nonzero rate to all the rates and renormalizing


      if(0 %in% trinuc_proportion_matrix[tumor_name,]){
        # finding the lowest nonzero rate
        lowest_rate <- min(trinuc_proportion_matrix[tumor_name,
                                                    -which(trinuc_proportion_matrix[tumor_name,]==0)])

        # adding it to the rates
        trinuc_proportion_matrix[tumor_name,] <- trinuc_proportion_matrix[tumor_name,] + lowest_rate

        # renormalizing to 1
        trinuc_proportion_matrix[tumor_name,] <-
          trinuc_proportion_matrix[tumor_name,] / sum(trinuc_proportion_matrix[tumor_name,])
      }

    }



  }




  return(list(tumors_with_a_mutation_rate=tumors_with_a_mutation_rate,
              trinuc_proportion_matrix=trinuc_proportion_matrix,
              signatures_output_list=signatures_output_list,
              algorithm_choice=algorithm_choice)
  )


}
