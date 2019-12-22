#' Trinucleotide mutation weights
#'
#' Calculates
#' Runs Removes signatures from full potential signatures to minimize "signature bleeding", along the rationale proposed within the manuscript that originally calculated the v3 signature set doi: https://doi.org/10.1101/322859. Use `NULL` to keep all signatures. Signatures must correspond to the signature names in `signature_choice`.
#'
#' @param cesa # CESAnalysis object
#' @param algorithm_choice The choice of the algorithm used in determining trinucleotide weights.
#'   Defaults to `weighted`, where all tumors with >= 50 substitutions have mutation
#'   rates directly from deconstructSigs, and tumors with < 50 substitutions have that rate weighted
#'   by the average of all tumors >= 50 substitutions relative to the substitution count < 50 in those tumors.
#'   Other options currently include `all_average`, where all tumors have a mutation rate equal to the
#'   average of all tumors >= 50 substitutions, `nearest_neighbor` where all tumors <50 substitutions
#'   have the mutation rate of the tumor with >=50 substitutions closest to the mutational
#'   distribution profile of the low substitution tumor as determined by a distance matrix, and
#'   `all_calculated`, where all tumors have the mutation rate from the deconstructSigs output regardless of substitution number.
#' @param remove_recurrent if TRUE (default), removes recurrent variants from signature analysis
#' @param signature_choice "cosmic_v3" (default), "cosmic_v2", or a properly-formatted data frame with of trinucleotide signatures
#' @param v3_artifact_accounting when COSMIC v3 signatures associated with sequencing artifacts are detected, renormalizes to isolate "true" sources of mutational flux.
#' @param signatures_to_remove specify any signatures to exclude from analysis; some signatures automatically get excluded
#'     from COSMIC v3 analyses; set to signatures_to_remove="none" to prevent this behavior (and )
#' @param v3_exome_hypermutation_rules T/F on whether to follow the mutation count rules outlined in https://doi.org/10.1101/322859, the manuscript reported the v3 COSMIC signature set.
#'
#' @return
#' @export
#' @examples
#'
#'
#'
trinucleotide_mutation_weights <- function(cesa,
                                           algorithm_choice = "weighted",
                                           remove_recurrent = TRUE,
                                           signature_choice = "cosmic_v3",
                                           v3_artifact_accounting = TRUE,
                                           v3_exome_hypermutation_rules = TRUE,
                                           signatures_to_remove = "" # cosmic_v3 analysis gets signatures added here later unless "none"
                                           ){

  artifact_accounting = FALSE # gets set to true if v3_artifact_account = TRUE and signature_choice = cosmic_v3
  algorithms = c("weighted", "all_average", "nearest_neighbor", "all_calculated")
  if(! "character" %in% class(algorithm_choice) || ! algorithm_choice %in% algorithms ) {
    stop(paste0("algorithm_choice (for dealing with samples with fewer than 50 variants) must be one of ", paste(algorithms, collapse = ", ")))
  }
  if(! "character" %in% class(signatures_to_remove)) {
    stop("signatures_to_remove should be a character vector")
  }
  bad_sig_choice_msg = "signature_choice should be \"cosmic_v3\", \"cosmic_v2\", or a properly-formatted signatures data.frame"
  if("character" %in% class(signature_choice)) {
    signature_choice = tolower(signature_choice)
    if (length(signature_choice) != 1 || ! signature_choice %in% c("cosmic_v3", "cosmic_v2")) {
      stop(bad_sig_choice_msg)
    }
    if(signature_choice == "cosmic_v3") {
      signatures = signatures_cosmic_May2019
      message("Analyzing tumors for COSMIC v3 (May 2019) SNV mutational signatures....")
      if (v3_artifact_accounting) {
        message("Sequencing artifact signatures will be subtracted to isolate true mutational processes (disable with v3_artifact_accounting=FALSE).")
        artifact_accounting = TRUE
      }
      if (v3_exome_hypermutation_rules) {
        message("Samples with many variants will have special exome hypermutation rules applied (disable with v3_exome_hypermutation_rules=FALSE).")
      }
      if(length(signatures_to_remove) == 1 && signatures_to_remove == "") {
        signatures_to_remove = c("SBS25","SBS31","SBS32","SBS35")
        message("The following signatures will be excluded (to include all signatures, set signatures_to_remove=\"none\"):")
        removed_sigs = paste0("\tSBS25 (dubious and specific to Hodgkin's lymphoma cell lines)\n",
                              "\tSBS31 (associated with platinum drug chemotherapy)\n",
                              "\tSBS32 (associated with azathioprine treatment)\n",
                              "\tSBS35 (associated with platinum drug chemotherapy)")
        message(crayon::black(removed_sigs))
      }
    } else if (signature_choice == "cosmic_v2") {
      data("signatures.cosmic", package = "deconstructSigs")
      signatures = signatures.cosmic
    }

  } else if("data.frame" %in% class(signature_choice)) {
      if(! identical(colnames(signatures_cosmic_May2019), colnames(signature_choice))) {
        stop("Your signatures data.frame is not properly formatted. See the \"signatures_cosmic_May2019\" object for a template.")
      }
      signatures = signature_choice
  } else {
    stop(bad_sig_choice_msg)
  }

  # make sure all signatures that user has requested removed are actually valid signatures
  if(length(signatures_to_remove) != 1 || tolower(signatures_to_remove[1]) != "none") {
    if(any(! signatures_to_remove %in% rownames(signatures))) {
      stop("One or more signatures specified in signatures_to_remove are not valid signatures")
    }
  }



  #TODO: it is possible that the only mutation in a tumor is a recurrent variant,
  # and if we remove this variant, then the tumor is lost from the mutation
  # rate calculation. Should this be added back in if algorithm_choice
  # is something involving averages? ---- UPDATE: it is added in "weightd", need to
  # add in other options
  #
  # TODO: really, averages should be averages of weight, and then product is calculated
  # individually per tumor, as opposed to average product and average rate.


  # pre-process ----
  message("Calculating trinucleotide composition and signatures...")
  MAF = cesa@maf # calculate weights
  sample_ID_column = "Unique_Patient_Identifier"
  chr_column = "Chromosome"
  pos_column = "Start_Position"
  ref_column = "Reference_Allele"
  alt_column = "Tumor_Allele"

  ds_maf = MAF
  if(remove_recurrent){
    # take just SNVs
    ds_maf <- MAF[which(MAF[, ref_column] %in% c('A', 'T', 'C', 'G') & MAF[, alt_column] %in% c('A', 'T', 'C', 'G')),]

    # remove all recurrent SNVs (SNVs appearing in more than one sample)
    duplicated_vec_first <- duplicated(ds_maf[,c(pos_column,chr_column,alt_column)])
    duplicated_vec_last <- duplicated(ds_maf[,c(pos_column,chr_column,alt_column)],fromLast=T)
    duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
    if (length(duplicated_vec_pos) > 0) {
      ds_maf <- ds_maf[-duplicated_vec_pos,]
    }
  }


  tumors_with_a_mutation_rate <- unique(ds_maf[,sample_ID_column])
  substitution_counts <- table(ds_maf[,sample_ID_column])
  tumors_with_50_or_more <- names(which(substitution_counts>=50))
  tumors_with_less_than_50 <- setdiff(ds_maf[,sample_ID_column],tumors_with_50_or_more)


  # this will move elsewhere once support for arbitrary genome build is finished
  data("tri.counts.genome", package = "deconstructSigs")
  data("tri.counts.exome", package = "deconstructSigs")

  trinuc_breakdown_per_tumor <- deconstructSigs::mut.to.sigs.input(mut.ref = ds_maf,
                                                                   sample.id = sample_ID_column,
                                                                   chr = chr_column,
                                                                   pos = pos_column,
                                                                   ref = ref_column,
                                                                   alt = alt_column,
                                                                   bsg = cesa@genome)


  deconstructSigs_trinuc_string <- colnames(trinuc_breakdown_per_tumor)



  trinuc_proportion_matrix <- matrix(data = NA,
                                     nrow = nrow(trinuc_breakdown_per_tumor),
                                     ncol = ncol(trinuc_breakdown_per_tumor))
  rownames(trinuc_proportion_matrix) <- rownames(trinuc_breakdown_per_tumor)
  colnames(trinuc_proportion_matrix) <- colnames(trinuc_breakdown_per_tumor)

  signatures_output_list <- vector(mode = "list",length = nrow(trinuc_breakdown_per_tumor))
  names(signatures_output_list) <- rownames(trinuc_breakdown_per_tumor)


  # algorithms ----

  if(signature_choice == "signatures_cosmic_May2019"){
    signatures <- signatures_cosmic_May2019 # part of CES data
  }

  if(signature_choice == "signatures.cosmic"){
    signatures <- data("signatures.cosmic", package = "deconstructSigs")
  }


  if(artifact_accounting && signature_choice == "cosmic_v3"){
    data("signatures_names_matrix", package="cancereffectsizeR")
    possible_artifact_signatures <- signatures_names_matrix[startsWith(x = signatures_names_matrix[,2], prefix = "*"),1]
    tumors_zeroed_out <- NULL
  }

  main_signatures <- signatures



  if(algorithm_choice == "weighted"){

    for(tumor_name in 1:nrow(trinuc_breakdown_per_tumor)){
      signatures <- main_signatures
      these_signatures_to_remove <- signatures_to_remove




      # Assignment rules for hypermutator tumors found
      # in "SigProfiler_signature_assignment_rules" supp data here:
      # https://doi.org/10.1101/322859
      # (these rules are for whole-exome, not whole-genome.
      # Different rules apply then, also in that supp data)
      if(v3_exome_hypermutation_rules) {
        if(substitution_counts[rownames(trinuc_breakdown_per_tumor)[tumor_name]] < 2*10^3){

          these_signatures_to_remove <- c(these_signatures_to_remove, "SBS10a","SBS10b")

          if(substitution_counts[rownames(trinuc_breakdown_per_tumor)[tumor_name]] < 2*10^2){

            these_signatures_to_remove <- c(these_signatures_to_remove,
                                            "SBS6",
                                            "SBS14",
                                            "SBS15",
                                            "SBS20",
                                            "SBS21",
                                            "SBS26",
                                            "SBS44")

          }
        }
      }

      if(!is.null(these_signatures_to_remove)){
        signatures <- signatures[-which(rownames(signatures) %in% these_signatures_to_remove),]
      }


      signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                            signatures.ref = signatures,
                                                            sample.id = rownames(trinuc_breakdown_per_tumor)[tumor_name],
                                                            contexts.needed = TRUE,
                                                            tri.counts.method = 'exome2genome')


      # add the removed signatures back in
      if(!is.null(these_signatures_to_remove)){
        df_to_add <- as.data.frame(matrix(data = 0,nrow = 1,ncol = length(these_signatures_to_remove)),stringsAsFactors=F)
        colnames(df_to_add) <- these_signatures_to_remove
        signatures_output$weights <- cbind(signatures_output$weights,df_to_add)
        # match original ordering
        signatures_output$weights <-
          signatures_output$weights[,match(rownames(main_signatures),colnames(signatures_output$weights))]

        signatures <- main_signatures
      }


      # some of the COSMIC v3 signatures are likely artifacts,
      # once accounted for in deconstructSigs::whichSignatures,
      # we remove artifacts and renormalize so that we can
      # determine mutational flux in tumor from true sources
      if(artifact_accounting){


        if(any(signatures_output$weights[possible_artifact_signatures] > 0 )){

          # often, sum of best estimation of signatures is < 1.
          # will renormalize to this original total explanation of weight.
          current_sum <- sum(signatures_output$weights)

          # remove artifacts from mutational signatures flux
          signatures_output$weights[possible_artifact_signatures] <- 0

          # if other weights exist, renormalize
          # else, assign this tumor the average product and weight of tumors
          # with > 50 variants.
          if(any(signatures_output$weights>0)){

            signatures_output$weights <- signatures_output$weights/(sum(signatures_output$weights)/current_sum)

            signatures_output$product <- as.matrix(signatures_output$weights) %*% as.matrix(signatures)

          }else{

            signatures_output$weights <- signatures_output$weights*0

            signatures_output$product <- as.matrix(signatures_output$weights) %*% as.matrix(signatures)

            # if zeroed out, remove from downstream analyses
            if(rownames(trinuc_breakdown_per_tumor)[tumor_name] %in% tumors_with_50_or_more){
              tumors_with_50_or_more <- tumors_with_50_or_more[-which(tumors_with_50_or_more == rownames(trinuc_breakdown_per_tumor)[tumor_name])]
            }

            if(rownames(trinuc_breakdown_per_tumor)[tumor_name] %in% tumors_with_less_than_50){
              tumors_with_less_than_50 <- tumors_with_less_than_50[-which(tumors_with_less_than_50 == rownames(trinuc_breakdown_per_tumor)[tumor_name])]
            }

            tumors_zeroed_out <- c(tumors_zeroed_out, rownames(trinuc_breakdown_per_tumor)[tumor_name])


          }




        }

      }


      signatures_output_list[[rownames(trinuc_breakdown_per_tumor)[tumor_name]]] <-
        list(signatures_output = signatures_output,
             substitution_count = length(which(MAF[,sample_ID_column] ==
                                                 rownames(trinuc_breakdown_per_tumor)[tumor_name])))


      # if all 0, it means accounting for artifacts zeroed out weights,
      # and we will later assume weights for this tumor consistent with
      # average weights of the tissue.
      if(any(signatures_output$product!=0)){

        trinuc_proportion_matrix[rownames(trinuc_breakdown_per_tumor)[tumor_name],] <- signatures_output$product/sum( signatures_output$product) #need it to sum to 1.

        # Not all trinuc weights in the cosmic dataset are nonzero for certain signatures
        # This leads to the rare occasion where a certain combination of signatues leads to ZERO
        # rate for particular trinucleotide contexts.
        # True rate is nonzero, as we do see those variants in those tumors, so renormalizing
        # the rates by adding the lowest nonzero rate to all the rates and renormalizing


        if(0 %in% trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),]){
          # finding the lowest nonzero rate
          lowest_rate <- min(trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),
                                                      -which(trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),]==0)])

          # adding it to the rates
          trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),] <- trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),] + lowest_rate

          # renormalizing to 1
          trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),] <-
            trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),] / sum(trinuc_proportion_matrix[which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),])
        }
      }else{
        # if zeroed out, remove, add back later.
        trinuc_proportion_matrix <- trinuc_proportion_matrix[-which(rownames(trinuc_proportion_matrix) == rownames(trinuc_breakdown_per_tumor)[tumor_name]),]
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


    averaged_product_mat <- as.data.frame(matrix(data = averaged_product,nrow=1),stringsAsFactors=F)
    rownames(averaged_product_mat) <- "averaged_product"; colnames(averaged_product_mat) <- names(averaged_product)

    averaged_weight_deconstructed <- deconstructSigs::whichSignatures(tumor.ref = averaged_product_mat,
                                                                      signatures.ref = signatures,
                                                                      sample.id = "averaged_product",
                                                                      contexts.needed = TRUE,
                                                                      tri.counts.method = 'exome2genome')

    averaged_weight <- averaged_weight_deconstructed$weights



    if(artifact_accounting){


      if(any(averaged_weight[possible_artifact_signatures] > 0 )){

        current_sum <- sum(averaged_weight)

        # remove artifacts from mutational signatures flux
        averaged_weight[possible_artifact_signatures] <- 0

        # if other weights exist, renormalize
        # else, assign this tumor the average product and weight of tumors
        # with > 50 variants.
        if(any(signatures_output$weights>0)){

          averaged_weight <- averaged_weight/(sum(averaged_weight)/current_sum)

        }else{
          stop("After accounting for artifacts,
              all weights are zero in the average weights within tumors with
              >50 substitutions")
        }
      }
    }


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
          matrix(data=trinuc_proportion_matrix[tumors_with_less_than_50[tumor_name_index],], nrow=1)
        rownames(signatures_output_list[[tumors_with_less_than_50[tumor_name_index]]]$signatures_output$product) <-
          tumors_with_less_than_50[tumor_name_index]
        colnames(signatures_output_list[[tumors_with_less_than_50[tumor_name_index]]]$signatures_output$product) <-
          colnames(signatures_output_list[[tumors_with_less_than_50[tumor_name_index]]]$signatures_output$tumor)
      }
    }


    # If a tumor had all variants removed in preprocessing
    # (meaning it only had recurrent variants)
    # then add back in the mutation rates as the average of the > 50

    run_next_block <- F
    if((length(signatures_output_list) < length(unique(MAF[,sample_ID_column])))){
      run_next_block <- T
    }
    if(artifact_accounting & signature_choice == "signatures_cosmic_May2019" & !is.null(tumors_zeroed_out)){
      run_next_block <- T
    }

    if(run_next_block){

      if(length(signatures_output_list) < length(unique(MAF[,sample_ID_column])) ){
        tumors_to_add <- unique(MAF[,sample_ID_column])[which(!unique(MAF[,sample_ID_column]) %in% names(signatures_output_list))]
      }else{
        tumors_to_add <- NULL
      }

      if(artifact_accounting & signature_choice == "signatures_cosmic_May2019"){
        if(!is.null(tumors_zeroed_out)){

          tumors_to_add <- c(tumors_to_add, tumors_zeroed_out)

        }

      }

      matrix_to_add <- matrix(data = NA, nrow = length(tumors_to_add), ncol = ncol(trinuc_breakdown_per_tumor))
      rownames(matrix_to_add) <- tumors_to_add
      colnames(matrix_to_add) <- colnames(trinuc_breakdown_per_tumor)

      for(matrix_row in 1:nrow(matrix_to_add)){

        matrix_to_add[matrix_row,] <- as.numeric(averaged_product)
        signatures_output_list[[rownames(matrix_to_add)[matrix_row]]]$signatures_output$product <- as.data.frame(matrix(data=averaged_product,nrow=1),stringsAsFactors=F)
        rownames(signatures_output_list[[rownames(matrix_to_add)[matrix_row]]]$signatures_output$product) <- rownames(matrix_to_add)[matrix_row]
        colnames(signatures_output_list[[rownames(matrix_to_add)[matrix_row]]]$signatures_output$product) <- names(averaged_product)
        signatures_output_list[[rownames(matrix_to_add)[matrix_row]]]$signatures_output$weights <- as.data.frame(matrix(data=as.numeric(averaged_weight), nrow=1),stringsAsFactors=F)
        rownames(signatures_output_list[[rownames(matrix_to_add)[matrix_row]]]$signatures_output$weights) <- rownames(matrix_to_add)[matrix_row]
        colnames(signatures_output_list[[rownames(matrix_to_add)[matrix_row]]]$signatures_output$weights) <- names(averaged_weight)
      }
      trinuc_proportion_matrix <- rbind(trinuc_proportion_matrix, matrix_to_add)
      tumors_with_a_mutation_rate <- rownames(trinuc_proportion_matrix)
    }
  }


  # takes the average of all profiles with >= 50 substitutions
  if(algorithm_choice == "all_average"){

    for(tumor_name in 1:nrow(trinuc_breakdown_per_tumor)){

      signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                            signatures.ref = signatures,
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
                                                            signatures.ref = signatures,
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

    # now, trinuc_proportion_matrix has the proportion of all trinucs in every tumor.

  }



  if(algorithm_choice == "all_calculated"){

    for(tumor_name in 1:nrow(trinuc_breakdown_per_tumor)){

      signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                            signatures.ref = signatures,
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


  # Need to make sure we are only calculating the selection intensities from
  # tumors in which we are able to calculate a mutation rate

  cesa@maf = MAFdf(MAF[MAF$"Unique_Patient_Identifier" %in% tumors_with_a_mutation_rate,])
  no_mutation_rate = MAF[! MAF$"Unique_Patient_Identifier" %in% tumors_with_a_mutation_rate,]
  if (nrow(no_mutation_rate) > 0) {
    no_mutation_rate$Exclusion_Reason = "no_tumor_mutation_rate"
    cesa@excluded = MAFdf(rbind(cesa@excluded, no_mutation_rate))  
    message(paste("Note: Some tumor(s) have only recurrent mutations (mutations also appearing in other samples) ",
                  "so a baseline mutation rate cannot be calculated. MAF data for these have been excluded from further analysis."))
  }
  cesa@trinucleotide_mutation_weights = list(tumors_with_a_mutation_rate=tumors_with_a_mutation_rate,
                                             trinuc_proportion_matrix=trinuc_proportion_matrix,
                                             signatures_output_list=signatures_output_list,
                                             algorithm_choice=algorithm_choice)
  return(cesa)
}
