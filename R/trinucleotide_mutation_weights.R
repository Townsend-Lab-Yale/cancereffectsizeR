#' Trinucleotide mutation weights
#'
#' Calculates expected relative rates of SNV substitutions by trinucleotide context in tumors
#' using deconstructSigs to estimate weightings of known mutational signatures (e.g., COSMIC v3)
#'
#' The purpose of this function is to calculate relative rates of trinucleotide-context-specific SNV
#' mutations within tumors that can be attributed to known tumor-specific mutational processes. Since
#' recurrent mutations (those that appear in multiple tumor samples) are more likely to be undergoing
#' selection than other mutations, by default only non-recurrent mutations are considered in order to
#' get a closer estimation of mutation rates independent of selection. This function currently uses
#' the deconstructSigs to assign mutational signature weightings to tumors.
#' As deconstructSigs suggests that a tumor must have at least 50 mutations for signature attribution
#' to be reliable, this function offers several methods for dealing with tumors with fewer mutations.
#' The default method, "weighted".... (to be continued)
#'
#' @param cesa # CESAnalysis object
#' @param cores how many cores to use to process tumor samples in parallel (requires parallel package)
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
#'        (if using cosmic_v3, just leave option default instead of passing your own data frame, or you'll get improper behavior)
#' @param v3_artifact_accounting when COSMIC v3 signatures associated with sequencing artifacts are detected, renormalizes to isolate "true" sources of mutational flux.
#' @param signatures_to_remove specify any signatures to exclude from analysis; some signatures automatically get excluded
#'     from COSMIC v3 analyses; set to signatures_to_remove="none" to prevent this behavior (and )
#' @param v3_exome_hypermutation_rules T/F on whether to follow the mutation count rules outlined in https://doi.org/10.1101/322859, the manuscript reported the v3 COSMIC signature set.
#'
#' @export
#'
#'
#'
trinucleotide_mutation_weights <- function(cesa,
                                           cores = 1,
                                           algorithm_choice = "weighted",
                                           remove_recurrent = TRUE,
                                           signature_choice = "cosmic_v3",
                                           v3_artifact_accounting = TRUE,
                                           v3_exome_hypermutation_rules = TRUE,
                                           signatures_to_remove = "" # cosmic_v3 analysis gets signatures added here later unless "none"
                                           ){

  algorithms = c("weighted", "all_average", "nearest_neighbor", "all_calculated")
  if(! "character" %in% class(algorithm_choice) || ! algorithm_choice %in% algorithms ) {
    stop(paste0("algorithm_choice (for dealing with samples with fewer than 50 variants) must be one of ", paste(algorithms, collapse = ", ")))
  }
  if(! "character" %in% class(signatures_to_remove)) {
    stop("signatures_to_remove should be a character vector")
  }
  bad_sig_choice_msg = "signature_choice should be \"cosmic_v3\", \"cosmic_v2\", or a properly-formatted signatures data.frame"
  running_cosmic_v3 = FALSE # gets set to true when user chooses
  if("character" %in% class(signature_choice)) {
    signature_choice = tolower(signature_choice)
    if (length(signature_choice) != 1 || ! signature_choice %in% c("cosmic_v3", "cosmic_v2")) {
      stop(bad_sig_choice_msg)
    }
    if(signature_choice == "cosmic_v3") {
      running_cosmic_v3 = TRUE
      signatures = signatures_cosmic_May2019
      message("Analyzing tumors for COSMIC v3 (May 2019) SNV mutational signatures....")
      if (v3_artifact_accounting) {
        message(crayon::black("Sequencing artifact signatures will be subtracted to isolate true mutational processes (disable with v3_artifact_accounting=FALSE)."))
        v3_artifact_accounting = TRUE
      }
      if (v3_exome_hypermutation_rules) {
        message(crayon::black("Samples with many variants will have special exome hypermutation rules applied (disable with v3_exome_hypermutation_rules=FALSE)."))
      }
      if(length(signatures_to_remove) == 1 && signatures_to_remove == "") {
        signatures_to_remove = c("SBS25","SBS31","SBS32","SBS35")
        message(crayon::black("The following signatures will be excluded (to include all signatures, set signatures_to_remove=\"none\"):"))
        removed_sigs = paste0("\tSBS25 (dubious and specific to Hodgkin's lymphoma cell lines)\n",
                              "\tSBS31 (associated with platinum drug chemotherapy)\n",
                              "\tSBS32 (associated with azathioprine treatment)\n",
                              "\tSBS35 (associated with platinum drug chemotherapy)")
        message(crayon::black(removed_sigs))
      }
    } else if (signature_choice == "cosmic_v2") {
      data("signatures.cosmic", package = "deconstructSigs")
      signatures = signatures.cosmic
    } else {
        stop(bad_sig_choice_msg)
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
  if(length(signatures_to_remove) == 1 && tolower(signatures_to_remove[1]) == "none") {
    signatures_to_remove = character()
  }

  # can't apply our COSMIC v3 enhancements unless we're using those signatures, naturally
  if(! running_cosmic_v3) {
    v3_exome_hypermutation_rules = FALSE
    v3_artifact_accounting = FALSE
  }

  #
  # TODO: really, averages should be averages of weight, and then product is calculated
  # individually per tumor, as opposed to average product and average rate.


  MAF = cesa@maf
  sample_ID_column = "Unique_Patient_Identifier"
  chr_column = "Chromosome"
  pos_column = "Start_Position"
  ref_column = "Reference_Allele"
  alt_column = "Tumor_Allele"

  ds_maf = MAF
  if(remove_recurrent){
    # take just SNVs
    bases = c('A', 'T', 'C', 'G')
    ds_maf <- MAF[Reference_Allele %in% bases & Tumor_Allele %in% bases]

    # remove all recurrent SNVs (SNVs appearing in more than one sample)
    duplicated_vec_first <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)])
    duplicated_vec_last <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)],fromLast=T)
    duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
    if (length(duplicated_vec_pos) > 0) {
      ds_maf <- ds_maf[-duplicated_vec_pos,]
    }
  }

  all_tumors = unique(MAF$Unique_Patient_Identifier)
  tumors_with_a_mutation_rate <- unique(ds_maf$Unique_Patient_Identifier)
  tumors_with_only_recurrent_mutations = setdiff(all_tumors, tumors_with_a_mutation_rate)
  substitution_counts = table(ds_maf$Unique_Patient_Identifier)
  tumors_with_50_or_more = names(which(substitution_counts>=50))
  tumors_with_less_than_50 = setdiff(tumors_with_a_mutation_rate, tumors_with_50_or_more)

  if (any(substitution_counts < 50)) {
    message(crayon::black(paste0("Note: ", length(tumors_with_less_than_50), " of ", length(tumors_with_a_mutation_rate),
      " tumor samples have fewer than 50 mutations. Mutational signature weightings for these samples ",
      "will be determined using the \"", algorithm_choice, "\" method (see documentation).")))
  }


  # this will move elsewhere once support for arbitrary genome build is finished
  data("tri.counts.genome", package = "deconstructSigs")
  data("tri.counts.exome", package = "deconstructSigs")

  withCallingHandlers(
  {
    trinuc_breakdown_per_tumor = deconstructSigs::mut.to.sigs.input(mut.ref = ds_maf,
                                                                 sample.id = sample_ID_column,
                                                                 chr = chr_column,
                                                                 pos = pos_column,
                                                                 ref = ref_column,
                                                                 alt = alt_column,
                                                                 bsg = cesa@genome)
  }, warning = function(w) {
    # user was already notified about samples with <50 mutations above
    if (startsWith(conditionMessage(w), "Some samples have fewer than 50 mutations"))
      invokeRestart("muffleWarning")
  })


  # algorithms ----
  if(signature_choice == "signatures_cosmic_May2019"){
    signatures <- signatures_cosmic_May2019 # part of CES data
  }

  if(signature_choice == "signatures.cosmic"){
    data("signatures.cosmic", package = "deconstructSigs")
    signatures <- signatures.cosmic # v2 signatures from deconstructSigs
  }

  message("Calculating mutational signature weightings...")

  artifact_signatures = NULL
  if (v3_artifact_accounting) {
    # COSMIC v3 artifact signatures are read in from package data 
    artifact_signatures = signatures_names_matrix[startsWith(x = signatures_names_matrix[,2], prefix = "*"),1] 
    if(any(artifact_signatures %in% signatures_to_remove)) {
      warning(paste0("Warning: You are have chosen to remove at least one sequencing-artifact-associated signature from analysis,
                      which will change how artifact accounting (which has been left on) behaves."))
    }
  }

  # hypermutation signatures from https://doi.org/10.1101/322859
  cosmic_v3_highly_hm_sigs = c("SBS10a","SBS10b")
  cosmic_v3_modest_hm_sigs = c("SBS6","SBS14","SBS15","SBS20","SBS21","SBS26","SBS44")



  # for parallelization, a function to process each tumor
  process_tumor = function(tumor_name) {
    tumor_trinuc_counts = as.data.frame(trinuc_breakdown_per_tumor[tumor_name,]) # deconstructSigs requires a data.frame
    num_variants = sum(tumor_trinuc_counts)

    # Apply exome hypermuation rules if using
    ### If tumor is evidently not hypermutated, remove hypermutation-associated signatures from consideration
    ### Note: Assignment rules for hypermutator tumors found in
    ### "SigProfiler_signature_assignment_rules" supp data here: https://doi.org/10.1101/322859
    ### (These rules are for whole-exome, not whole-genome. Different rules apply then, also in that supp data.)
    current_sigs_to_remove = signatures_to_remove
    
    # To-do: for WGS data, thresholds are 10^5 and 10^4, respectively
    if (v3_exome_hypermutation_rules) {
      if(num_variants < 2000) {
        current_sigs_to_remove = union(current_sigs_to_remove, cosmic_v3_highly_hm_sigs)
      }
      if(num_variants < 200) {
        current_sigs_to_remove = union(current_sigs_to_remove, cosmic_v3_modest_hm_sigs)
      }
    }
    ds = cancereffectsizeR:::run_deconstructSigs(tumor_trinuc_counts = tumor_trinuc_counts, signatures_df = signatures,
                                                 signatures_to_remove = current_sigs_to_remove, artifact_signatures = artifact_signatures)
    return(list(list(signatures_output = ds[[1]], substitution_count = num_variants), ds[[2]]))
  }

  tumor_names = rownames(trinuc_breakdown_per_tumor)
  ds_output = pbapply::pblapply(tumor_names, process_tumor, cl = cores)

  # store results
  # matrix with rows = tumors, columns = relative rates of mutation for each trinuc SNV (starts empty)
  trinuc_proportion_matrix <- matrix(data = NA, nrow = nrow(trinuc_breakdown_per_tumor), ncol = ncol(trinuc_breakdown_per_tumor),
                                     dimnames = list(rownames(trinuc_breakdown_per_tumor), colnames(trinuc_breakdown_per_tumor)))
  signatures_output_list = list()
  zeroed_out_tumors = character() # for tumors that get all zero signature weights (rare)
  for(i in 1:length(ds_output)) {
    tumor_name = tumor_names[i]
    ds = ds_output[[i]]
    if(all(ds[[1]]$signatures_output$product == 0)) {
      zeroed_out_tumors = c(zeroed_out_tumors, tumor_name)
    }
    signatures_output_list[[tumor_name]] = ds[[1]]
    trinuc_proportion_matrix[i,] = ds[[2]]
  }


  # handle rare occurrence of zeroed-out tumors
  tumors_with_50_or_more = setdiff(tumors_with_50_or_more, zeroed_out_tumors)
  tumors_with_less_than_50 = setdiff(tumors_with_less_than_50, zeroed_out_tumors)
  tumors_with_a_mutation_rate = setdiff(tumors_with_a_mutation_rate, zeroed_out_tumors)


  # If any samples have <50 mutations, determine weightings for those samples using method specified by user
  # Note there's nothing more to do if choice was "all_calculated"
  if(length(tumors_with_less_than_50) > 0) {
    # nearest_neighbor method falls back to using a 50+ average weighting for zeroed-out tumors and those with
    # only recurrent mutations
    nearest_neighbor_needs_mean = FALSE
    if(algorithm_choice == "nearest_neighbor") {
      if(length(zeroed_out_tumors) > 0 || length(tumors_with_only_recurrent_mutations) > 0) {
        nearest_neighbor_needs_mean = TRUE
      }
    }
    if (algorithm_choice %in% c("weighted", "all_average") || nearest_neighbor_needs_mean) {
      mean_trinuc_prop = colMeans(trinuc_proportion_matrix[tumors_with_50_or_more, , drop = F])

      # convert to data frame because that's what deconstructSigs wants
      mean_trinuc_prop = as.data.frame(t(mean_trinuc_prop))
      
      
      # Hypermutation signatures are presumed absent from tumors with fewer than 50 variants, but may be present in the 50+ mutation tumors.
      # Therefore, they will be treated like artifact signatures: included in deconstructSigs signature weight calculation, but normalized
      # out when calculating the "true" relative trinucleotide SNV mutation rates
      # however, they are left in for the "all_average" method in the spirit of assuming all tumors have same mutational processes
      mean_calc_artifact_signatures = artifact_signatures
      if (v3_exome_hypermutation_rules && algorithm_choice == "weighted") {
        mean_calc_artifact_signatures = unique(c(artifact_signatures, cosmic_v3_modest_hm_sigs, cosmic_v3_highly_hm_sigs))
      }
      
      mean_ds <- cancereffectsizeR:::run_deconstructSigs(tumor_trinuc_counts = mean_trinuc_prop, signatures_df = signatures, 
                                                        signatures_to_remove = signatures_to_remove,
                                                        artifact_signatures = mean_calc_artifact_signatures)[[1]] # just need signatures_output element
      mean_weights <- mean_ds$weights
      # this should never happen
      if(all(mean_weights == 0)) {
        stop("Somehow, mean signature weights across all tumors with >50 mutations are all zero.")
      }
      if (algorithm_choice == "weighted") {
         for (tumor in tumors_with_less_than_50) {
          own_weighting = substitution_counts[tumor] / 50
          group_weighting = 1 - own_weighting
          trinuc_proportion_matrix[tumor, ] = as.numeric(trinuc_proportion_matrix[tumor, ] * own_weighting + mean_trinuc_prop * group_weighting)
          new_ds = signatures_output_list[[tumor]]$signatures_output
          new_ds$weights = new_ds$weights * own_weighting + mean_weights * group_weighting
          new_ds$product = trinuc_proportion_matrix[tumor, , drop=F] 
          signatures_output_list[[tumor]]$signatures_output = new_ds
        }       
      } else if (algorithm_choice == "all_average") {
        mean_trinuc_prop = as.numeric(mean_trinuc_prop)
        for(i in 1:nrow(trinuc_proportion_matrix)) {
          trinuc_proportion_matrix[i,] = mean_trinuc_prop
        }
      }

    }
    if (algorithm_choice == "nearest_neighbor") {
      distance_matrix <- as.matrix(dist(trinuc_breakdown_per_tumor))
      #find closest tumor that has at least 50 mutations
      for(tumor in tumors_with_less_than_50) {
        tumor_dist  = distance_matrix[tumor, tumors_with_50_or_more, drop = F]
        closest_tumor = colnames(tumor_dist)[which(tumor_dist == sort(tumor_dist)[1])] # get name of closest tumor
        trinuc_proportion_matrix[tumor,] = trinuc_proportion_matrix[closest_tumor,]
        signatures_output_list[[tumor]]$signatures_output = signatures_output_list[[closest_tumor]]$signatures_output
        signatures_output_list[[tumor]]$tumor_signatures_used = closest_tumor
      }
    }
  }

  # Tumors with no recurrent mutations and  tumors with zeroed-out signatures both get assigned the 50+ average weightings
  # Exception: if user chose "all_calculated", these samples will get thrown out
  tumors_to_add = c(tumors_with_only_recurrent_mutations, zeroed_out_tumors)
  for (tumor in tumors_to_add) {
    # quirky handling to match exact output of previous version of script 
    if (tumor %in% names(signatures_output_list)) {
      signatures_output_list[[tumor]] = NULL
      trinuc_proportion_matrix = trinuc_proportion_matrix[rownames(trinuc_proportion_matrix) != tumor,]
    }
  }
  if (algorithm_choice != "all_calculated" && length(tumors_to_add) > 0) {
    matrix_to_add <- matrix(data = NA, nrow = length(tumors_to_add), ncol = ncol(trinuc_breakdown_per_tumor),
                            dimnames = list(tumors_to_add, colnames(trinuc_breakdown_per_tumor)))
    mean_trinuc_prop = as.numeric(mean_trinuc_prop)
    for(tumor in tumors_to_add) {
      matrix_to_add[tumor,] = mean_trinuc_prop
      signatures_output_list[[tumor]]$signatures_output$product = mean_trinuc_prop
      signatures_output_list[[tumor]]$signatures_output$weights = mean_weights

    }
    trinuc_proportion_matrix <- rbind(trinuc_proportion_matrix, matrix_to_add)
    tumors_with_a_mutation_rate = c(tumors_with_a_mutation_rate, tumors_to_add)
  }

  # Need to make sure we are only calculating the selection intensities from
  # tumors in which we are able to calculate a mutation rate
  cesa@maf = MAF[MAF$"Unique_Patient_Identifier" %in% tumors_with_a_mutation_rate,]
  no_mutation_rate = MAF[! MAF$"Unique_Patient_Identifier" %in% tumors_with_a_mutation_rate,]
  if (nrow(no_mutation_rate) > 0) {
    no_mutation_rate$Exclusion_Reason = "no_tumor_mutation_rate"
    cesa@excluded = rbind(cesa@excluded, no_mutation_rate) 
    message(paste("Note: Some tumor(s) have only recurrent mutations (mutations also appearing in other samples) ",
                  "so a baseline mutation rate cannot be calculated. MAF data for these have been excluded from further analysis."))
  }
  cesa@trinucleotide_mutation_weights = list(tumors_with_a_mutation_rate=tumors_with_a_mutation_rate,
                                             trinuc_proportion_matrix=trinuc_proportion_matrix,
                                             signatures_output_list=signatures_output_list,
                                             algorithm_choice=algorithm_choice)
  return(cesa)
}


