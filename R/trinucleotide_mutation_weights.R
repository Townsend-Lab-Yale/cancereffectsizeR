#' Trinucleotide mutation weights
#'
#'
#' This function calculates relative rates of trinucleotide-context-specific SNV
#' mutations within tumors by attributing SNVs to mutational processes represented
#' in mutation signature sets (such as COSMIC v3). This function currently uses
#' deconstructSigs to assign mutational signature weightings to tumors. Tumor samples
#' with targeted sequencing data are assigned the average trinucleotide
#' mutation rates calculated across all exome/genome data, which means that you need at
#' least some exome or genome data to run.
#' 
#'
#' To reduce the influence of selection on the estimation of relative trinucleotide mutation
#' rates, this function only uses non-recurrent SNVs (i.e., those that do not appear in more than one 
#' sample in the data set).
#'
#' @param cesa # CESAnalysis object
#' @param cores how many cores to use to process tumor samples in parallel (requires parallel package)
#' @param signature_choice "cosmic_v3" (default), "cosmic_v2", or a properly-formatted data frame with of trinucleotide signatures 
#'        (if using cosmic_v3, just leave option default instead of passing your own data frame, or you'll get improper behavior)
#' @param assume_identical_mutational_processes (default FALSE) instead of assigning different signature weights to each tumor (reflective of
#'   tumor-specific mutational processes), use well-mutated tumors (those with number of eligible mutations meeting sig_averaging_threshold)
#'   calculate group average signature weights and assign these to all tumors
#' @param sig_averaging_threshold Mutational threshold (default 50) that determines which tumors inform the
#'   calculation of group-average signature weightings. When assume_identical_mutational_processes == FALSE (the default), 
#'   these group averages are blended into the signature weightings of tumors with few mutations (those below the threshold).
#' @param v3_artifact_accounting when COSMIC v3 signatures associated with sequencing artifacts are detected, renormalizes to isolate "true" sources of mutational flux.
#' @param signatures_to_remove specify any signatures to exclude from analysis; some signatures automatically get excluded
#'     from COSMIC v3 analyses; set signatures_to_remove = "none" to prevent this behavior
#' @param v3_hypermutation_rules T/F on whether to follow the mutation count rules outlined in https://doi.org/10.1101/322859, the manuscript reported the v3 COSMIC signature set.
#'
#' @export
#'
#'
#'
trinucleotide_mutation_weights <- function(cesa,
                                           cores = 1,
                                           signature_choice = "cosmic_v3",
                                           assume_identical_mutational_processes = FALSE,
                                           sig_averaging_threshold = 50,
                                           v3_artifact_accounting = TRUE,
                                           v3_hypermutation_rules = TRUE,
                                           signatures_to_remove = "" # cosmic_v3 analysis gets signatures added here later unless "none"
                                           ){

  if(is.null(cesa) || ! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be a CESAnalysis object")
  }
  if(cesa@maf[, .N] == 0) {
    stop("No MAF data in the CESAnalysis")
  }
  
  if(all(cesa@samples$coverage == "targeted")) {
    stop("We can't estimate relative trinucleotide mutation rates without some exome/genome data in the CESAnalysis (all data is targeted sequencing).")
  }
  
  if(! is.logical(assume_identical_mutational_processes) || length(assume_identical_mutational_processes) != 1 || is.na(assume_identical_mutational_processes)) {
    stop("Expected assume_identical_mutational_processes to be TRUE/FALSE (default is FALSE).")
  }
  
  if(! is.numeric(sig_averaging_threshold) || length(sig_averaging_threshold) != 1 || is.na(sig_averaging_threshold) || sig_averaging_threshold < 0) {
    stop("Expected positive integer for sig_averaging_threshold")
  }
  sig_averaging_threshold = as.integer(sig_averaging_threshold) # below, we test that at least one tumor meets this threshold
  
  if(! "character" %in% class(signatures_to_remove)) {
    stop("signatures_to_remove should be a character vector")
  }
  bad_sig_choice_msg = "signature_choice should be \"cosmic_v3\", \"cosmic_v2\", or a properly-formatted signatures data.frame"
  running_cosmic_v3 = FALSE # gets set to true when user chooses
  signature_set_name = "custom" # gets set appropriately below for cosmic v2/v3
  if("character" %in% class(signature_choice)) {
    if(! GenomeInfoDb::providerVersion(cesa@genome) %in% c("hg19", "hg38")) {
      stop("When not running with the human genome (hg38 or hg19), signatures must be supplied as a data frame (see docs).")
    }
    signature_choice = tolower(signature_choice)
    if (length(signature_choice) != 1 || ! signature_choice %in% c("cosmic_v3", "cosmic_v2")) {
      stop(bad_sig_choice_msg)
    }
    if(signature_choice == "cosmic_v3") {
      signature_set_name = "COSMIC v3"
      running_cosmic_v3 = TRUE
      signatures = signatures_cosmic_May2019
      message("Analyzing tumors for COSMIC v3 (May 2019) SNV mutational signatures....")
      if (v3_artifact_accounting) {
        message(crayon::black("Sequencing artifact signatures will be subtracted to isolate true mutational processes (disable with v3_artifact_accounting=FALSE)."))
        v3_artifact_accounting = TRUE
      }
      if (v3_hypermutation_rules) {
        message(crayon::black("Samples with many variants will have special exome hypermutation rules applied (disable with v3_hypermutation_rules=FALSE)."))
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
      signature_set_name = "COSMIC v2"
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
    v3_hypermutation_rules = FALSE
    v3_artifact_accounting = FALSE
  }

  # keeping TGS data in the ds_maf until after recurrency testing
  ds_maf = cesa@maf[Variant_Type == "SNV"]

  # remove all recurrent SNVs (SNVs appearing in more than one sample)
  duplicated_vec_first <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)])
  duplicated_vec_last <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)],fromLast=T)
  duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
  if (length(duplicated_vec_pos) > 0) {
    ds_maf <- ds_maf[-duplicated_vec_pos,]
  }


  all_tumors = cesa@samples$Unique_Patient_Identifier # may include some tumors with no SNVs, or with only recurrent SNVs
  
  # can only find signatures in tumors with exome/genome data that have non-recurrent SNVs
  tumors_eligible_for_trinuc_calc = intersect(cesa@samples[coverage != "targeted", Unique_Patient_Identifier], unique(ds_maf$Unique_Patient_Identifier))
  tumors_needing_group_average_rates = setdiff(all_tumors, tumors_eligible_for_trinuc_calc)
  
  # subset to just the tumors that will be run through deconstructSigs
  ds_maf = ds_maf[Unique_Patient_Identifier %in% tumors_eligible_for_trinuc_calc]
  
  substitution_counts = table(ds_maf$Unique_Patient_Identifier)
  
  
  tumors_above_threshold = names(which(substitution_counts>= sig_averaging_threshold))
  if(length(tumors_above_threshold) == 0) {
    stop(paste0("No tumors have enough mutations to confidently assess mutational signatures. To run anyway, lower\n",
                "sig_averaging_threshold (possibly to 0), and consider setting assume_identical_mutational_processes = TRUE"))
  }
  # identify tumors that will have average weights blended into their weights (normally)
  # when assume_identical_mutational_processes = TRUE, data from these tumors isn't used at all
  tumors_below_threshold = names(which(substitution_counts < sig_averaging_threshold))


  tri.counts.genome = get_genome_data(cesa, "tri.counts.genome")
  
  # for each exome coverage gr (besides default generic, which is pre-calculated), tabulate trinucs
  exome_counts_by_gr = list()
  exomes_to_calc = cesa@samples[coverage == "exome", unique(covered_regions)]
  for (exome_name in exomes_to_calc) {
    if (exome_name %in% c("exome", "exome+")) {
      exome_counts_by_gr[[exome_name]] = get_genome_data(cesa, "tri.counts.exome")
    } else {
      exome_seq = getSeq(cesa@genome, cesa@coverage[[exome_name]])
      exome_tri_contexts = Biostrings::trinucleotideFrequency(exome_seq)
      exome_tri_contexts = colSums(exome_tri_contexts)
      
      # deconstructSigs_trinuc_string is internal in cancereffectsizeR
      # here, we need unique trinucleotide contexts (without mutations) in deconstructSigs ordering,
      # which is why we produce this by converting from their column headings
      context_names = unique(sub("\\[([ACTG]).*\\]", "\\1", deconstructSigs_trinuc_string))
      context_names = sort(context_names)
      reverse_complement_names = unique(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(context_names))))
      
      # reorder the counts as desired, then save as a data frame since that's what deconstructSigs wants
      exome_counts = exome_tri_contexts[context_names] + exome_tri_contexts[reverse_complement_names]
      exome_counts = data.frame(x = exome_counts)
      exome_counts_by_gr[[exome_name]] = exome_counts
    }
  }
  
  
  # build the data.frame required by deconstructSigs (and probably similar to what is required by most other SNV signature software)
  # rows are samples, columns are counts of each of 96 trinuc-context-specific mutations in the order expected by deconstructSigs
  trinuc = BSgenome::getSeq(cesa@genome, ds_maf$Chromosome, ds_maf$Start_Position - 1, ds_maf$Start_Position + 1, as.character = T)
  
  # internal dict converts trinuc/mut (e.g., GTA:C) into deconstructSigs format ("G[T>C]A")
  ds_muts = factor(trinuc_translator[paste0(trinuc, ":", ds_maf$Tumor_Allele), "deconstructSigs_format"], levels = deconstructSigs_trinuc_string)
  
  # mysteriously convert two-way table to data frame
  tmp = table(ds_maf$Unique_Patient_Identifier, ds_muts)
  counts = apply(tmp, 2, rbind)
  rownames(counts) = rownames(tmp)
  trinuc_breakdown_per_tumor = as.data.frame(counts)

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

    # Apply hypermuation rules if using
    ### If tumor is evidently not hypermutated, remove hypermutation-associated signatures from consideration
    ### Note: Assignment rules for hypermutator tumors found in
    ### "SigProfiler_signature_assignment_rules" supp data here: https://doi.org/10.1101/322859
    current_sigs_to_remove = signatures_to_remove
    
    # To-do: for WGS data, thresholds are 10^5 and 10^4, respectively
    if (v3_hypermutation_rules) {
      # apply hypermuation rules
      if (cesa@samples[tumor_name, coverage] == "exome") {
        if(num_variants < 2000) {
          current_sigs_to_remove = union(current_sigs_to_remove, cosmic_v3_highly_hm_sigs)
        }
        if(num_variants < 200) {
          current_sigs_to_remove = union(current_sigs_to_remove, cosmic_v3_modest_hm_sigs)
        }
      } else {
        # if not exome, must be genome, since TGS data not used in this function
        if(num_variants < 10^5) {
          current_sigs_to_remove = union(current_sigs_to_remove, cosmic_v3_highly_hm_sigs)
        }
        if(num_variants < 10^4) {
          current_sigs_to_remove = union(current_sigs_to_remove, cosmic_v3_modest_hm_sigs)
        }
      }
    }
    if (cesa@samples[tumor_name, coverage] == "genome") {
      normalization = "default" # this actually means no normalization (since signatures and MAF coverage are both whole-genome)
    } else {
      normalization = tri.counts.genome / exome_counts_by_gr[[cesa@samples[tumor_name, covered_regions]]]
    }
    
    return(cancereffectsizeR::run_deconstructSigs(tumor_trinuc_counts = tumor_trinuc_counts, tri.counts.method = normalization,
                                                signatures_df = signatures, signatures_to_remove = current_sigs_to_remove, artifact_signatures = artifact_signatures))
  }

  ds_output = pbapply::pblapply(tumors_eligible_for_trinuc_calc, process_tumor, cl = cores)
  
  # store results
  # matrix with rows = tumors, columns = relative rates of mutation for each trinuc SNV (starts empty)
  trinuc_proportion_matrix <- matrix(data = NA, nrow = length(all_tumors), ncol = ncol(trinuc_breakdown_per_tumor),
                                     dimnames = list(all_tumors, colnames(trinuc_breakdown_per_tumor)))
  signatures_output_list = list()
  zeroed_out_tumors = character() # for tumors that get all zero signature weights (rare)
  
  
  for(i in 1:length(ds_output)) {
    tumor_name = tumors_eligible_for_trinuc_calc[i]
    signatures_output = ds_output[[i]]
    signatures_output_list[[tumor_name]] = signatures_output
    
    # if all signature weights are zero (very unlikely, but would probably be due to artifact accounting),
    # then nothing to do but return the output and let user or package deal with it downstream
    if (all(signatures_output$product == 0)) {
      zeroed_out_tumors = c(zeroed_out_tumors, tumor_name)
    } else {
      # Some trinuc SNVs have substitution rates of zero under certain signatures.
      # In rare cases, a tumor's fitted combination of signatures can therefore also
      # have a substitution rate of zero for particular trinucleotide contexts.
      # If this happens, we add the lowest nonzero rate to all rates and renormalize.
      trinuc_prop = signatures_output$product/sum(signatures_output$product)
      if(any(trinuc_prop == 0)) {
        lowest_nonzero_rate = min(trinuc_prop[trinuc_prop != 0])
        trinuc_prop = trinuc_prop + lowest_nonzero_rate
        # renormalize so rates sum to 1
        trinuc_prop = trinuc_prop / sum(trinuc_prop)
      }
      trinuc_proportion_matrix[tumor_name,] = trinuc_prop
    }
  }

  # in the (very rare) occurrence of zeroed-out tumors, set them aside to be assigned group average signature weightings
  tumors_above_threshold = setdiff(tumors_above_threshold, zeroed_out_tumors)
  tumors_below_threshold = setdiff(tumors_below_threshold, zeroed_out_tumors)
  tumors_needing_group_average_rates = union(tumors_needing_group_average_rates, zeroed_out_tumors)


  # If any samples have <50 mutations, determine weightings for those samples using method specified by user
  if(length(tumors_below_threshold) > 0 || length(tumors_needing_group_average_rates) > 0) {
    mean_trinuc_prop = colMeans(trinuc_proportion_matrix[tumors_above_threshold, , drop = F])

    # convert to data frame because that's what deconstructSigs wants
    mean_trinuc_prop = as.data.frame(t(mean_trinuc_prop))
    
    
    # Hypermutation signatures are presumed absent from tumors subthreshold number of mutations, but may be present in the above-threshold
    # tumors. Therefore, they will be treated like artifact signatures: included in deconstructSigs signature weight calculation, but normalized
    # out when calculating the "true" relative trinucleotide SNV mutation rates.
    # However, they are left in for the assume_identical_mutational_processes method in the spirit of assuming all 
    # tumors have same mutational processes.
    mean_calc_artifact_signatures = artifact_signatures
    if (v3_hypermutation_rules && assume_identical_mutational_processes == FALSE) {
      mean_calc_artifact_signatures = unique(c(artifact_signatures, cosmic_v3_modest_hm_sigs, cosmic_v3_highly_hm_sigs))
    }
    
    rownames(mean_trinuc_prop) = "mean" # deconstructSigs crashes unless a rowname is supplied here
    mean_ds <- cancereffectsizeR::run_deconstructSigs(tumor_trinuc_counts = mean_trinuc_prop, signatures_df = signatures, 
                                                      signatures_to_remove = signatures_to_remove, tri.counts.method = "default",
                                                      artifact_signatures = mean_calc_artifact_signatures)
    mean_weights <- mean_ds$weights
    mean_trinuc_prop = as.numeric(mean_trinuc_prop) # convert back to numeric for insertion into trinuc_proportion_matrix
    # this should never happen
    if(all(mean_weights == 0)) {
      stop("Somehow, mean signature weights across all well-mutated tumors are all zero.")
    }
    if (assume_identical_mutational_processes) {
      for(i in 1:nrow(trinuc_proportion_matrix)) {
        trinuc_proportion_matrix[i,] = mean_trinuc_prop
      }
    } else {
      for (tumor in tumors_below_threshold) {
        # handle tumors with 0 non-recurrent SNVs
        if(is.na(substitution_counts[tumor])) {
          trinuc_proportion_matrix[tumor, ] = mean_trinuc_prop
        } else {
          own_weighting = substitution_counts[tumor] / sig_averaging_threshold
          group_weighting = 1 - own_weighting
          trinuc_proportion_matrix[tumor, ] = trinuc_proportion_matrix[tumor, ] * own_weighting + mean_trinuc_prop * group_weighting
          new_ds = signatures_output_list[[tumor]]
          new_ds$weights = new_ds$weights * own_weighting + mean_weights * group_weighting
          new_ds$product = trinuc_proportion_matrix[tumor, , drop=F] 
          signatures_output_list[[tumor]] = new_ds
        }
      }  
    }
  }

  # TGS tumors, tumors with zeroed-out weights, and tumors with no non-recurrent SNVs get assigned group-average weights
  for (tumor in tumors_needing_group_average_rates) {
    trinuc_proportion_matrix[tumor, ] = mean_trinuc_prop
  }

  # Update CESAnalysis status and add in trinuc results
  trinuc_method = paste0("Calculated using ", signature_set_name, " signatures")
  if (assume_identical_mutational_processes) {
    trinuc_method = paste0(trinuc_method, " (assume_identical_mutation_processes = TRUE)")
  }
  
  cesa@status[["trinucleotide mutation rates"]] = trinuc_method
  cesa@status[["gene mutation rates"]] = "uncalculated (run gene_level_mutation_rates)"
  cesa@trinucleotide_mutation_weights = list(trinuc_proportion_matrix=trinuc_proportion_matrix,
                                             signatures_output_list=signatures_output_list)
  return(cesa)
}


