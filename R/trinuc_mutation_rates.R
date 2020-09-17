#' Calculate relative rates of trinucleotide-context-specific mutations by extracting underlying mutational processes
#'
#' This function calculates expected relative rates of trinucleotide-context-specific SNV
#' mutations within tumors by attributing SNVs to mutational processes represented
#' in mutation signature sets (such as COSMIC v3). deconstructSigs is used to extract mutational 
#' signature weights in each tumor. Tumors with targeted sequencing data are assigned the
#' average trinucleotide mutation rates calculated across all exome/genome data, 
#' which means that you need at least some exome or genome data to run.
#' 
#' To reduce the influence of selection on the estimation of relative trinucleotide mutation
#' rates, only non-recurrent SNVs (those that do not appear in more than one 
#' sample in the data set) are used.
#' 
#' A custom signature set should be given as a named three-item list, where "signatures"
#' is a pure data.frame with signature definitions, "name" is a 1-length character naming
#' the set, and "metadata" is a data.table with a "Signature" column that matches rownames
#' in the signature definitions. For artifact accounting to work, there should be a
#' TRUE/FALSE (logical) Likely_Artifact column; additional columns are optional. If you
#' don't have artifact information, you can also just supply an empty data.table. For a
#' template signature set object, run \code{sig_set = get_ces_signature_set("ces_hg19_v1",
#' "COSMIC_v3.1")}.
#'
#' @param cesa CESAnalysis object
#' @param signature_set name of built-in signature set (see \code{list_ces_signature_sets()}), or a custom signature set (see details)
#' @param signatures_to_remove specify any signatures to exclude from analysis; use \code{suggest_cosmic_signatures_to_remove()} for advice on COSMIC signatures
#' @param cores how many cores to use for processing tumors in parallel
#' @param assume_identical_mutational_processes use well-mutated tumors (those with number of eligible mutations meeting sig_averaging_threshold) 
#'   to calculate group average signature weights, and assign these to all tumors
#' @param sig_averaging_threshold Mutational threshold (default 50) that determines which tumors inform the
#'   calculation of group-average signature weightings. When assume_identical_mutational_processes == FALSE (the default), 
#'   these group averages are blended into the signature weightings of tumors with few mutations (those below the threshold).
#' @param cosmic_hypermutation_rules T/F on whether to follow mutation count rules outlined in https://doi.org/10.1101/322859 (COSMIC v3 manuscript)
#'   (only applies when running with COSMIC v3/v3.1 signatures)
#' @param artifact_accounting set false to disable special handling of artifact signatures (rarely recommended)
#' @param use_dS_exome2genome internal dev option (don't use)
#' @return CESAnalysis with expected relative trinucleotide mutation rates ($trinuc_rates) and a table of tumor-specific 
#'         signature weights ($mutational_signatures). Note that tumors with group_avg_blended == TRUE have signature
#'         weights influenced by the average weights of well-mutated tumors; you may want to exclude these from some analyses.
#' @export
#'
#'
#'
trinuc_mutation_rates <- function(cesa,
                                  signature_set = NULL,
                                  signatures_to_remove = character(),
                                  cores = 1,
                                  assume_identical_mutational_processes = FALSE,
                                  sig_averaging_threshold = 50,
                                  cosmic_hypermutation_rules = TRUE,
                                  artifact_accounting = TRUE,
                                  use_dS_exome2genome = FALSE
                                  ){  
  if(is.null(cesa) || ! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be a CESAnalysis object", call. = F)
  }
  cesa = update_cesa_history(cesa, match.call())
  
  # If for some reason reference data is no longer in the package environment, restore it
  if (! cesa@ref_key %in% ls(.ces_ref_data)) {
    preload_ref_data(cesa@ref_data_dir)
  }
  
  bsg = .ces_ref_data[[cesa@ref_key]]$genome
  
  if(cesa@maf[, .N] == 0) {
    stop("No MAF data in the CESAnalysis", call. = F)
  }
  
  if(all(cesa@samples$coverage == "targeted")) {
    stop("We can't estimate relative trinucleotide mutation rates without some exome/genome data in the CESAnalysis (all data is targeted sequencing).", call. = F)
  }
  
  if(! is.logical(assume_identical_mutational_processes) || length(assume_identical_mutational_processes) != 1 || is.na(assume_identical_mutational_processes)) {
    stop("Expected assume_identical_mutational_processes to be TRUE/FALSE (default is FALSE).", call. = F)
  }
  
  if(! is.numeric(sig_averaging_threshold) || length(sig_averaging_threshold) != 1 || is.na(sig_averaging_threshold) || sig_averaging_threshold < 0) {
    stop("Expected positive integer for sig_averaging_threshold", call. = F)
  }
  sig_averaging_threshold = as.integer(sig_averaging_threshold) # below, we test that at least one tumor meets this threshold
  
  if(! is(signatures_to_remove, "character")) {
    stop("signatures_to_remove should be a character vector", call. = F)
  }

  
  running_cosmic_v3 = FALSE # gets set to true when user loads COSMIC v3 or v3.1, so that hypermutation rules can be applied
  if(is(signature_set, "character")) {
    if (length(signature_set) != 1) {
      stop("signature_set should be 1-length character; run list_ces_signature_sets() for options.\n(Or, it can a custom signature set; see docs.)", call. = F)
    }
    
    signature_set_data = get_ces_signature_set(cesa@ref_key, signature_set)
    signature_set_name = signature_set_data$name
    signatures = signature_set_data$signatures
    signature_metadata = signature_set_data$meta

    if(signature_set_name %in% c("COSMIC v3", "COSMIC v3.1")) {
      running_cosmic_v3 = TRUE
      if (cosmic_hypermutation_rules) {
        message(crayon::black(paste0("Samples with many mutations will have hypermutation rules applied.\n",
                                     "(Disable with cosmic_hypermutation_rules=FALSE.)")))
      }
      
      if (signature_set_name == "COSMIC v3") {
        # automatically strip out v3.1 signatures when running v3
        new_in_3_1 = setdiff(rownames(get_ces_signature_set(cesa@ref_key, "COSMIC_v3.1")$signatures), rownames(signatures))
        signatures_to_remove = signatures_to_remove[! signatures_to_remove %in% new_in_3_1]
      }
    }
  } else if(is(signature_set, "list")) {
      signature_set_data = signature_set
      signature_set_name = signature_set_data$name
      signatures = signature_set_data$signatures
      signature_metadata = signature_set_data$meta
      if (! is(signature_set_name, "character") || length(signature_set_name) != 1 || ! is(signatures, "data.frame") ||
          ! is(signature_metadata, "data.table")) {
        stop("Improperly formatted custom signature set; see documentation.", call. = F)
      }
      if(is(signatures, "data.table") || is(signatures, "tbl")) {
        stop("For compatibility with deconstructSigs, signature definitions must be given as a pure data.frame (see docs).", call. = F)
      }
      if(! identical(sort(deconstructSigs_trinuc_string), sort(colnames(signatures)))) {
        tmp = paste0("\n", '"', paste(deconstructSigs_trinuc_string, collapse = '", "'), '"')
        cat("Expected signature definition column names:\n")
        writeLines(strwrap(tmp, indent = 4, exdent = 4))
        stop("Your signature definition data frame has improper column names.", call. = F)
      }
      # Validate signature metadata if it's not empty
      if (signature_metadata[, .N] > 0) {
        if (is.null(signature_metadata$Signature)) {
          stop("Signature metadata incorrectly formatted (see docs).")
        }
        if (any(! rownames(signatures) %in% signature_metadata$Signature)) {
          stop("Improperly formatted signature set: Some signatures in your signature definitions are missing from the metadata table.")
        }
        if(length(signature_metadata$Signature) != length(unique(signature_metadata$Signature))) {
          stop("Improperly formatted signature set: Some signatures are repeated in your signature metadata table")
        }
      }
  } else {
    stop("signature_set should be type character; run list_ces_signature_sets() for options.\n",
    "(Or, it can a custom signature set; see docs.)", call. = F)
  }
  
  # Save signature set to CESAnalysis (leaving out other attributes for now)
  cesa@advanced$snv_signatures = signature_set_data[c("name", "meta", "signatures")]
  
  # Put columns of sgnature data.frame into canonical deconstructSigs order (to match up with exome count order, etc.)
  signatures = signatures[, deconstructSigs_trinuc_string]
  
  # make sure all signatures that user has requested removed are actually valid signatures
  if(length(signatures_to_remove) > 0) {
    if(any(! signatures_to_remove %in% rownames(signatures))) {
      stop("One or more signatures in signatures_to_remove are not in the signature set.", call. = F)
    }
  }
  
  # can't apply our COSMIC v3/3.1 enhancements unless we're using those signatures, naturally
  if(! running_cosmic_v3) {
    cosmic_hypermutation_rules = FALSE
  }
  

  # keeping TGS data in the ds_maf until after recurrency testing
  ds_maf = cesa@maf[variant_type == "snv"]

  # remove all recurrent SNVs (SNVs appearing in more than one sample)
  duplicated_vec_first <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)])
  duplicated_vec_last <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)],fromLast=T)
  duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
  if (length(duplicated_vec_pos) > 0) {
    ds_maf <- ds_maf[-duplicated_vec_pos,]
  }

  setkey(cesa@samples, "Unique_Patient_Identifier")
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


  tri.counts.genome = get_ref_data(cesa, "tri.counts.genome")
  
  # for each exome coverage gr (besides default generic, which is pre-calculated), tabulate trinucs
  exome_counts_by_gr = list()
  exomes_to_calc = cesa@samples[coverage == "exome", unique(covered_regions)]
  for (exome_name in exomes_to_calc) {
    if (exome_name %in% c("exome", "exome+")) {
      if (use_dS_exome2genome) {
        data("tri.counts.exome", package = "deconstructSigs")
        exome_counts_by_gr[[exome_name]] = tri.counts.exome
      } else {
        exome_counts_by_gr[[exome_name]] = get_ref_data(cesa, "tri.counts.exome")
      }
    } else {
      exome_seq = getSeq(bsg, cesa@coverage$exome[[exome_name]])
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
  trinuc = BSgenome::getSeq(bsg, ds_maf$Chromosome, ds_maf$Start_Position - 1, ds_maf$Start_Position + 1, as.character = T)
  
  # internal table converts trinuc/mut (e.g., GTA:C) into deconstructSigs format ("G[T>C]A")
  ds_muts = factor(deconstructSigs_notations[.(trinuc, ds_maf$Tumor_Allele), deconstructSigs_ID], levels = deconstructSigs_trinuc_string)
  
  # mysteriously convert two-way table to data frame
  tmp = table(ds_maf$Unique_Patient_Identifier, ds_muts)
  counts = apply(tmp, 2, rbind)
  
  # edge case: when just 1 sample in data set, counts comes back as integer vector but need matrix
  if(! is.matrix(counts)) {
    counts = t(as.matrix(counts))
  }
  rownames(counts) = rownames(tmp)
  trinuc_breakdown_per_tumor = as.data.frame(counts)

  artifact_signatures = NULL
  if (artifact_accounting) {
    # COSMIC v3 artifact signatures are read in from package data 
    if (! "Likely_Artifact" %in% colnames(signature_metadata)) {
      message("Note: There is no Likely_Artifact column in the signature set metadata, so\n",
              "artifact accounting can't be done. If you know some signatures reflect\n",
              "sequencing error or other artifacts, you should fix this.")
    } else if (! is(signature_metadata$Likely_Artifact, "logical")) {
      stop("Improperly formatted signature set metadata: column Likely_Artifact should be logical.", call. = F)
    } else {
      artifact_signatures = signature_metadata[Likely_Artifact == TRUE, Signature]
      if(any(artifact_signatures %in% signatures_to_remove)) {
        warning("Warning: You are have chosen to remove at least one sequencing-artifact-associated signature from analysis,\n",
                "which will change how artifact accounting behaves (usually, all artifact signatures should be left in).")
      }
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
    if (cosmic_hypermutation_rules) {
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
    return(run_deconstructSigs(tumor_trinuc_counts = tumor_trinuc_counts, tri.counts.method = normalization,
                                                signatures_df = signatures, signatures_to_remove = current_sigs_to_remove, artifact_signatures = artifact_signatures))
  }

  
  
  # store results
  # matrix with rows = tumors, columns = relative rates of mutation for each trinuc SNV (starts empty)
  trinuc_proportion_matrix <- matrix(data = NA, nrow = length(all_tumors), ncol = ncol(trinuc_breakdown_per_tumor),
                                     dimnames = list(all_tumors, colnames(trinuc_breakdown_per_tumor)))
  signatures_output_list = list()
  
  
  zeroed_out_tumors = character() # for tumors that get all zero signature weights (rare)
  message("Extracting mutational signatures from SNVs...")
  ds_output = pbapply::pblapply(tumors_eligible_for_trinuc_calc, process_tumor, cl = cores)
  for(i in 1:length(ds_output)) {
    tumor_name = tumors_eligible_for_trinuc_calc[i]
    signatures_output = ds_output[[i]]
    signatures_output_list[[tumor_name]] = signatures_output
    
    # rarely, weights will come out all 0, so trinuc_prop will be NULL; these tumors are "zeroed-out"
    if(is.null(signatures_output$adjusted_sig_output$trinuc_prop)) {
      zeroed_out_tumors = c(zeroed_out_tumors, tumor_name)
      substitution_counts[[tumor_name]] = 0 # no SNVs inform trinuc rates in these tumors
    } else {
      trinuc_proportion_matrix[tumor_name, ] = signatures_output$adjusted_sig_output$trinuc_prop
    }
  }
  # in the (very rare) occurrence of zeroed-out tumors, set them aside to be assigned group average signature weightings
  tumors_above_threshold = setdiff(tumors_above_threshold, zeroed_out_tumors)
  tumors_below_threshold = union(tumors_below_threshold, zeroed_out_tumors)
  tumors_needing_group_average_rates = union(tumors_needing_group_average_rates, zeroed_out_tumors)

  if(assume_identical_mutational_processes == TRUE) {
    tumors_needing_group_average_rates = union(tumors_needing_group_average_rates, tumors_eligible_for_trinuc_calc)
  }
  
  # If any samples have subthreshold number of mutations, determine weightings for those samples using method specified by user
  mean_ds = NULL
  if(length(tumors_below_threshold) > 0 || length(tumors_needing_group_average_rates) > 0) {
    if(length(tumors_above_threshold) == 0) {
      stop(paste0("No tumors have enough mutations to inform group-average signature extraction\n", 
                  "(or, rarely, all such tumors have all their mutations attributed to artifacts)."))
    }
    mean_trinuc_prop = colMeans(trinuc_proportion_matrix[tumors_above_threshold, , drop = F])

    # convert to data frame because that's what deconstructSigs wants
    mean_trinuc_prop = as.data.frame(t(mean_trinuc_prop))
    
    
    # Hypermutation signatures are presumed absent from tumors subthreshold number of mutations, but may be present in the above-threshold
    # tumors. Therefore, they will be treated like artifact signatures: included in deconstructSigs signature weight calculation, but normalized
    # out when calculating the "true" relative trinucleotide SNV mutation rates.
    # However, they are left in for the assume_identical_mutational_processes method in the spirit of assuming all 
    # tumors have the same mutational processes.
    mean_calc_artifact_signatures = artifact_signatures
    if (cosmic_hypermutation_rules && assume_identical_mutational_processes == FALSE) {
      mean_calc_artifact_signatures = unique(c(artifact_signatures, cosmic_v3_modest_hm_sigs, cosmic_v3_highly_hm_sigs))
    }
    
    rownames(mean_trinuc_prop) = "mean" # deconstructSigs crashes unless a rowname is supplied here
    message("Determining group-average signatures from well-mutated tumors....")
    mean_ds <- run_deconstructSigs(tumor_trinuc_counts = mean_trinuc_prop, signatures_df = signatures, 
                                                      signatures_to_remove = signatures_to_remove, tri.counts.method = "default",
                                                      artifact_signatures = mean_calc_artifact_signatures)
    
    mean_weights <- mean_ds$adjusted_sig_output$weights
    mean_trinuc_prop = as.numeric(mean_trinuc_prop) # convert back to numeric for insertion into trinuc_proportion_matrix
    
    # TGS tumors, tumors with zeroed-out weights, and tumors with no non-recurrent SNVs get assigned group-average rates
    for (tumor in tumors_needing_group_average_rates) {
      trinuc_proportion_matrix[tumor, ] = mean_trinuc_prop
    }
    
    # this should never happen
    if(all(mean_weights == 0)) {
      stop("Somehow, mean signature weights across all well-mutated tumors are all zero.")
    }
    if (assume_identical_mutational_processes) {
      # when assume_identical_mutational_processes is TRUE, this is where trinuc_proportion_matrix gets built
      for(i in 1:nrow(trinuc_proportion_matrix)) {
        trinuc_proportion_matrix[i,] = mean_trinuc_prop
      }
    } else {
      for (tumor in tumors_below_threshold) {
        # handle tumors with 0 non-recurrent SNVs
        if(is.na(substitution_counts[tumor])) {
          trinuc_proportion_matrix[tumor, ] = mean_trinuc_prop
        } else {
          # zeroed-out tumors are sub-threshold even when sig_averaging_threshold is 0; take no weight from these
          if (sig_averaging_threshold == 0) {
            own_weighting = 0
          } else {
            own_weighting = substitution_counts[tumor] / sig_averaging_threshold
          }
          group_weighting = 1 - own_weighting
          mean_blended_trinuc_prop = trinuc_proportion_matrix[tumor, ] * own_weighting + mean_trinuc_prop * group_weighting
          ds_output = signatures_output_list[[tumor]]
          mean_blended_weights = ds_output$adjusted_sig_output$weights * own_weighting + mean_weights * group_weighting
          ds_output$mean_blended = list(weights = mean_blended_weights, trinuc_prop = mean_blended_trinuc_prop)
          trinuc_proportion_matrix[tumor, ] = mean_blended_trinuc_prop
          signatures_output_list[[tumor]] = ds_output
        }
      }  
    }
  }

  cesa@trinucleotide_mutation_weights = list(trinuc_proportion_matrix=trinuc_proportion_matrix,
                                             signatures_output_list=signatures_output_list)
  
  if (! assume_identical_mutational_processes) {
    # identify tumors with weights informed by well-mutated (above SNV threshold) tumors
    blended_tumors = names(which(sapply(signatures_output_list, function(x) ! is.null(x$mean_blended))))
    
    blended_weights = data.table(t(sapply(blended_tumors, function(x) as.numeric(signatures_output_list[[x]]$mean_blended$weights))), 
                                 keep.rownames = "Unique_Patient_Identifier")
    
    nonblended_tumors = setdiff(names(signatures_output_list), blended_tumors)
    above_threshold_weights = data.table(t(sapply(nonblended_tumors, function(x) as.numeric(signatures_output_list[[x]]$adjusted_sig_output$weights))), 
                                         keep.rownames = "Unique_Patient_Identifier")
    
    # get all signature weights into a data table (one row per sample)
    sig_table = rbind(above_threshold_weights, blended_weights)
    
    # use first sample to set column names to signature names
    colnames(sig_table)[2:ncol(sig_table)] = colnames(signatures_output_list[[1]]$adjusted_sig_output$weights)
    sig_table[, group_avg_blended := Unique_Patient_Identifier %in% blended_tumors]
    sig_table[, sig_extraction_snvs := as.numeric(substitution_counts[Unique_Patient_Identifier])] # otherwise will be "table" class
    
    total_snv_counts = cesa@maf[variant_type == "snv"][sig_table, .(total_snvs = .N), on = "Unique_Patient_Identifier", by = "Unique_Patient_Identifier"]
    sig_table = sig_table[total_snv_counts, on = "Unique_Patient_Identifier"]
    
    
    tumors_without_data = setdiff(cesa@samples$Unique_Patient_Identifier, sig_table$Unique_Patient_Identifier)
    num_to_add = length(tumors_without_data)
    if (num_to_add > 0) {
      group_avg_weights = as.numeric(mean_ds$adjusted_sig_output$weights)
      new_rows = matrix(nrow = num_to_add, data = rep.int(group_avg_weights, num_to_add), byrow = T)
      colnames(new_rows) = colnames(mean_ds$adjusted_sig_output$weights)
      total_snvs = cesa@maf[variant_type == "snv"][, .N, keyby = "Unique_Patient_Identifier"][tumors_without_data, N]
      total_snvs[is.na(total_snvs)] = 0
      new_table = data.table(Unique_Patient_Identifier = tumors_without_data, total_snvs = total_snvs, 
                             sig_extraction_snvs = 0, group_avg_blended = T)
      new_table = cbind(new_table, new_rows)
      sig_table = rbind(sig_table, new_table)
    }
    setcolorder(sig_table, c("Unique_Patient_Identifier", "total_snvs", "sig_extraction_snvs", "group_avg_blended"))
    cesa@trinucleotide_mutation_weights[["signature_weight_table"]] = sig_table
  }
  
  if(! is.null(mean_ds)) {
    cesa@trinucleotide_mutation_weights[["group_average_dS_output"]] = mean_ds
  }
  return(cesa)
}


