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
#' in the signature definitions. The following columns allow special functionality:
#' \itemize{ 
#' \item Etiology: Known or hypothesized mutational processes underlying the signature. Used 
#' for human-readable tables and plots, so best to enter something like "Unknown" rather
#' than leaving any entries empty or NA
#' \item Likely_Artifact (logical T/F): Marks signatures that are believed to
#' derive from sample processing, sequencing, calling error, or other non-biological
#' sources. cancereffectsizeR adjusts for artifact signatures when inferring relative
#' trinucleotide mutation rates.
#' \item Exome_Min: Minimum number of mutations a WES sample
#' must have for the presence of the signature to be plausible. This information is used to
#' prevent hypermutation signatures from being found in tumors with few mutations. Can be
#' left NA or 0 for non-hypermutation signatures. If this column is present, Genome_Min
#' must be present and always greater than or equal to Exome_Min. 
#' \item Genome_Min:
#' Minimum number of mutations a WGS sample must have for the presence of the signature to
#' be plausible. This information is used to prevent hypermutation signatures from being found
#' in tumors with few mutations. Can be left NA or 0 for non-hypermutation signatures. If
#' this column is present, Exome_Min must be present and always less than or equal to
#' Genome_Min.
#' }
#' If you don't have any metadata available for your signature set, an empty data table
#' can also be supplied. For a template signature set object, run 
#' \code{sig_set = get_ces_signature_set("ces.refset.hg19", "COSMIC_v3.2")}.
#'
#' @param cesa CESAnalysis object
#' @param signature_set name of built-in signature set (see \code{list_ces_signature_sets()}), or a custom signature set (see details)
#' @param signatures_to_remove specify any signatures to exclude from analysis; use \code{suggest_cosmic_signatures_to_remove()} for advice on COSMIC signatures
#' @param sample_group vector of sample group(s) to calculate rates for (default all)
#' @param cores how many cores to use for processing tumors in parallel
#' @param sig_averaging_threshold Mutational threshold (default 50) that determines which
#'   tumors inform the calculation of group-average signature weights. When
#'   assume_identical_mutational_processes == FALSE (the default), these group averages
#'   are blended into the signature weights of sub-threshold tumors.
#' @param assume_identical_mutational_processes use well-mutated tumors (those with number
#'   of eligible mutations meeting sig_averaging_threshold) to calculate group average
#'   signature weights, and assign these to all tumors
#' @param use_dS_exome2genome historical dev option (don't use)
#' @return CESAnalysis with sample-specific signature weights and inferred
#'   trinucleotide-context-specific relative mutation rates. Note that tumors with few
#'   mutations (group_avg_blended == TRUE in the signature weights tables) have weights
#'   influenced by the average weights of well-mutated tumors; you may want to exclude
#'   these from some mutational signature analyses.
#' @export
#'
#'
#'
trinuc_mutation_rates <- function(cesa,
                                  signature_set = NULL,
                                  signatures_to_remove = character(),
                                  sample_group = NULL,
                                  cores = 1,
                                  assume_identical_mutational_processes = FALSE,
                                  sig_averaging_threshold = 50,
                                  use_dS_exome2genome = FALSE
                                  ){  
  if(is.null(cesa) || ! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be a CESAnalysis object", call. = F)
  }
  cesa = update_cesa_history(cesa, match.call())
  
  if(cesa@maf[, .N] == 0) {
    stop("No MAF data in the CESAnalysis", call. = F)
  }
  
  if(! is.null(sample_group) & ! is(sample_group, "character")) {
    stop("sample_group should be character")
  }
  sample_group = unique(sample_group)
  if (! all(sample_group %in% cesa@groups)) {
    possible_groups = setdiff(cesa@groups, "stageless") # historical
    if (length(possible_groups) == 0) {
      stop("The CESAnalysis has no user-defined sample groups, so you can't use the sample_group parameter.")
    }
    stop("Unrecognized sample groups. Those defined in the CESAnalysis are ", paste(possible_groups, sep = ", "), ".")
  }
  
  if (! is.null(sample_group)) {
    curr_sample_info = cesa@samples[group %in% sample_group]
  } else {
    curr_sample_info = cesa@samples
  }
  if (curr_sample_info[, .N] == 0) {
    stop("No selected samples to run (no samples in the chosen group?).")
  }
  
  bsg = get_cesa_bsg(cesa)
  
  if(all(curr_sample_info$coverage == "targeted")) {
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
  if(is(signature_set, "character")) {
    if (length(signature_set) != 1) {
      stop("signature_set should be 1-length character; run list_ces_signature_sets() for options.\n(Or, it can a custom signature set; see docs.)", call. = F)
    }
    signature_set_data = get_ces_signature_set(cesa@ref_key, signature_set)
  } else if(is(signature_set, "list")) {
    validate_signature_set(signature_set)
    signature_set_data = signature_set
  } else {
    stop("signature_set should be type character; run list_ces_signature_sets() for options.\n",
         "(Or, it can be a custom signature set; see docs.)")
  }
  signature_set_name = signature_set_data$name
  signatures = signature_set_data$signatures
  signature_metadata = signature_set_data$meta
  
  # columns Exome_Min and Genome_Min always go together in metadata, per signature set validation rules
  # If they're present, we'll enforce their signature mutation count minimums for each tumor
  mutation_count_rules = "Exome_Min" %in% colnames(signature_metadata)

  # If running with v3.0/v3.1 COSMIC signature set, help out the user by dropping signatures from future releases
  if(signature_set_name %in% c("COSMIC v3", "COSMIC v3.1")) {
    if (require("ces.refset.hg19", quietly = T)) {
      latest_set = ces.refset.hg19$signatures$COSMIC_v3.2$signatures
      if (is.null(latest_set)) {
        # If refset package is < v1.1.1
        latest_set = ces.refset.hg19$signatures$COSMIC_v3.1$signatures
      }
      future_signatures = setdiff(rownames(latest_set), rownames(signatures))
      signatures_to_remove = signatures_to_remove[! signatures_to_remove %in% future_signatures]
    }
  }
  
  # Save signature set to CESAnalysis
  # If already ran with different sample group, make sure signature set data is the same
  if (! is.null(cesa@advanced$snv_signatures)) {
    if (! all.equal(cesa@advanced$snv_signatures, signature_set_data, check.attributes = F)) {
      stop("Signature set does not exactly match the one previously used with this CESAnalysis.")
    }
  } else {
    cesa@advanced$snv_signatures = signature_set_data
  }
  
  
  
  # Put columns of sgnature data.frame into canonical deconstructSigs order (to match up with exome count order, etc.)
  signatures = signatures[, deconstructSigs_trinuc_string]
  
  # make sure all signatures that user has requested removed are actually valid signatures
  if(length(signatures_to_remove) > 0) {
    if(any(! signatures_to_remove %in% rownames(signatures))) {
      stop("One or more signatures in signatures_to_remove are not in the signature set.", call. = F)
    }
  }

  # keeping TGS data (and data from other groups, if applicable) in the ds_maf until after recurrency testing
  ds_maf = cesa@maf[variant_type == "snv"]

  # remove all recurrent SNVs (SNVs appearing in more than one sample)
  duplicated_vec_first <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)])
  duplicated_vec_last <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)],fromLast=T)
  duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
  if (length(duplicated_vec_pos) > 0) {
    ds_maf <- ds_maf[-duplicated_vec_pos,]
  }

  setkey(curr_sample_info, "Unique_Patient_Identifier")
  all_tumors = curr_sample_info$Unique_Patient_Identifier # may include some tumors with no SNVs, or with only recurrent SNVs
  
  # Check that tumors aren't already present in trinuc rates / weights tables
  # Since you can't have a signature weight entry without also being in rates matrix, just check rates matrix
  if(any(all_tumors %in% rownames(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix))) {
    stop("Trinucleotide-context-specific mutation rates have already been calculated for some or all of these tumors.")
  }
  
  
  # can only find signatures in tumors with exome/genome data that have non-recurrent SNVs
  tumors_eligible_for_trinuc_calc = intersect(curr_sample_info[coverage != "targeted", Unique_Patient_Identifier], unique(ds_maf$Unique_Patient_Identifier))
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
  tri_counts_by_gr = list()
  exome_cov_to_calc = curr_sample_info[coverage == "exome", unique(covered_regions)]
  genome_cov_to_calc = curr_sample_info[coverage == "genome", setdiff(unique(covered_regions), "genome")]
  grs_to_calc = c(cesa@coverage$exome[exome_cov_to_calc], cesa@coverage$genome[genome_cov_to_calc])
  gr_names = names(grs_to_calc)
  for (gr_name in gr_names) {
    if (gr_name %in% c("exome", "exome+")) {
      if (use_dS_exome2genome) {
        data("tri.counts.exome", package = "deconstructSigs", envir = environment())
        tri_counts_by_gr[[gr_name]] = tri.counts.exome
      } else {
        tri_counts_by_gr[[gr_name]] = get_ref_data(cesa, "tri.counts.exome")
      }
    } else {
      covered_seq = BSgenome::getSeq(bsg, grs_to_calc[[gr_name]])
      tri_contexts = Biostrings::trinucleotideFrequency(covered_seq)
      tri_contexts = colSums(tri_contexts)
      
      # deconstructSigs_trinuc_string is internal in cancereffectsizeR
      # here, we need unique trinucleotide contexts (without mutations) in deconstructSigs ordering,
      # which is why we produce this by converting from their column headings
      context_names = unique(sub("\\[([ACTG]).*\\]", "\\1", deconstructSigs_trinuc_string))
      context_names = sort(context_names)
      reverse_complement_names = unique(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(context_names))))
      
      # reorder the counts as desired, then save as a data frame since that's what deconstructSigs wants
      tri_counts = tri_contexts[context_names] + tri_contexts[reverse_complement_names]
      tri_counts = data.frame(x = tri_counts)
      tri_counts_by_gr[[gr_name]] = tri_counts
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
  if (! "Likely_Artifact" %in% colnames(signature_metadata)) {
    msg = paste("Note: There is no Likely_Artifact column in the signature set metadata, so",
                "artifact accounting can't be done. If you know that some signatures reflect",
                "sequencing error or other artifacts, you should fix this.")
    pretty_message(msg)
  } else {
    artifact_signatures = signature_metadata[Likely_Artifact == TRUE, Signature]
    if(any(artifact_signatures %in% signatures_to_remove)) {
      warning("Warning: You have chosen to remove at least one artifact signature from analysis,\n",
              "which will change how artifact accounting behaves (usually, all artifact signatures should be left in).")
    }
  }

  # for parallelization, a function to process each tumor
  process_tumor = function(tumor_name) {
    tumor_trinuc_counts = as.data.frame(trinuc_breakdown_per_tumor[tumor_name,]) # deconstructSigs requires a data.frame
    num_variants = sum(tumor_trinuc_counts)


    current_sigs_to_remove = signatures_to_remove
    
    # Hypermutation signatures have Exome_Min and Genome_Min values in metadata that give the smallest
    # number of mutations a tumor can have and still reasonably have the mutational process present.
    # Here, remove signatures that require more mutations than the tumor has.
    curr_sample_cov = curr_sample_info[tumor_name, coverage]
    if (mutation_count_rules) {
      if (curr_sample_cov == "exome") {
        current_sigs_to_remove = union(current_sigs_to_remove, signature_metadata[num_variants < Exome_Min, Signature])
      } else if (curr_sample_cov == "genome") {
        current_sigs_to_remove = union(current_sigs_to_remove, signature_metadata[num_variants < Genome_Min, Signature])
      }
    }

    # Set normalization argument for deconstructSigs based on coverage
    covered_regions = curr_sample_info[tumor_name, covered_regions]
    if (covered_regions == "genome") {
      normalization = "default" # this actually means no normalization (since signatures and MAF coverage are both whole-genome)
    } else {
      normalization = tri.counts.genome / tri_counts_by_gr[[covered_regions]]
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
    
    
    # Signatures that "expect" more mutations than are present in any of the subthreshold
    # tumors will be treated like artifact signatures: included in deconstructSigs
    # signature weight calculation, but normalized out when calculating the "true"
    # relative trinucleotide SNV mutation rates. However, they are left in for the
    # assume_identical_mutational_processes method in the spirit of assuming all tumors
    # have the same mutational processes.
    mean_calc_artifact_signatures = artifact_signatures
    if (mutation_count_rules && assume_identical_mutational_processes == FALSE) {
      mean_calc_artifact_signatures = union(artifact_signatures, signature_metadata[sig_averaging_threshold < Exome_Min, Signature])
    }
    
    rownames(mean_trinuc_prop) = "mean" # deconstructSigs crashes unless a rowname is supplied here
    message("Determining group-average signatures from samples with at least ", sig_averaging_threshold, " SNVs...")
    mean_ds <- run_deconstructSigs(tumor_trinuc_counts = mean_trinuc_prop, signatures_df = signatures, 
                                                      signatures_to_remove = signatures_to_remove, tri.counts.method = "default",
                                                      artifact_signatures = mean_calc_artifact_signatures)
    
    mean_weights = mean_ds$adjusted_sig_output$weights
    mean_weights_raw = mean_ds$adjusted_sig_output$raw_weights
    mean_trinuc_prop = as.numeric(mean_trinuc_prop) # convert back to numeric for insertion into trinuc_proportion_matrix
    
    # TGS tumors, tumors with zeroed-out weights, and tumors with no non-recurrent SNVs get assigned group-average rates
    for (tumor in tumors_needing_group_average_rates) {
      trinuc_proportion_matrix[tumor, ] = mean_trinuc_prop
    }
    
    # this should never happen
    if(all(mean_weights == 0)) {
      stop("Somehow, mean signature weights across all samples are all zero.")
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
          if (sig_averaging_threshold == 0 || tumor %in% zeroed_out_tumors) {
            own_weighting = 0
          } else {
            own_weighting = substitution_counts[tumor] / sig_averaging_threshold
          }
          group_weighting = 1 - own_weighting
          mean_blended_trinuc_prop = trinuc_proportion_matrix[tumor, ] * own_weighting + mean_trinuc_prop * group_weighting
          ds_output = signatures_output_list[[tumor]]
          mean_blended_weights = ds_output$adjusted_sig_output$weights * own_weighting + mean_weights * group_weighting
          raw_mean_blended = ds_output$adjusted_sig_output$raw_weights * own_weighting + mean_weights_raw * group_weighting
          ds_output$mean_blended = list(weights = mean_blended_weights, trinuc_prop = mean_blended_trinuc_prop, raw_weights = raw_mean_blended)
          trinuc_proportion_matrix[tumor, ] = mean_blended_trinuc_prop
          signatures_output_list[[tumor]] = ds_output
        }
      }  
    }
  }
  
  # If this is the first setting of trinuc rates, create lists
  if (length(cesa@trinucleotide_mutation_weights) == 0) {
    cesa@trinucleotide_mutation_weights = list(trinuc_proportion_matrix=trinuc_proportion_matrix,
                                               signatures_output_list=signatures_output_list)
  } else {
    new_mat = rbind(trinuc_proportion_matrix, cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix)
    cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix = new_mat
    new_sig_out_list = c(signatures_output_list, cesa@trinucleotide_mutation_weights$signatures_output_list)
    cesa@trinucleotide_mutation_weights$signatures_output_list = new_sig_out_list
  }
  
  # Put together signature weight tables
  if (! assume_identical_mutational_processes) {
    # identify tumors with weights informed by well-mutated (above SNV threshold) tumors
    blended_tumors = names(which(sapply(signatures_output_list, function(x) ! is.null(x$mean_blended))))
    
    blended_weights = data.table(t(sapply(blended_tumors, function(x) as.numeric(signatures_output_list[[x]]$mean_blended$weights))), 
                                 keep.rownames = "Unique_Patient_Identifier")
    blended_raw_weights = data.table(t(sapply(blended_tumors, function(x) as.numeric(signatures_output_list[[x]]$mean_blended$raw_weights))), 
                                     keep.rownames = "Unique_Patient_Identifier")
    
    nonblended_tumors = setdiff(names(signatures_output_list), blended_tumors)
    above_threshold_weights = data.table(t(sapply(nonblended_tumors, function(x) as.numeric(signatures_output_list[[x]]$adjusted_sig_output$weights))), 
                                         keep.rownames = "Unique_Patient_Identifier")
    above_threshold_raw_weights = data.table(t(sapply(nonblended_tumors, function(x) as.numeric(signatures_output_list[[x]]$adjusted_sig_output$raw_weights))), 
                                             keep.rownames = "Unique_Patient_Identifier")
    
    # get all signature weights into a data table (one row per sample)
    # sig_table_raw includes artifact signature_weights
    sig_table = rbind(above_threshold_weights, blended_weights)
    sig_table_raw = rbind(above_threshold_raw_weights, blended_raw_weights)
    
    
    # use first sample to set column names to signature names
    colnames(sig_table)[2:ncol(sig_table)] = colnames(signatures_output_list[[1]]$adjusted_sig_output$weights)
    colnames(sig_table_raw) = colnames(sig_table)
    sig_table[, group_avg_blended := Unique_Patient_Identifier %in% blended_tumors]
    sig_table[, sig_extraction_snvs := as.numeric(substitution_counts[Unique_Patient_Identifier])] # otherwise will be "table" class
    total_snv_counts = cesa@maf[variant_type == "snv"][sig_table, .(total_snvs = .N), on = "Unique_Patient_Identifier", by = "Unique_Patient_Identifier"]
    sig_table = sig_table[total_snv_counts, on = "Unique_Patient_Identifier"]
    
    
    sig_table_raw = sig_table_raw[sig_table[,.(Unique_Patient_Identifier, group_avg_blended, sig_extraction_snvs, total_snvs)], ,
                                  on = "Unique_Patient_Identifier"]
    
    # to reflect that no SNVs informed inferred weights of zeroed-out tumors, set sig_extraction_snvs to 0
    # note that sig_table_raw keeps the actual counts used by dS
    sig_table[Unique_Patient_Identifier %in% zeroed_out_tumors, sig_extraction_snvs := 0] 
    
    
    tumors_without_data = setdiff(curr_sample_info$Unique_Patient_Identifier, sig_table$Unique_Patient_Identifier)
    num_to_add = length(tumors_without_data)
    if (num_to_add > 0) {
      group_avg_weights = as.numeric(mean_ds$adjusted_sig_output$weights)
      new_rows = matrix(nrow = num_to_add, data = rep.int(group_avg_weights, num_to_add), byrow = T)
      colnames(new_rows) = colnames(mean_ds$adjusted_sig_output$weights)
      total_snvs = cesa@maf[variant_type == "snv"][, .N, keyby = "Unique_Patient_Identifier"][tumors_without_data, N]
      total_snvs[is.na(total_snvs)] = 0
      new_table = data.table(Unique_Patient_Identifier = tumors_without_data, total_snvs = total_snvs, 
                             sig_extraction_snvs = 0, group_avg_blended = T)
      sig_table = rbind(sig_table, cbind(new_table, new_rows))
      
      # repeat for table that includes artifacts
      group_avg_weights_raw = as.numeric(mean_ds$adjusted_sig_out$raw_weights)
      new_rows_raw = matrix(nrow = num_to_add, data = rep.int(group_avg_weights_raw, num_to_add), byrow = T)
      colnames(new_rows_raw) = colnames(mean_ds$adjusted_sig_out$raw_weights)
      new_table[, sig_extraction_snvs := as.numeric(substitution_counts[Unique_Patient_Identifier])]
      new_table[is.na(sig_extraction_snvs), sig_extraction_snvs := 0]
      sig_table_raw = rbind(sig_table_raw, cbind(new_table, new_rows_raw))
    }
    setcolorder(sig_table, c("Unique_Patient_Identifier", "total_snvs", "sig_extraction_snvs", "group_avg_blended"))
    setcolorder(sig_table_raw, c("Unique_Patient_Identifier", "total_snvs", "sig_extraction_snvs", "group_avg_blended"))
    
    new_sig_table = rbind(cesa@trinucleotide_mutation_weights[["signature_weight_table"]], sig_table)
    cesa@trinucleotide_mutation_weights[["signature_weight_table"]] = new_sig_table
    new_sig_raw_table = rbind(cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]], sig_table_raw)
    cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]] = new_sig_raw_table
  }
  
  if(! is.null(mean_ds)) {
    if (is.null(sample_group)) {
      cesa@trinucleotide_mutation_weights[["group_average_dS_output"]][["all"]] = mean_ds
    } else {
      cesa@trinucleotide_mutation_weights[["group_average_dS_output"]][[paste(sample_group, collapse = ',')]] = mean_ds
    }
  }
  cesa@advanced$locked = T
  
  if(nrow(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix) == cesa@samples[, .N]) {
    cesa@advanced$trinuc_done = T
  }
  return(cesa)
}


