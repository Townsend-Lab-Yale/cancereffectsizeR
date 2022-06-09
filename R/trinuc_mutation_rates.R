#' Calculate relative rates of trinucleotide-context-specific mutations by extracting underlying mutational processes
#'
#' This function calculates expected relative rates of trinucleotide-context-specific SNV
#' mutations within tumors by attributing SNVs to mutational processes represented in
#' mutation signature sets (such as "COSMIC v3.2"). Signature extraction can be done with
#' MutationalPatterns (default) or deconstructSigs. Tumors with targeted sequencing data
#' are assigned the average trinucleotide mutation rates calculated across all
#' exome/genome data, which means that you need at least some exome or genome data to run.
#' 
#' To reduce the influence of selection on the estimation of relative trinucleotide mutation
#' rates, only non-recurrent SNVs (those that do not appear in more than one 
#' sample in the current run) are used.
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
#' @param signature_set Name of built-in signature set (see
#'   \code{list_ces_signature_sets()}), or a custom signature set (see details)
#' @param signature_exclusions Specify any signatures to exclude from analysis; use
#'   \code{suggest_cosmic_signature_exclusions()} for advice on COSMIC signatures
#' @param samples Which samples to include in the current run. Defaults to all samples.
#'   Can be a vector of Unique_Patient_Identifiers, or a data.table containing rows from
#'   the CESAnalysis sample table.
#' @param cores How many cores to use for processing tumors in parallel.
#' @param signature_extractor One of "MutationalPatterns" (default) or "deconstructSigs".
#' @param mp_strict_args Named list of arguments to pass to MutationalPatterns'
#'   fit_to_signatures_strict function. Note that mut_matrix and signatures arguments are
#'   generated automatically, and that if you'd rather not use the strict method, you can
#'   emulate fit_to_signatures() by setting max_delta = 0.
#' @param bootstrap_mutations T/F (default FALSE). Instead of using actual SNV counts for
#'   the samples, do a single bootstrap sampling of each sample (in other words, run
#'   MutationalPatterns::fit_to_signatures_bootstrapped() with \code{n_boot=1}). This can
#'   be useful if you intend to run this function multiple times to get a distribution of
#'   signature attributions and trinuc rates (and downstream cancer effect sizes). This
#'   option may be replaced with more thorough support for bootstrapping in the future.
#' @param sig_averaging_threshold Mutation prevalence threshold (default 50) that
#'   determines which tumors inform the calculation of group-average signature weights.
#'   When assume_identical_mutational_processes == FALSE (the default), these group
#'   averages are blended into the signature weights of sub-threshold tumors.
#' @param assume_identical_mutational_processes use well-mutated tumors (those with number
#'   of eligible mutations meeting sig_averaging_threshold) to calculate group average
#'   signature weights, and assign these (and implied trinucleotide mutation rates) to all tumors
#' @param sample_group (Deprecated; use samples.) Vector of sample group(s) to calculate rates for.
#' @return CESAnalysis with sample-specific signature weights and inferred
#'   trinucleotide-context-specific relative mutation rates. The snv_counts matrix gives
#'   the counts of SNVs in each trinucleotide context for all samples in the CESAnalysis.
#'   (While recurrent mutations are excluded from signature analysis in
#'   trinuc_mutation_rates(), they are present in snv_counts for completeness.) The
#'   snv_counts matrix, produced by `trinuc_snv_counts()`, can be fed directly into
#'   MutationalPatterns if you wish to run your own extended signature analysis.
#'   
#'   The raw_attributions table contains signature attributions as produced by
#'   MutationalPatterns or deconstructSigs. The biological_weights table has several
#'   differences:
#'   \itemize{ 
#'   \item Weights for signatures associated with artifactual (as opposed to biological) processes are set to zero.
#'   \item The remaining weights are normalized to sum to 1. 
#'   \item Tumors with few mutations (defined by `sig_averaging_threshold`, default = 50)
#'   have their weights redefined using a blend of their original weights and weights
#'   derived from running signature extraction en masse on tumors with above-threshold
#'   mutation counts. These samples are identifiable by filtering the table on
#'   `group_avg_blended == TRUE`, and we recommend excluding them from most downstream
#'   signature analysis. (These weights are useful as a best
#'   guess of mutational processes, but they shouldn't be reported in any way that
#'   implies their independence from the group-average weights.)
#'   
#'   Biological weights can be interpreted as follows: Out of all the mutations caused by
#'   biological processes represented in the signatures, the proportion of mutations
#'   attributed to given signature is its weight.
#'   
#'   Either signature attributions table can be converted into the matrix format used by MutationalPatterns 
#'   with `convert_signature_weights_for_mp()`.
#'   }
#' @export
#'
trinuc_mutation_rates <- function(cesa,
                                  signature_set = NULL,
                                  signature_exclusions = character(),
                                  samples = character(),
                                  cores = 1,
                                  signature_extractor = "MutationalPatterns",
                                  mp_strict_args = list(),
                                  bootstrap_mutations = FALSE,
                                  assume_identical_mutational_processes = FALSE,
                                  sample_group = NULL,
                                  sig_averaging_threshold = 50,
                                  signatures_to_remove = NULL) {  
  
  # Documented argumented name has changed from signatures_to_remove to signatures_to_exclude
  if (is.null(signatures_to_remove)) {
    signatures_to_remove = signature_exclusions
  } else {
    message("FYI, signatures_to_remove has been renamed signature_exclusions, and will be removed eventually.")
  }
  if (! is(signature_extractor, "character") || length(signature_extractor) != 1 || 
      ! tolower(signature_extractor) %in% c('mutationalpatterns', 'deconstructsigs')) {
    stop("signature_extractor must be 'MutationalPatterns' or 'deconstructSigs'. (Or, use set_signature_weights to ",
         "load signature weights determined with some other method.)")
  }
  signature_extractor = ifelse(tolower(signature_extractor) == "mutationalpatterns", 'MutationalPatterns', 'deconstructSigs')
  
  if (! requireNamespace(signature_extractor, quietly = T)) {
    stop(signature_extractor, " is the chosen signature method, but it doesn't seem to be installed..")
  }
  
  # fit_to_signatures_strict and many other MP features appear just before v3 release; for simplicity, we won't support older.
  if(signature_extractor == "MutationalPatterns" && packageVersion(signature_extractor) < as.package_version("2.99.4")) {
    stop("Please update MutationalPatterns to v2.99.4 or later.")
  }
  
  if (! is(mp_strict_args, "list") || uniqueN(names(mp_strict_args)) != length(mp_strict_args)) {
    stop("mp_strict_args should a named list of arguments to pass.")
  }
  if (any(c("signatures", "mut_matrix") %in% names(mp_strict_args))) {
    stop("mp_strict_args: You can't supply mut_matrix or signatures because this function generates them automatically.")
  }
  
  if('n_boot' %in% names(mp_strict_args)) {
    stop("Sorry, you can't set n_boot with mp_strict_args.")
  }
  
  if(! is.logical(bootstrap_mutations) || length(bootstrap_mutations) != 1) {
    stop("bootstrap_mutations should be TRUE/FALSE.")
  }
  if(bootstrap_mutations == TRUE && signature_extractor != 'MutationalPatterns') {
    stop("bootstrap_mutations requires use of MutationalPatterns.")
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
  
  
  if(is.null(cesa) || ! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be a CESAnalysis object", call. = F)
  }
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  
  if(cesa@maf[, .N] == 0) {
    stop("No MAF data in the CESAnalysis", call. = F)
  }
  
  # get validated subset of cesa@samples for the input samples
  curr_sample_info = select_samples(cesa, samples)
  
  # handle (deprecated) sample_group method of sample selection
  if(! is.null(sample_group)) {
    if(! identical(samples, character())) {
      stop("Use just one of samples and sample_group (use samples, as sample_group is deprecated).")
    }
    warning("sample_group is deprecated and will be removed; \"samples\" is more flexible")
    if(! is(sample_group, "character")) {
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
    curr_sample_info = cesa@samples[group %in% sample_group]
    if (curr_sample_info[, .N] == 0) {
      stop("No selected samples to run (no samples in the chosen group?).")
    }
  }


  if(! all(is.na(curr_sample_info$sig_analysis_grp))) {
    msg = pretty_message(paste0("Some samples in current input have already been run. Use clear_trinuc_rates_and_signatures() if ",
                                "you want to re-run them."), emit = F)
    stop(msg)
  }
  
  if(all(curr_sample_info$coverage == "targeted")) {
    stop("We can't estimate relative trinucleotide mutation rates without some exome/genome data in the CESAnalysis (all data is targeted sequencing).", call. = F)
  }
  
  signature_set_name = signature_set_data$name
  signatures = signature_set_data$signatures
  signature_metadata = signature_set_data$meta
  
  # columns Exome_Min and Genome_Min always go together in metadata, per signature set validation rules
  # If they're present, we'll enforce their signature mutation count minimums for each tumor
  mutation_count_rules = "Exome_Min" %in% names(signature_metadata)

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
  
  # Save a copy of signature set in the CESAnalysis
  # If already ran with different sample group, make sure signature set data is the same
  previously_saved = cesa@advanced$snv_signatures[[signature_set_name]]
  if (is.null(previously_saved)) {
    cesa@advanced$snv_signatures[[signature_set_name]] = copy(signature_set_data)
  } else {
    # Okay if attributes don't match, probably (e.g., one set's table may be keyed from interactive use)
    if(! all.equal(previously_saved, signature_set_data, check.attributes = F)) {
      stop("A signature set with the same name has already been used in the CESAnalysis, but it is not ",
           "identical to the input set.")
    }
  }
  
  # Put columns of signature data.frame into canonical deconstructSigs order (to match up with exome count order, etc.)
  signatures = signatures[, deconstructSigs_trinuc_string]
  
  # make sure all signatures that user has requested removed are actually valid signatures
  if(length(signatures_to_remove) > 0) {
    if(any(! signatures_to_remove %in% rownames(signatures))) {
      stop("One or more signatures in signatures_to_remove are not in the signature set.", call. = F)
    }
  }

  setkey(curr_sample_info, "Unique_Patient_Identifier")
  all_tumors = curr_sample_info$Unique_Patient_Identifier # may include some tumors with no SNVs, or with only recurrent SNVs
  
  # Check that tumors aren't already present in trinuc rates / weights tables
  # Since you can't have a signature weight entry without also being in rates matrix, just check rates matrix
  if(any(all_tumors %in% rownames(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix))) {
    stop("Trinucleotide-context-specific mutation rates have already been calculated for some or all of these tumors.")
  }
  
  # Get SNV counts. Yes, always getting MutationalPatterns style, and then will convert to deconstructSigs format
  # later if needed.
  # First, save the full counts across all data for user access. Re-running if user has loaded more MAF data.
  full_trinuc_snv_counts = cesa@trinucleotide_mutation_weights$trinuc_snv_counts
  if (is.null(full_trinuc_snv_counts) || 
      ncol(full_trinuc_snv_counts) < cesa@maf[variant_type == 'snv', uniqueN(Unique_Patient_Identifier)]) {
    full_trinuc_snv_counts = trinuc_snv_counts(cesa@maf, genome = get_cesa_bsg(cesa), 
                                               exclude_recurrent = FALSE, style = 'MutationalPatterns')
  }
  
  
  # Recalculate using just tumors in current run. Tumors with no non-recurrent SNVs won't be included in output.
  trinuc_breakdown_per_tumor = trinuc_snv_counts(cesa@maf[all_tumors, on = "Unique_Patient_Identifier"],
                                                 genome = get_cesa_bsg(cesa), exclude_recurrent = TRUE,
                                                 style = 'MutationalPatterns')
  
  # can only find signatures in tumors with exome/genome data that have non-recurrent SNVs
  tumors_eligible_for_trinuc_calc = intersect(curr_sample_info[coverage != "targeted", Unique_Patient_Identifier], unique(colnames(trinuc_breakdown_per_tumor)))
  tumors_needing_group_average_rates = setdiff(all_tumors, tumors_eligible_for_trinuc_calc)
  
  substitution_counts = colSums(trinuc_breakdown_per_tumor)
  substitution_counts = substitution_counts[tumors_eligible_for_trinuc_calc]
  
  tumors_above_threshold = names(which(substitution_counts >= sig_averaging_threshold))
  if(length(tumors_above_threshold) == 0) {
    stop(paste0("No tumors have enough mutations to confidently assess mutational signatures. To run anyway, lower\n",
                "sig_averaging_threshold (possibly to 0), and consider setting assume_identical_mutational_processes = TRUE"))
  }
  # identify tumors that will have average weights blended into their weights (normally)
  # when assume_identical_mutational_processes = TRUE, data from these tumors isn't used at all
  tumors_below_threshold = names(which(substitution_counts < sig_averaging_threshold))
  
  tri.counts.genome = get_ref_data(cesa, "tri.counts.genome")
  
  # only used for deconstructSigs since MutationalPatterns doesn't support normalization
  if (signature_extractor == 'deconstructSigs') {
    # for each exome coverage gr (besides default generic, which is pre-calculated), tabulate trinucs
    tri_counts_by_gr = list()
    exome_cov_to_calc = curr_sample_info[coverage == "exome", unique(covered_regions)]
    genome_cov_to_calc = curr_sample_info[coverage == "genome", setdiff(unique(covered_regions), "genome")]
    grs_to_calc = c(cesa@coverage$exome[exome_cov_to_calc], cesa@coverage$genome[genome_cov_to_calc])
    gr_names = names(grs_to_calc)
    for (gr_name in gr_names) {
      if (gr_name %in% c("exome", "exome+")) {
          tri_counts_by_gr[[gr_name]] = get_ref_data(cesa, "tri.counts.exome")
      } else {
        bsg = get_cesa_bsg(cesa)
        covered_seq = BSgenome::getSeq(bsg, grs_to_calc[[gr_name]])
        tri_contexts = Biostrings::trinucleotideFrequency(covered_seq)
        tri_contexts = colSums(tri_contexts)
        
        # deconstructSigs_trinuc_string is internal in cancereffectsizeR
        # here, we need unique trinucleotide contexts (without mutations) in deconstructSigs ordering,
        # which is why we produce this by converting from their column headings
        context_names = unique(sub("\\[([ACTG]).*\\]", "\\1", deconstructSigs_trinuc_string))
        context_names = sort(context_names)
        reverse_complement_names = unique(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(context_names))))
        
        # reorder the counts as desired, then save as a data.frame since that's what deconstructSigs wants
        tri_counts = tri_contexts[context_names] + tri_contexts[reverse_complement_names]
        tri_counts = data.frame(x = tri_counts)
        tri_counts_by_gr[[gr_name]] = tri_counts
      }
    } 
  }

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
  
  signature_names = rownames(signatures)

  # for parallelization, a function to process each tumor
  process_tumor = function(tumor_name) {
    tumor_trinuc_counts = trinuc_breakdown_per_tumor[, tumor_name]
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
    
    if (signature_extractor == 'MutationalPatterns') {
      # MutationalPatterns requires a matrix where columns are samples
      tumor_trinuc_counts = as.matrix(tumor_trinuc_counts)
      all_weights = run_mutational_patterns(tumor_trinuc_counts = tumor_trinuc_counts, signatures_df = signatures, 
                                       signatures_to_remove = current_sigs_to_remove, mp_strict_args = mp_strict_args,
                                       bootstrap_mutations = bootstrap_mutations)
    } else {
      # Set normalization argument for deconstructSigs based on coverage
      covered_regions = curr_sample_info[tumor_name, covered_regions]
      if (covered_regions == "genome") {
        # this actually means no normalization (since signatures and MAF coverage are both whole-genome)
        normalization = "default"
      } else  {
        normalization = tri.counts.genome / tri_counts_by_gr[[covered_regions]]
      }
      
      # deconstructSigs requires a data.frame where the rows are samples
      tumor_trinuc_counts = as.data.frame(t(trinuc_breakdown_per_tumor[,tumor_name]))
      # must provide column names (trinucleotide mutations) and row names (tumors) otherwise
      # deconstructSigs crashes
      colnames(tumor_trinuc_counts) = rownames(trinuc_breakdown_per_tumor)
      rownames(tumor_trinuc_counts) = tumor_name

      all_weights = run_deconstructSigs(tumor_trinuc_counts = tumor_trinuc_counts, tri.counts.method = normalization,
                                        signatures_df = signatures, signatures_to_remove = current_sigs_to_remove)
    }
    return(all_weights)
  }

  
  
  # store results
  # matrix with rows = tumors, columns = relative rates of mutation for each trinuc SNV (starts empty)
  trinuc_proportion_matrix <- matrix(data = NA, nrow = length(all_tumors), ncol = nrow(trinuc_breakdown_per_tumor),
                                     dimnames = list(all_tumors, rownames(trinuc_breakdown_per_tumor)))
  
  
  message("Extracting mutational signatures from SNVs...")
  raw_signature_output = rbindlist(pbapply::pblapply(tumors_eligible_for_trinuc_calc, process_tumor, cl = cores))
  raw_signature_output$Unique_Patient_Identifier = tumors_eligible_for_trinuc_calc

  rel_bio_weights = artifact_account(raw_signature_output, signature_names, artifact_signatures, 
                                              fail_if_zeroed = FALSE)
  trinuc_rates = calculate_trinuc_rates(weights = as.matrix(rel_bio_weights[, -c("Unique_Patient_Identifier")]),
                                        signatures = as.matrix(signatures), tumor_names = tumors_eligible_for_trinuc_calc)
  
  setkey(raw_signature_output, 'Unique_Patient_Identifier')
  setkey(rel_bio_weights, "Unique_Patient_Identifier")
  blended_weights = copy(raw_signature_output)
  
  # Set aside tumors that get all-zero biological signature weights (rare).
  # For the rest, put trinuc rates in output matrix.
  zeroed_out_tumors = rownames(trinuc_rates)[rowSums(trinuc_rates) == 0]
  not_zeroed_tumors = setdiff(rownames(trinuc_rates), zeroed_out_tumors)
  trinuc_proportion_matrix[not_zeroed_tumors, ] = trinuc_rates[not_zeroed_tumors,]

  # in the (very rare) occurrence of zeroed-out tumors, set them aside to be assigned group average signature weightings
  tumors_above_threshold = setdiff(tumors_above_threshold, zeroed_out_tumors)
  tumors_below_threshold = union(tumors_below_threshold, zeroed_out_tumors)
  tumors_needing_group_average_rates = union(tumors_needing_group_average_rates, zeroed_out_tumors)

  if(assume_identical_mutational_processes == TRUE) {
    tumors_needing_group_average_rates = union(tumors_needing_group_average_rates, tumors_eligible_for_trinuc_calc)
  }
  
  # If any samples have subthreshold number of mutations, determine weightings for those samples using method specified by user
  aggregate_extraction = NULL
  if(length(tumors_below_threshold) > 0 || length(tumors_needing_group_average_rates) > 0) {
    if(length(tumors_above_threshold) == 0) {
      stop(paste0("No tumors have enough mutations to inform group-average signature extraction\n", 
                  "(or, nearly impossibly, all such tumors have all their mutations attributed to artifacts)."))
    }
    obs_trinuc_prop = colMeans(trinuc_proportion_matrix[tumors_above_threshold, , drop = F])

    
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
    
    message("Determining group-average signatures from samples with at least ", sig_averaging_threshold, " SNVs...")
    if (signature_extractor == "MutationalPatterns") {
      # supply sample data in columns and as matrix as required by MutationalPatterns
      # Even with the bootstrap_mutations = TRUE method, not going to bootstrap here since tumor_trinuc_counts
      # were already generated off of bootstrapped samples.
      obs_trinuc_prop = as.matrix(obs_trinuc_prop)
      mean_all_weights = run_mutational_patterns(tumor_trinuc_counts = obs_trinuc_prop, signatures_df = signatures, 
                                              signatures_to_remove = signatures_to_remove)
    } else {
      # convert to data.frame and supply rowname as required by deconstructSigs
      obs_trinuc_prop = as.data.frame(t(obs_trinuc_prop))
      rownames(obs_trinuc_prop) = 'mean'
      
      mean_all_weights = run_deconstructSigs(tumor_trinuc_counts = obs_trinuc_prop, signatures_df = signatures, 
                                              signatures_to_remove = signatures_to_remove, tri.counts.method = "default")
    }
    
    mean_all_weights = as.data.table(mean_all_weights)
    mean_all_weights$Unique_Patient_Identifier = 'mean'
    mean_rel_bio_weights = artifact_account(weights = mean_all_weights, signature_names = signature_names,
                                                     artifact_signatures = artifact_signatures, fail_if_zeroed = FALSE)
    mean_all_weights[, Unique_Patient_Identifier := NULL]
    mean_rel_bio_weights[, Unique_Patient_Identifier := NULL]
    
    # calculate rates using just the mean weights for use in assume_identical_mutational_processes (and tumors_needing_group_average_rates)
    # converting to numeric for insertion into trinuc_proportion_matrix
    agg_run_trinuc_prop = as.numeric(calculate_trinuc_rates(weights = as.matrix(mean_rel_bio_weights), signatures = as.matrix(signatures), tumor_names = 'mean'))
    
    mean_all_weights = as.data.frame(mean_all_weights)
    mean_rel_bio_weights = as.data.frame(mean_rel_bio_weights)
    rownames(mean_rel_bio_weights) = rownames(mean_all_weights) = 'mean'
    
    aggregate_extraction = list(rel_bio_weights = mean_rel_bio_weights, all_weights = mean_all_weights)
    
    # TGS tumors, tumors with zeroed-out weights, and tumors with no non-recurrent SNVs get assigned group-average rates
    for (tumor in tumors_needing_group_average_rates) {
      trinuc_proportion_matrix[tumor, ] = agg_run_trinuc_prop
    }
    
    # this should never happen
    if(all(mean_rel_bio_weights == 0)) {
      stop("Somehow, mean signature weights across all samples are all zero.")
    }
    if (assume_identical_mutational_processes) {
      # when assume_identical_mutational_processes is TRUE, this is where trinuc_proportion_matrix gets built
      for(i in 1:nrow(trinuc_proportion_matrix)) {
        trinuc_proportion_matrix[i,] = agg_run_trinuc_prop
      }
    } else {
      for (tumor in tumors_below_threshold) {
        # handle tumors with 0 non-recurrent SNVs
        if(is.na(substitution_counts[tumor])) {
          trinuc_proportion_matrix[tumor, ] = agg_run_trinuc_prop
        } else {
          # zeroed-out tumors are sub-threshold even when sig_averaging_threshold is 0; take no weight from these
          if (sig_averaging_threshold == 0 || tumor %in% zeroed_out_tumors) {
            own_weighting = 0
          } else {
            own_weighting = substitution_counts[tumor] / sig_averaging_threshold
          }
          group_weighting = 1 - own_weighting
          mean_blended_trinuc_prop = trinuc_proportion_matrix[tumor, ] * own_weighting + agg_run_trinuc_prop * group_weighting
          rel_bio_weights[tumor, (signature_names) :=  as.numeric(.SD) * own_weighting + mean_rel_bio_weights * group_weighting, .SDcols = signature_names]
          blended_weights[tumor, (signature_names) := as.numeric(.SD) * own_weighting + mean_all_weights * group_weighting, .SDcols = signature_names]
          trinuc_proportion_matrix[tumor, ] = mean_blended_trinuc_prop
        }
      }  
    }
  }
  
  # If this run is the first setting of trinuc rates, create lists
  if (length(cesa@trinucleotide_mutation_weights) == 0) {
    cesa@trinucleotide_mutation_weights = list(trinuc_snv_counts = full_trinuc_snv_counts,
                                               trinuc_proportion_matrix=trinuc_proportion_matrix)
  } else {
    new_mat = rbind(trinuc_proportion_matrix, cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix)
    cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix = new_mat
    cesa@trinucleotide_mutation_weights$trinuc_snv_counts = full_trinuc_snv_counts
  }
  
  # Put together signature weight tables
  if (! assume_identical_mutational_processes) {
    # get all signature weights into a data table (one row per sample)
    # adjusted_sig_table includes artifact signature_weights
    bio_sig_table = rel_bio_weights
    adjusted_sig_table = blended_weights
    bio_sig_table[, group_avg_blended := Unique_Patient_Identifier %in% tumors_below_threshold] # will include zeroed tumors
    bio_sig_table[, sig_extraction_snvs := as.numeric(substitution_counts[Unique_Patient_Identifier])] # otherwise will be "table" class
    total_snv_counts = cesa@maf[variant_type == "snv"][bio_sig_table, .(total_snvs = .N), on = "Unique_Patient_Identifier", by = "Unique_Patient_Identifier"]
    bio_sig_table = bio_sig_table[total_snv_counts, on = "Unique_Patient_Identifier"]
    
    adjusted_sig_table[bio_sig_table, 
                  c("group_avg_blended", "sig_extraction_snvs", "total_snvs") := list(group_avg_blended, sig_extraction_snvs, total_snvs),
                  on = "Unique_Patient_Identifier"]
    
    # to reflect that no SNVs informed inferred weights of zeroed-out tumors, set sig_extraction_snvs to 0
    # note that adjusted_sig_table keeps the actual counts used by signature extractor
    bio_sig_table[Unique_Patient_Identifier %in% zeroed_out_tumors, sig_extraction_snvs := 0] 
    
    
    tumors_without_data = setdiff(curr_sample_info$Unique_Patient_Identifier, bio_sig_table$Unique_Patient_Identifier)
    num_to_add = length(tumors_without_data)
    if (num_to_add > 0) {
      group_avg_weights = as.numeric(aggregate_extraction$rel_bio_weights)
      new_rows = matrix(nrow = num_to_add, data = rep.int(group_avg_weights, num_to_add), byrow = T)
      colnames(new_rows) = colnames(aggregate_extraction$rel_bio_weights)
      total_snvs = cesa@maf[variant_type == "snv"][, .N, keyby = "Unique_Patient_Identifier"][tumors_without_data, N]
      total_snvs[is.na(total_snvs)] = 0
      new_table = data.table(Unique_Patient_Identifier = tumors_without_data, total_snvs = total_snvs, 
                             sig_extraction_snvs = 0, group_avg_blended = T)
      bio_sig_table = rbind(bio_sig_table, cbind(new_table, new_rows))
      
      # repeat for table that includes artifacts
      group_avg_weights_raw = as.numeric(aggregate_extraction$all_weights)
      new_rows_raw = matrix(nrow = num_to_add, data = rep.int(group_avg_weights_raw, num_to_add), byrow = T)
      colnames(new_rows_raw) = colnames(aggregate_extraction$all_weights)
      new_table[, sig_extraction_snvs := as.numeric(substitution_counts[Unique_Patient_Identifier])]
      new_table[is.na(sig_extraction_snvs), sig_extraction_snvs := 0]
      adjusted_sig_table = rbind(adjusted_sig_table, cbind(new_table, new_rows_raw))
    }
    setcolorder(bio_sig_table, c("Unique_Patient_Identifier", "total_snvs", "sig_extraction_snvs", "group_avg_blended"))
    setcolorder(adjusted_sig_table, c("Unique_Patient_Identifier", "total_snvs", "sig_extraction_snvs", "group_avg_blended"))
    
    new_bio_sig_table = rbind(cesa@trinucleotide_mutation_weights[["signature_weight_table"]], bio_sig_table, fill = T)
    cesa@trinucleotide_mutation_weights[["signature_weight_table"]] = new_bio_sig_table
    new_adjusted_sig_table = rbind(cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]], adjusted_sig_table, fill = T)
    cesa@trinucleotide_mutation_weights[["signature_weight_table_with_artifacts"]] = new_adjusted_sig_table
    
    setcolorder(raw_signature_output, "Unique_Patient_Identifier")
    
    new_raw_table = rbind(cesa@trinucleotide_mutation_weights[["raw_signature_weights"]], raw_signature_output, fill = T)
    setkey(new_raw_table, "Unique_Patient_Identifier")
    cesa@trinucleotide_mutation_weights[["raw_signature_weights"]] = new_raw_table
    
    
  }
  
  if(! is.null(aggregate_extraction)) {
    if (is.null(sample_group)) {
      cesa@trinucleotide_mutation_weights[["group_average_dS_output"]][["all"]] = aggregate_extraction
    } else {
      cesa@trinucleotide_mutation_weights[["group_average_dS_output"]][[paste(sample_group, collapse = ',')]] = aggregate_extraction
    }
  }

  
  curr_grp_num = max(c(0, na.omit(cesa@samples$sig_analysis_grp))) + 1
  cesa@samples[curr_sample_info$Unique_Patient_Identifier, sig_analysis_grp := curr_grp_num, on = "Unique_Patient_Identifier"]
  return(cesa)
}

#' Clear mutational signature attributions and related mutation rate information
#' 
#' Removes all data calculated or supplied via trinuc_mutation_rates,
#' set_signature_weights, set_trinuc_rates, etc. This function can be used if you want to
#' re-run signature analysis with different sample groupings or parameters.
#' 
#' @param cesa CESAnalysis
#' @export
clear_trinuc_rates_and_signatures = function(cesa) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be CESAnalysis")
  }
  if(all(is.na(cesa@samples$sig_analysis_grp))) {
    stop("The analysis has no signature extraction or trinucleotide-context-specific mutation rate data.")
  }
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  cesa@trinucleotide_mutation_weights = list()
  cesa@samples[, sig_analysis_grp := NA_integer_]
  cesa@advanced$snv_signatures = NULL
  return(cesa)
}


#' Tabulate SNVs by trinucleotide context
#' 
#' This function produces trinucleotide-context-specific SNV counts from MAF data for
#' input to mutational signature extraction tools. Output can be tailored to meet
#' formatting requirements of MutationalPatterns or deconstructSigs, which are probably
#' similar to formats used by other tools.
#' 
#' @param maf a cancereffectsizeR-style MAF data table
#' @param genome BSgenome reference genome (for looking up trinucleotide contexts)
#' @param exclude_recurrent Default FALSE. When TRUE, only mutations private to each sample are included in counts, in order to
#' reduce the influence of selection. (If you load more MAF data into the CESAnalysis later, recurrency may change.)
#' @param style "MutationalPatterns" or "deconstructSigs"
#' @return Matrix or data frame of SNV counts, suitable for use with MutationalPatterns or
#'   deconstructSigs. Samples with zero passing SNVs will not appear.
#' @export
#' 
#' 
trinuc_snv_counts = function(maf,
                             genome,
                             exclude_recurrent = FALSE,
                             style = "MutationalPatterns"
                             ) {
  if(! is(maf, "data.table")) {
    "maf should be an MAF-like data.table."
  }
  
  if(maf[, .N] == 0) {
    stop("No MAF data in the CESAnalysis", call. = F)
  }
  
  if(! all(c("variant_type", "Unique_Patient_Identifier") %in% names(maf))) {
    stop("Expected columns Unique_Patient_Identifier and variant_type in maf, as seen in MAF tables validated by preload_maf().")
  }
  
  bsg = genome
  if (! is(bsg, "BSgenome")) {
    stop("genome should be a BSgenome object.")
  }
  
  if(! is.logical(exclude_recurrent) || length(exclude_recurrent) != 1) {
    stop("exclude_recurrent should be TRUE/FALSE.")
  }
  
  if (! is(style, "character") || length(style) != 1 || 
      ! tolower(style) %in% c('mutationalpatterns', 'deconstructsigs')) {
    stop("style must be 'MutationalPatterns' or 'deconstructSigs'.")
  }
  style = ifelse(tolower(style) == "mutationalpatterns", 'MutationalPatterns', 'deconstructSigs')
  
  # All data in MAF (even if it's TGS, etc.) impacts recurrency testing
  ds_maf = maf[variant_type == "snv"]
  
  if("problem" %in% names(ds_maf)) {
    problems = ds_maf[! is.na(problem), which = T]
    if (length(problems) > 0) {
      ds_maf = ds_maf[! problems]
      message("Found a preload_maf()-style \"problem\" column and removed ", length(problems), " 46 problematic records.")
    }
  }
  if (exclude_recurrent) {
    # remove all recurrent SNVs (SNVs appearing in more than one sample)
    duplicated_vec_first <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)])
    duplicated_vec_last <- duplicated(ds_maf[,.(Chromosome, Start_Position, Tumor_Allele)],fromLast=T)
    duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
    if (length(duplicated_vec_pos) > 0) {
      ds_maf <- ds_maf[-duplicated_vec_pos,]
    }
  }
  
  # build the data.frame required by MutationalPatterns (similar to deconstructSigs)
  # rows are samples, columns are counts of each of the 96 trinuc-context-specific mutations in the 
  # order expected by deconstructSigs
  trinuc = BSgenome::getSeq(bsg, ds_maf$Chromosome, ds_maf$Start_Position - 1, ds_maf$Start_Position + 1, as.character = T)
  
  # internal table converts trinuc/mut (e.g., GTA:C) into deconstructSigs format ("G[T>C]A")
  trinuc_snv = factor(deconstructSigs_notations[.(trinuc, ds_maf$Tumor_Allele), deconstructSigs_ID], levels = deconstructSigs_trinuc_string)
  
  # mysteriously convert two-way table to data.frame
  tmp = table(ds_maf$Unique_Patient_Identifier, trinuc_snv)
  counts = apply(tmp, 2, rbind)
  
  # edge case: when just 1 sample in data set, counts comes back as integer vector but need matrix
  if(! is.matrix(counts)) {
    counts = t(as.matrix(counts))
  }

  rownames(counts) = rownames(tmp)
  
  if (style == 'MutationalPatterns') {
    return(t(counts))
  } else {
    return(as.data.frame(counts))
  }
}




