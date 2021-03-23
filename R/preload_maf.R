#' Read and verify MAF somatic mutation data 
#' 
#' Reads MAF-formatted data from a text file or data table, checks for problems, and
#' provides a few quality check annotations (if available). If core MAF columns don't have
#' standard names (Chromosome, Start_Position, etc., with Tumor_Sample_Barcode used as the
#' sample ID column), you can supply your own column names. If the data you are loading is
#' from a different genome build than the chosen reference data set (refset) you can use
#' the "chain_file" option to supply a UCSC-style chain file, and your MAF coordinates
#' will be automatically converted with rtracklayer's version of liftOver.
#' 
#' The default ces.refset.hg19 provides three annotations that may consider using for
#' quality filtering of MAF records:
#' \itemize{
#' \item cosmic_site_tier Indicates if the variant's position overlaps a mutation in
#' COSMIC v92's Cancer Mutation Census. Mutations are classified as Tier 1, Tier 2, Tier
#' 3, and Other. Note that the MAF mutation itself is not necessarily in the census. See
#' COSMIC's website for tier definitions.
#' \item germline_variant_site The variant's position overlaps a site of common germline
#' variation. Roughly, this means that gnomAD 2.1.1 shows an overlapping germline variant at
#' greater than 1% prevalence in some population.
#' \item repetitive_region The variant overlaps a site marked as repetitive sequence by
#' the RepeatMasker tool (data taken from UCSC Table Browser). Variant calls in repetitive
#' sites frequently reflect calling error.
#' }
#' 
#' @param maf Path of tab-delimited text file in MAF format, or a data.table/data.frame with MAF data
#' @param refset name of reference data set (refset) to use; run \code{list_ces_refsets()} for
#'   available refsets. Alternatively, the path to a custom reference data directory.
#' @param sample_col column name with sample ID data (Tumor_Sample_Barcode or Unique_Patient_Identifier)
#' @param chr_col column name with chromosome data  (Chromosome)           
#' @param start_col column name with start position (Start_Position)
#' @param ref_col column name with reference allele data (Reference_Allele)
#' @param tumor_allele_col column name with alternate allele data; by default,
#'   values from Tumor_Seq_Allele2 and Tumor_Seq_Allele1 columns are used
#' @param keep_extra_columns TRUE/FALSE to load data columns not needed by cancereffectsizeR,
#' or a vector of column names to keep.
#' @param coverage_intervals_to_check If available, a BED file or GRanges object
#'   represented the expected coverage intervals of the sequencing method used to generate
#'   the MAF data. Unless the coverage intervals are incorrect, most records will be
#'   covered. Output will show how far away uncovered records are from covered regions,
#'   which can inform whether to use the covered_regions_padding option in load_maf().
#'   (For example, some variant callers will identify variants up to 100bp out of the
#'   target regions, and you may want to pad the covered intervals to allow these variants
#'   to remain in your data. Alternatively, if all records are already covered, then the
#'   calls have probably already be trimmed to the coverage intervals, which means no
#'   padding should be added.)
#' @return a data.table of MAF data, with any problematic records flagged and a few
#'   quality-control annotations (if available with the chosen refset data).
#' @export
preload_maf = function(maf = NULL, refset = "ces.refset.hg19", coverage_intervals_to_check = NULL,
                    chain_file = NULL, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", keep_extra_columns = FALSE) {
  
  if (is(refset, "environment")) {
    refset_name = as.character(substitute(refset))
  } else {
    refset_name = refset
  }
  # Check for and load reference data for the chosen genome/transcriptome data
  if (! is(refset_name, "character")) {
    stop("refset should be a refset object, the name of an installed refset package, or a path to custom refset directory.")
  }
  using_custom_refset = TRUE
  if (refset_name %in% names(.official_refsets)) {
    using_custom_refset = FALSE
    if(file.exists(refset_name)) {
      stop("You've given the name of a CES reference data set package, but a file/folder with the same name is in your working directory. Stopping to avoid confusion.")
    }
    if(! require(refset_name, character.only = T)) {
      if(refset_name == "ces.refset.hg19") {
        message("Install ces.refset.hg19 like this:\n",
                "remotes::install_github(\"Townsend-Lab-Yale/ces.refset.hg19@*release\")")
      }
      stop("CES reference data set ", refset_name, " not installed.")
    }
    req_version = .official_refsets[[refset_name]]
    actual_version = packageVersion(refset_name)
    if (actual_version < req_version) {
      stop("CES reference data set ", refset_name, " is version ", actual_version, ", but your version of cancereffectsizeR requires at least ",
           "version ", req_version, ".\nRun this to update:\n",
           "remotes::install_github(\"Townsend-Lab-Yale/", refset_name, "\")")
    }
    ref_data_version = actual_version
    data_dir = system.file("refset", package = refset_name)
  } else {
    if (! dir.exists(refset_name)) {
      if (grepl('/', refset_name)) {
        stop("Could not find reference data at ", refset_name)
      } else {
        stop("Invalid reference set name. Check spelling, or view available data sets with list_ces_refsets().")
      }
    }
    
    data_dir = refset_name
    refset_name = basename(refset_name)
    if(refset_name %in% names(.official_refsets)) {
      stop("Your custom reference data set has the same name (", refset_name, ") as a CES reference data package. Please rename it.")
    }
  }
  
  if (refset_name %in% ls(.ces_ref_data)) {
    refset_env = .ces_ref_data[[refset_name]]
  } else if(using_custom_refset) {
    refset_env = preload_ref_data(data_dir)
    .ces_ref_data[[refset_name]] = refset_env
  } else {
    refset_env = get(refset_name, envir = as.environment(paste0('package:', refset_name)))
    .ces_ref_data[[refset_name]] = refset_env
  }
  
  # By default, only load core columns to save time and memory
  more_cols = NULL
  if (identical(keep_extra_columns, TRUE)) {
    more_cols = 'all'
  } else if (! identical(keep_extra_columns, FALSE)) {
    if (is.character(keep_extra_columns)) {
      more_cols = keep_extra_columns
    } else {
      stop("keep_extra_columns should be T/F or names of extra columns to include.")
    }
  }
  maf = read_in_maf(maf = maf, refset_env = refset_env, chr_col = chr_col, start_col = start_col, ref_col = ref_col,
                    tumor_allele_col = tumor_allele_col, sample_col = sample_col, more_cols = more_cols, chain_file = chain_file)
  
  coverage_gr = NULL
  if(! is.null(coverage_intervals_to_check)) {
    if (is.character(coverage_intervals_to_check)) {
      if (length(coverage_intervals_to_check) != 1) {
        stop("coverage_intervals_to_check should be a BED filename or GRanges object")
      }
      if (! file.exists(coverage_intervals_to_check)) {
        stop("BED file not found; check path?", call. = F)
      }
      coverage_gr = rtracklayer::import.bed(coverage_intervals_to_check)
    } else if (is(coverage_intervals_to_check, 'GRanges')) {
      coverage_gr = coverage_intervals_to_check
    } else {
      stop("coverage_intervals_to_check should be a BED filename or GRanges object")
    }
    coverage_gr = clean_granges_for_cesa(refset_env = refset_env, gr = coverage_gr)
  }
  
  if ('preload_anno' %in% ls(.ces_ref_data[[refset_name]])) {
    anno_grs = .ces_ref_data[[refset_name]][['preload_anno']]
  } else {
    # If there are no preload annotation grs already loaded, check if any exist and load
    # them, or return if none exist
    preload_anno_files = list.files(paste0(data_dir, "/maf_preload_anno"), pattern = '\\.rds$', full.names = T)
    anno_grs = lapply(preload_anno_files, readRDS)
    if(length(anno_grs) > 0) {
      .ces_ref_data[[refset_name]][['preload_anno']] = anno_grs
    }
  }
  
  
  # Make MAF-based gr if needed
  which_no_problem = maf[is.na(problem), which = T]
  if (! is.null(coverage_gr) || length(anno_grs) > 0) {
    maf_gr = makeGRangesFromDataFrame(maf[which_no_problem], start.field = "Start_Position", end.field = "Start_Position",
                             seqnames.field = "Chromosome")
  }
  
  
  if(! is.null(coverage_gr)) {
    dist = as.data.table(distanceToNearest(maf_gr, coverage_gr))
    maf[which_no_problem, dist_to_coverage_intervals := dist$distance]
  }


  for (gr in anno_grs) {
    # can either be a GRanges or a list of them
    if (is(gr, "GRanges")) {
      anno_colname = attr(gr, "anno_col_name", exact = T)
      maf[which_no_problem, (anno_colname) := maf_gr %within% gr]
    } else if (is(gr, "list") && unique(sapply(gr, function(x) is(x, "GRanges"))) == T) {
      # will go in reverse order since ranges are listed in order of precedence (first gr's overlaps should always appear)
      anno_colname = attr(gr, "anno_col_name", exact = T)
      gr = rev(gr)
      labels = names(gr)
      for (i in 1:length(labels)) {
        curr_label = labels[i]
        has_overlap = maf_gr %within% gr[[i]]
        maf[which_no_problem[has_overlap], (anno_colname) := curr_label]
      }
    } else {
      warning("A misformatted annotation source was skipped.")
    }
  }
  problem_summary = maf[! which_no_problem, .(num_records = .N), by = "problem"]
  
  if(problem_summary[, .N] > 0) {
    pretty_message("Some MAF records have problems:")
    # this is how to print a table nicely in the message stream
    message(crayon::black(paste0(capture.output(print(problem_summary, row.names = F)), collapse = "\n"))) 
    pretty_message("You can remove or fix these records, or let load_maf() exclude them automatically.")
  }
  return(maf[])
}


#' Catch duplicate samples
#' 
#' Takes in a data.table of MAF data (produced, typically, with \code{preload_maf()}) and
#' identifies samples with relatively high proportions of shared SNV mutations. Some
#' flagged sample pairs may reflect shared driver mutations or chance overlap of variants
#' in SNV or sequencing error hotspots. Very high overlap may indicate sample duplication,
#' re-use of samples across data sources, or within-experiment sample contamination. To limit
#' the influence of shared calling error, it's recommended to run this function after
#' any quality filtering of MAF records, as a final step.
#' 
#' Sample pairs are flagged when...
#' \itemize{
#' \item Both samples have <6 total SNVs and any shared SNVs.
#' \item Both samples have <21 total SNVs and >1 shared mutation.
#' \item One sample has just 1 or 2 total SNVs and has any overlaps with the other sample.
#' \item The samples have >2 shared SNVs and at least one percent of SNVs are shared (in the sample with fewer SNVs).
#' }
#' In addition to reporting SNV counts, this function divides the genome into 1000-bp
#' windows and reports the following:
#' \itemize{
#' \item variant_windows_A: Number of windows in which sample A has a variant.
#' \item variant_windows_B: Same for B.
#' \item windows_shared: Number of windows that contain a variant shared between both samples.
#' }
#' Sometimes, samples have little overlap except for a few hotspots that may derive from
#' shared calling error or highly mutable regions. These window counts can help
#' distinguish such samples from those with more pervasive SNV overlap.
#' 
#' @param maf_table A data.table with MAF data and cancereffectsizeR-style column names,
#'   as generated by \code{preload_maf()}.
#' @return a data.table with overlap statistics
#' @export
check_sample_overlap = function(maf_table) {
  maf = maf_table
  maf_cols = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele", "Unique_Patient_Identifier")
  if(! all(maf_cols %in% names(maf))) {
    stop("Missing some required MAF columns (consider running your MAF data through preload_maf() first.")
  }
  if("problem" %in% names(maf) & ! all(is.na(maf$problem))) {
    message("Excluding MAF records with problems from this check (if problems were fixed, stop and re-run preload_maf)...")
    maf = maf[is.na(problem)]
  }
  maf = maf[, ..maf_cols]
  nt = c('A', 'C', 'T', 'G')
  maf = maf[Reference_Allele %in% nt & Tumor_Allele %in% nt]
  maf[, mut_id := paste(Chromosome, Start_Position, Tumor_Allele, sep = "_")]
  
  
  # For efficiency, only need to check recurrent mutations
  maf[, window := paste(Chromosome, Start_Position %% 1000, sep = '.')]
  recurrent_muts = maf[, .N, by = "mut_id"][N > 1, mut_id]
  records_to_check = maf[recurrent_muts, on = "mut_id"]
  
  # Also will count shared mutations by 1000bp window
  setkey(records_to_check, "Unique_Patient_Identifier")
  setindex(records_to_check, "mut_id")
  
  # Only need to count samples with recurrent mutations
  samples_to_check = maf[recurrent_muts, unique(Unique_Patient_Identifier), on = 'mut_id']
  
  # Function to handle overlap counting
  count_overlaps = function(samples_to_check, dt) {
    sample_pairs = utils::combn(samples_to_check, 2, simplify = F)
    
    # Organize mutations by sample
    mutations <- new.env(hash=TRUE)
    for (sample in samples_to_check) {
      mutations[[sample]] = dt[sample, mut_id]
    }
    
    ## Assess overlap of mutations between every pair of samples
    count = 0
    with_inter = 0
    num_samples = length(samples_to_check)
    possible_dups = list()
    i = 0
    for (pair in sample_pairs) {
      s1_muts = mutations[[pair[1]]]
      s2_muts = mutations[[pair[2]]]
      muts_in_both = intersect(s1_muts, s2_muts)
      intersection = length(muts_in_both)
      if (intersection > 0) {
        with_inter = with_inter + 1
        intersect_windows = uniqueN(dt[muts_in_both, window, on = "mut_id"])
        possible_dups[[with_inter]] = c(pair[1], pair[2], intersection, intersect_windows)
      }
      i = i + 1
      if (i %% 10000 == 0) {
        message("Finished ", format(i, big.mark = ',', scientific = F), " pairs.")
      }
    }
    
    possible_dups = rbindlist(lapply(possible_dups, as.list))
    colnames(possible_dups) = c("sample_A", "sample_B", "variants_shared", "windows_shared")
    possible_dups[, variants_shared := as.integer(variants_shared)]
    possible_dups[, windows_shared := as.integer(windows_shared)]
    return(possible_dups)
  }
  
  # Run overlap count function
  num_pairs = format(choose(length(samples_to_check), 2), big.mark = ",", scientific = F)
  pretty_message(paste0("Comparing ", num_pairs, " pairs of samples for suspicious variant overlap..."))
  dups = count_overlaps(samples_to_check = samples_to_check, dt = records_to_check)
  
  # Add total variant counts
  counts_by_sample = maf[, .N, by = "Unique_Patient_Identifier"]
  dups[counts_by_sample, variants_A := N, on = c(sample_A = "Unique_Patient_Identifier")]
  dups[counts_by_sample, variants_B := N, on = c(sample_B = "Unique_Patient_Identifier")]
  
  
  ## Flag cases where both samples have <6 total mutations and any shared mutations
  ## Flag cases where both samples have <21 total mutations and >1 shared mutation
  ## Flag cases where one sample has just 1 or 2 mutations and has any overlaps with other samples
  ## Flag cases where there are >2 shared mutations and at least 1% overlap
  dups[variants_A < 21 & variants_B < 21 & variants_shared > 1, to_examine := T]
  dups[variants_A < 6 & variants_B < 6 & variants_shared > 0, to_examine := T]
  dups[(variants_A < 3 | variants_B < 3) & variants_shared > 0, to_examine := T] 
  
  # Mark pairs with high overlap
  dups[, greater_overlap := pmax((variants_shared / variants_A), (variants_shared / variants_B))]
  dups[(variants_shared > 2 & greater_overlap > .01), to_examine := T]
  
  # Filter down to the pairs to examine manually, add variant window counts, format for output
  to_examine = dups[to_examine == T, -"to_examine"][order(-greater_overlap)]
  to_examine[, greater_overlap := round(greater_overlap, 3)]
  window_counts_by_sample = maf[, .(N = uniqueN(window)), by = "Unique_Patient_Identifier"]
  to_examine[window_counts_by_sample, variant_windows_A := N, on = c(sample_A = "Unique_Patient_Identifier")]
  to_examine[window_counts_by_sample, variant_windows_B := N, on = c(sample_B = "Unique_Patient_Identifier")]
  setcolorder(to_examine, c("sample_A", "variants_A", "variant_windows_A", 
                            "sample_B", "variants_B", "variant_windows_B", "windows_shared",
                            "variants_shared", "greater_overlap"))
  if(to_examine[, .N] == 0) {
    message("No sample pairs met overlap reporting thresholds!")
  } else {
    message("Returning ", to_examine[, .N], " sample pairs for review.")
  }
  return(to_examine)
}
