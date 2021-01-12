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
#' @param more_cols Names of data columns to keep in addition to the core required
#'   columns. Choose "all" to keep all columns.
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
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", more_cols = NULL) {
  
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
           "remotes::install_github(\"Townsend-Lab-Yale/ces-reference-data/", refset_name, "\")")
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
    print(problem_summary, row.names = F)
    pretty_message("You can remove or fix these records, or let load_maf() exclude them automatically.")
  }
  return(maf[])
}



# load_maf = function(cesa = NULL, maf = NULL, annotate = TRUE, , coverage = "exome", coverage_intervals_to_check = NULL,
#                     coverage_intervals_to_check_name = NULL, coverage_intervals_to_check_padding = 0, group_col = NULL, chain_file = NULL, enforce_default_exome_coverage = FALSE) {
#   
#   
#   
#   
#CoveredRegions 