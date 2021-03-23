#' Load MAF somatic mutation data
#' 
#' Load MAF data from a text file or data table into your CESAnalysis. If column names
#' don't match MAF format specifications (Chromosome, Start_Position, etc., with
#' Tumor_Sample_Barcode used as the sample ID column), you can supply your own column
#' names. When your CESAnalysis has defined sample groups (see \code{?CESAnalysis}),
#' specify "group_col". By default, data is assumed to be derived from whole-exome
#' sequencing. Whole-genome data and targeted sequencing data are also supported when the
#' "coverage" option is specified. If the data you are loading is from a different genome
#' build than your CESAnalysis, you can use the "chain_file" option to supply a UCSC-style
#' chain file, and your MAF coordinates will be automatically converted with
#' rtracklayer's version of liftOver.
#' 
#' @param cesa the CESAnalysis object to load the data into
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.table or data.frame format
#' @param sample_col column name with sample ID data (Tumor_Sample_Barcode or Unique_Patient_Identifier)
#' @param chr_col column name with chromosome data  (Chromosome)           
#' @param start_col column name with start position (Start_Position)
#' @param ref_col column name with reference allele data (Reference_Allele)
#' @param tumor_allele_col column name with alternate allele data; by default,
#'   values from Tumor_Seq_Allele2 and Tumor_Seq_Allele1 columns are used.
#' @param group_col column in MAF with sample group labels (see \code{?CESAnalysis})
#' @param coverage exome, genome, or targeted (default exome)
#' @param covered_regions optional for exome, required for targeted: a GRanges object or a
#'   BED file of covered intervals matching the CESAnalysis genome
#' @param covered_regions_name a name describing the covered regions (e.g.,
#'   "my_custom_targeted_regions"); required when covered_regions are supplied
#' @param covered_regions_padding How many bases (default 0) to expand start and end of
#'   each covered_regions interval, to include variants called just outside of targeted
#'   regions. Consider setting from 0-100bp, or up to the sequencing read length. If the
#'   input data has been trimmed to the targeted regions, leave set to 0.
#' @param chain_file a LiftOver chain file (text format, name ends in .chain) to convert MAF
#'   records to the genome build used in the CESAnalysis.
#' @param enforce_default_exome_coverage When loading default exome data, exclude records
#'   that aren't covered in the default exome capture intervals included with
#'   CES genome reference data (default FALSE).
#' @return CESAnalysis with the specified MAF data loaded
#' @export
load_maf = function(cesa = NULL, maf = NULL, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", coverage = "exome", covered_regions = NULL,
                    covered_regions_name = NULL, covered_regions_padding = 0, group_col = NULL, chain_file = NULL, enforce_default_exome_coverage = FALSE) {
  
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  
  if(cesa@advanced$locked) {
    stop("You can't load more MAF data since you've already calculated some mutation rates. Create a new one if necessary.", call. = F)
  }
  cesa = update_cesa_history(cesa, match.call())
  bsg = get_cesa_bsg(cesa)
  
  # validate chain_file (presence means liftOver must run)
  need_to_liftOver = FALSE
  if(! is.null(chain_file)) {
    need_to_liftOver = TRUE
    if(!is(chain_file, "character") || length(chain_file) != 1 || !(endsWith(chain_file, ".chain"))) {
      stop("Argument chain_file expected to be the filename/path of a text file ending in .chain")
    }
    if(! file.exists(chain_file)) {
      stop("The chain_file specified could not be found; check the file path.")
    }
  }
  
  if (is.null(maf)) {
    stop("Supply MAF data via maf=[file path or data.table/data.frame].")
  }
  
  if (is.null(group_col) & length(cesa@groups) != 1) {
    stop("You must specify group_col in the MAF since this CESAnalysis specifies sample groups.")
  }
  
  if (! is.null(group_col) & length(cesa@groups) == 1) {
    stop(paste0("This CESAnalysis does not specify sample groups, so you can't use the \"group_col\" argument.\n",
                "Create a new CESAnalysis with \"sample_groups\" specified to include this information."))
  }
  
  # give a warning if interval padding is really high
  if (covered_regions_padding > 1000) {
    warning(sprintf("covered_regions_padding=%d is awfully high!", covered_regions_padding), call. = F)
  }
  
  # Validate covered_regions
  previous_covered_regions_names = c(names(cesa@coverage$exome), names(cesa@coverage$targeted)) # may be NULL
  if (! is.character(coverage) || ! coverage %in% c("exome", "genome", "targeted") || length(coverage) > 1) {
    stop("Argument coverage must be \"exome\", \"genome\", or \"targeted\"")
  }
  
  if (! coverage %in% names(cesa@coverage)) {
    cesa@coverage[[coverage]] = list()
  }
  
  
  # validate covered_regions_name (supplied if and only if covered_regions is)
  if (! is.null(covered_regions) & is.null(covered_regions_name)) {
    stop("You must supply a name for your covered_regions using covered_regions_name = ...")
  }
  if (is.null(covered_regions) && ! is.null(covered_regions_name)) {
    stop("covered_regions_name was supplied, but covered_regions wasn't.")
  }
  if (! is.null(covered_regions_name)) {
    if (! is(covered_regions_name, "character") || length(covered_regions_name) > 1) {
      stop("covered_regions_name should be a 1-length character vector, just a name to use for your covered_regions")
    }
  }
  
  if(is.null(covered_regions) && (length(covered_regions_padding) != 1 || covered_regions_padding != 0)) {
    stop("You can't use covered_regions_padding without supplying covered_regions.", call. = F)
  }
  
  # validate covered_regions (required for targeted, optional for exome/genome)
  if (coverage == "genome" & is.null(covered_regions)) {
    covered_regions_name = "genome"
  }
  
  if (coverage == "targeted" & is.null(covered_regions) ) {
    stop("can't load targeted data without covered_regions (see docs)")
  }
  
  if (! is(enforce_default_exome_coverage, "logical") || length(enforce_default_exome_coverage) != 1) {
    stop("enforce_default_exome_coverage should be T/F")
  }
  
  if (coverage == "exome" & is.null(covered_regions)) {
    covered_regions_name = "exome"
    if(! check_for_ref_data(cesa, "default_exome_gr")) {
      stop("This genome has no default exome intervals, so to load exome data you must supply covered_regions (see docs)")
    }
  }
  
  # Determine whether we're using "exome" or "exome+" for default exome data
  # Usually, the lenient exome+ option is used, but the choice must be consistent throughout a CESAnalysis
  if(covered_regions_name == "exome") {
    # Whether or not CESAnalysis uses "exome+", once the first default exome data has been loaded,
    # there will always be default exome intervals under "exome" in the coverage list
    if (! "exome" %in% names(cesa@coverage[["exome"]])) {
      cesa@advanced$using_exome_plus = ! enforce_default_exome_coverage
      covered_regions = get_ref_data(cesa, "default_exome_gr")
      cesa@coverage$exome[["exome"]] = covered_regions
      
    } else if (enforce_default_exome_coverage == TRUE & cesa@advanced$using_exome_plus) {
      stop("You can't have enforce_default_exome_coverage = TRUE because you previously loaded default exome data\n",
           "with enforce_default_exome_coverage = FALSE")
    } else if (enforce_default_exome_coverage == FALSE & ! cesa@advanced$using_exome_plus) {
      stop("You can't have enforce_default_exome_coverage == FALSE because you previously loaded default exome data\n",
           "with enforce_default_exome_coverage = TRUE.")
    } else {
      covered_regions = cesa@coverage$exome[["exome"]]
    }
    
    # For exome+, use GRanges if available
    # Otherwise, for exome or exome+, load the default exome intervals
    if (! enforce_default_exome_coverage) {
      covered_regions_name = "exome+"
      if("exome+" %in% names(cesa@coverage$exome)) {
        covered_regions = cesa@coverage$exome[["exome+"]]
      } else {
        cesa@coverage[["exome"]][["exome+"]] = covered_regions # that is, the default exome regions
      }
    } else {
      if("exome" %in% names(cesa@coverage$exome)) {
        covered_regions = cesa@coverage$exome[["exome"]]
      } else {
        covered_regions = get_ref_data(cesa, "default_exome_gr")
        cesa@coverage$exome[["exome"]] = covered_regions
      }
    }
    pretty_message("Assuming this data has default exome coverage (it's better to supply covered intervals if you have them; see docs)...")
  } else if (! is.null(covered_regions)) {
    # Use internal function to avoid updating covered_in now (annotate_variants step will handle this larger)
    cesa = .add_covered_regions(cesa = cesa, covered_regions = covered_regions, covered_regions_padding = covered_regions_padding, 
                               coverage_type = coverage, covered_regions_name = covered_regions_name, update_anno = FALSE)
  }
  
  refset_env = .ces_ref_data[[cesa@ref_key]]
  read_args = list(maf = maf, refset_env = refset_env,
                   sample_col = sample_col, chr_col = chr_col, start_col = start_col,
                   ref_col = ref_col, tumor_allele_col = tumor_allele_col,
                   chain_file = chain_file)
  if (! is.null(group_col)) {
    read_args = c(read_args, list(more_cols = group_col))
  }
  maf = do.call(read_in_maf, args = read_args)
  

  # Set aside records with problems and notify user
  initial_num_records = maf[, .N]
  excluded = maf[! is.na(problem), .(Unique_Patient_Identifier, Chromosome, Start_Position, 
                                     Reference_Allele, Tumor_Allele, problem)]
  maf = maf[is.na(problem), -"problem"]
  
  nt = c("A", "T", "C", "G")
  maf[Reference_Allele %in% nt & Tumor_Allele %in% nt, variant_type := "snv"]
  maf[is.na(variant_type), variant_type := "indel"]

  num_excluded = excluded[, .N]
  if(num_excluded > 0) {
    msg = paste0(num_excluded, " of ", initial_num_records, " MAF records (", 
                 sprintf("%.1f", 100 * num_excluded / initial_num_records), '%) ',
                 "had problems and were excluded: ")
    problem_summary = excluded[, .(num_records = .N), by = "problem"]
    message(crayon::black(paste0(capture.output(print(problem_summary, row.names = F)), collapse = "\n")))
    
    if(num_excluded / initial_num_records > .05) {
      warning("More than 5% of input records had problems.")
    }
    
    snv_mismatch = excluded[problem == "reference_mismatch" & Reference_Allele %in% nt & Tumor_Allele %in% nt, .N]
    if(snv_mismatch > 0) {
      msg = paste0(snv_mismatch, " SNV variants were excluded for having reference alleles that do not match the reference genome. You should probably figure out why",
                          " and make sure that the rest of your data set is okay to use before continuing.")
      warning(pretty_message(msg, emit = F))
    }
  }
  setnames(excluded, "problem", "Exclusion_Reason")
  
  # collect sample group information
  if (! is.null(group_col)) {
    if (is.factor(maf[[group_col]])) {
      warning("You supplied a sample group column as a factor, but it was converted to character.")
    }
    sample_groups = as.character(maf[[group_col]])
    if(any(is.na(sample_groups))) {
      stop("Error: There are NA values in your sample groups column.")
    }
    maf[, (group_col) := NULL]
  } else {
    sample_groups = cesa@groups[1] # indicates a stageless analysis
  }
  
  new_samples = data.table(Unique_Patient_Identifier = maf$Unique_Patient_Identifier, group = sample_groups)
  new_samples = new_samples[, .(group = unique(group)), by = "Unique_Patient_Identifier"]

  new_samples[, coverage := coverage]
  new_samples[, covered_regions := covered_regions_name]
  
  # ensure no sample has an illegal group
  bad_groups = setdiff(new_samples[, unique(group)], cesa@groups)
  if(length(bad_groups) > 0) {
    stop(paste0("The following groups were not declared in your CESAnalysis, but they were found in your MAF groups column:\n",
                paste(bad_groups, collapse = ", ")))
  }
  # see if any sample appears more than once in sample table (happens when one sample has multiple listed groups)
  repeated_samples = new_samples[duplicated(Unique_Patient_Identifier), unique(Unique_Patient_Identifier)]
  if(length(repeated_samples) > 0) {
    stop(paste0("The following samples are associated with multiple groups in the input data:\n", paste(repeated_samples, collapse=", ")))
  }
  
  # make sure no new samples were already in the previously loaded MAF data
  if(cesa@samples[, .N] > 0) {
    repeat_samples = intersect(cesa@samples[, Unique_Patient_Identifier], new_samples[, Unique_Patient_Identifier])
    if (length(repeat_samples) > 0) {
      stop(paste0("Error: Can't load MAF data because some sample IDs already appear in previously loaded data.\n",
                  "Either merge these data sets manually or remove duplicated samples: ",
                  paste(repeat_samples, collapse = ", ")))
    }
  }
  
  # notify the user if some of the declared groups don't appear in the data at all
  if(length(cesa@groups) > 1) {
    missing_groups = cesa@groups[! cesa@groups %in% new_samples[,unique(group)]]
    if (length(missing_groups) > 0) {
      msg = paste0("The following groups were declared in your CESAnalysis, but they weren't present in the MAF data: \n",
                     paste(missing_groups, collapse = ", "))
      pretty_message(msg)
    }    
  }

  # remove any MAF records that are not in the coverage, unless default exome with enforce_default_exome_coverage = FALSE
  if (covered_regions_name == "genome") {
    num_uncovered = 0
  } else {
    maf_grange = GenomicRanges::makeGRangesFromDataFrame(maf, seqnames.field = "Chromosome", start.field = "Start_Position", 
                                                         end.field = "Start_Position")
    
    # In an "exome+" data set, compare coverage to default exome instead of whatever the current exome+ intervals are
    if (covered_regions_name == "exome+") {
      # equivalent to %within%, but avoids importing
      is_uncovered = ! IRanges::overlapsAny(maf_grange, cesa@coverage[["exome"]][["exome"]], type = "within")
    } else {
      is_uncovered = ! IRanges::overlapsAny(maf_grange, cesa@coverage[[coverage]][[covered_regions_name]], type = "within")
    }
    num_uncovered = sum(is_uncovered)
  }

  if (num_uncovered > 0) {
    total = nrow(maf)
    percent = round((num_uncovered / total) * 100, 1)
    if (covered_regions_name == "exome+") {
      # merge previous exome+ covered_regions with the coverage of the new data
      covered_regions =  GenomicRanges::reduce(GenomicRanges::union(covered_regions, maf_grange[is_uncovered]))
      cesa@coverage[[coverage]][["exome+"]] = covered_regions

      # warn if a lot of records are in uncovered areas; may indicate whole-genome data or low quality exome data
      if (percent > 10) {
        warning(paste0("More than 10% of MAF records are not within the CES genome's default exome intervals.\n",
                       "Could this be whole-genome data? Or if you know the true covered regions, supply them\n",
                       "with the covered_regions argument."))
      }
      msg = paste0("Note: ", num_uncovered, " MAF records (", percent, 
                   "%) are at loci not covered in the default exome intervals ",
                   "in this CESAnalysis's reference data. The samples in this MAF (along with along other MAFs ",
                   "loaded as default exome data), will be assumed to all share the same coverage, ",
                   "which we'll call \"exome+\".")
      pretty_message(msg)
    } else {
      uncovered.maf = maf[is_uncovered,]
      uncovered.maf$Exclusion_Reason = paste0("uncovered_in_", covered_regions_name)
      maf = maf[!is_uncovered,]
      excluded = rbind(excluded, uncovered.maf)
      message(paste0("Note: ", num_uncovered, " MAF records out of ", total, " (", percent, 
                     "%) are at loci not covered in the input covered_regions. ",
                     "\nThese mutations will be excluded from analysis."))
    }
  }
  
  # drop any samples that had all mutations excluded
  new_samples = new_samples[Unique_Patient_Identifier %in% maf$Unique_Patient_Identifier]
  cesa@samples = rbind(cesa@samples, new_samples)
  setcolorder(cesa@samples, c("Unique_Patient_Identifier", "coverage", "covered_regions", "group"))
  setkey(cesa@samples, "Unique_Patient_Identifier")
  
  if (nrow(excluded) > 0) {
    colnames(excluded) = c(colnames(maf)[1:5], "Exclusion_Reason")
    cesa@excluded = rbind(cesa@excluded, excluded) 
  }
  
  cesa@maf = rbind(cesa@maf, maf, fill = T)
  message("Annotating variants...")
  cesa = annotate_variants(cesa)

  current_snv_stats = maf[variant_type == "snv", .(num_samples = length(unique(Unique_Patient_Identifier)), num_snv = .N)]
  message(paste0("Loaded ", current_snv_stats$num_snv, " SNVs from ", current_snv_stats$num_samples, " samples into CESAnalysis."))
  
  return(cesa)
}

