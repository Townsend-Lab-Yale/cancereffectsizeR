#' Load MAF-formatted somatic mutation data
#' 
#' Load MAF data from a text file or data.frame/data.table into a CESAnalysis object. If
#' column names don't match MAF format specifications (Chromosome, Start_Position, etc.,
#' with Tumor_Sample_Barcode used as the sample ID column), you can supply your own column
#' names. When your CESAnalysis has defined sample groups specify "group_col". By default,
#' data is assumed to be derived from whole-exome sequencing. Whole-genome data and
#' targeted sequencing data are also supported when the "coverage" option is specified. If
#' the data you are loading is from a different genome build than your CESAnalysis, you
#' can use the "chain_file" option to supply a UCSC-style chain file, and then your MAF
#' coordinates will be automatically converted with liftOver.
#' 
#' @param cesa the CESAnalysis object to load the data into
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.table or data.frame format
#' @param annotate Annotate mutations with gene and other reference information (required for effect size analysis)
#' @param sample_col column name with sample ID data (Tumor_Sample_Barcode or Unique_Patient_Identifier)
#' @param chr_col column name with chromosome data  (Chromosome)           
#' @param start_col column name with start position (Start_Position)
#' @param ref_col column name with reference allele data (Reference_Allele)
#' @param tumor_allele_col column name with alternate allele data; by default,
#'   values from Tumor_Seq_Allele2 and Tumor_Seq_Allele1 columns are used
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
#' @param chain_file a chain file (text format, name ends in .chain) to convert MAF
#'   records to the genome build used in the CESAnalysis
#' @param enforce_generic_exome_coverage when loading generic exome data, exclude records
#'   that aren't covered in the (somewhat arbitrary) generic exome intervals included with
#'   CES genome reference data (default FALSE)
#' @return CESAnalysis with the specified MAF data loaded
#' @export
load_maf = function(cesa = NULL, maf = NULL, annotate = TRUE, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", coverage = "exome", covered_regions = NULL,
                    covered_regions_name = NULL, covered_regions_padding = 0, group_col = NULL, chain_file = NULL, enforce_generic_exome_coverage = FALSE) {
  
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  
  cesa = update_cesa_history(cesa, match.call())
  
  # Need RefCDS if annotating
  if (! cesa@ref_key %in% ls(.ces_ref_data) & annotate == T) {
    preload_ref_data(cesa@ref_data_dir)
  }
  bsg = get_cesa_bsg(cesa)
  genome_info = GenomeInfoDb::seqinfo(bsg)
  
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
  
  if (! is(annotate, "logical") || length(annotate) != 1) {
    stop("annotate should be T/F", call. = F)
  }
  
  # We're not allowing a mix of annotated and unannotated records
  # The purpose of annotate = FALSE is to allow quick MAF loading for scratch work, etc.
  cesa_anno_status = cesa@advanced$annotated # starts false until first annotated data is loaded
  if(cesa_anno_status == F & annotate == T & cesa@samples[, .N] == 0) {
    cesa@advanced$annotated = TRUE
  } else if(cesa_anno_status == T & annotate == F){
    stop("The CESAnalysis already contains annotated variants, so you need to run with annotate = T.", call. = F)
  } else if (cesa_anno_status == F & annotate == T) {
    stop("The CESAnalysis already contains unannotated records. Either run annotate_variants() to annotate\n",
         "them, or re-run load_maf with annotate = F.", call. = F)
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
  
  # validate covered_regions (required for targeted, optional for exome, prohibited for genome)
  if (coverage == "genome") {
    if(! is.null(covered_regions)) {
      stop("covered_regions should be left NULL when coverage is \"genome\"")
    }
    covered_regions_name = "genome"
  }
  
  if (coverage == "targeted" & is.null(covered_regions) ) {
    stop("can't load targeted data without covered_regions (see docs)")
  }
  
  if (! is(enforce_generic_exome_coverage, "logical") || length(enforce_generic_exome_coverage) != 1) {
    stop("enforce_generic_exome_coverage should be T/F")
  }
  
  if (coverage == "exome" & is.null(covered_regions)) {
    covered_regions_name = "exome"
    if(! check_for_ref_data(cesa, "generic_exome_gr")) {
      stop("This genome has no generic exome intervals, so to load exome data you must supply covered_regions (see docs)")
    }
  }
  
  # Determine whether we're using "exome" or "exome+" for generic exome data
  # Usually, the lenient exome+ option is used, but the choice must be consistent throughout a CESAnalysis
  if(covered_regions_name == "exome") {
    # Whether or not CESAnalysis uses "exome+", once the first generic exome data has been loaded,
    # there will always be generic exome intervals under "exome" in the coverage list
    if (! "exome" %in% names(cesa@coverage[["exome"]])) {
      cesa@advanced$using_exome_plus = ! enforce_generic_exome_coverage
      covered_regions = get_ref_data(cesa, "generic_exome_gr")
      cesa@coverage$exome[["exome"]] = covered_regions
      
    } else if (enforce_generic_exome_coverage == TRUE & cesa@advanced$using_exome_plus) {
      stop("You can't have enforce_generic_exome_coverage = TRUE because you previously loaded generic exome data\n",
           "with enforce_generic_exome_coverage = FALSE")
    } else if (enforce_generic_exome_coverage == FALSE & ! cesa@advanced$using_exome_plus) {
      stop("You can't have enforce_generic_exome_coverage == FALSE because you previously loaded generic exome data\n",
           "with enforce_generic_exome_coverage = TRUE.")
    } else {
      covered_regions = cesa@coverage$exome[["exome"]]
    }
    
    # For exome+, use GRanges if available
    # Otherwise, for exome or exome+, load the generic exome intervals
    if (! enforce_generic_exome_coverage) {
      covered_regions_name = "exome+"
      if("exome+" %in% names(cesa@coverage$exome)) {
        covered_regions = cesa@coverage$exome[["exome+"]]
      } else {
        cesa@coverage[["exome"]][["exome+"]] = covered_regions # that is, the generic exome regions
      }
    } else {
      if("exome" %in% names(cesa@coverage$exome)) {
        covered_regions = cesa@coverage$exome[["exome"]]
      } else {
        covered_regions = get_ref_data(cesa, "generic_exome_gr")
        cesa@coverage$exome[["exome"]] = covered_regions
      }
    }
    pretty_message("Assuming this data has generic exome coverage (it's better to supply covered intervals if you have them; see docs)...")
  } else if (! is.null(covered_regions)) {
    # Use internal function to avoid updating covered_in now (annotate_variants step will handle this larger)
    cesa = .add_covered_regions(cesa = cesa, covered_regions = covered_regions, covered_regions_padding = covered_regions_padding, 
                               coverage_type = coverage, covered_regions_name = covered_regions_name, update_anno = FALSE)
  }
  
  
  select_cols = c(sample_col, chr_col, start_col, ref_col, group_col)
  if (tumor_allele_col == "guess") {
    select_cols = c(select_cols, "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Allele")
  } else {
    select_cols = c(select_cols, tumor_allele_col)
  }
  # when sample_col has default value, will also check for Unique_Patient_Identifier if the default isn't found (since CESAnalysis creates this column)
  if (sample_col == "Tumor_Sample_Barcode") {
    select_cols = c(select_cols, "Unique_Patient_Identifier")
  }
  
  bad_maf_msg = "Input MAF is expected to be a data frame or the filename of an MAF-formatted tab-delimited text file."
  if ("character" %in% class(maf)) {
    if(length(maf) > 1) {
      stop(bad_maf_msg)
    }
    if(file.exists(maf)) {
      # read in file and suppress warnings about missing columns since this function handles the validation
      withCallingHandlers(
      {
        
        maf = fread(maf, sep = "\t", quote ="", blank.lines.skip = T, select = select_cols)
      },
      error = function(e) {
        message("Unable to read specified MAF file:")
      },
      warning = function(w) {
          if (grepl("not found in column name header", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
      })
    } else {
      stop("Specified MAF not found.")
    }
  }
  if(! "data.frame" %in% class(maf)) {
    stop(bad_maf_msg)
  } else if(nrow(maf) == 0) {
    stop("Input MAF data set is empty.")
  } else {
    # user supplied pre-loaded data.frame or data.table
    # take just necessary columns and convert to data.table if necessary
    cols_to_keep = names(maf)[which(names(maf) %in% select_cols)]
    if ("data.table" %in% class(maf)) {
      maf = maf[, ..cols_to_keep]
    } else {
      maf = maf[, cols_to_keep]
      maf = data.table(maf)
    }
  }
  missing_cols = character()
  input_maf_cols = colnames(maf)
  
  if (sample_col == "Tumor_Sample_Barcode" & ! sample_col %in% input_maf_cols & "Unique_Patient_Identifier" %in% input_maf_cols) {
    sample_col = "Unique_Patient_Identifier"
    pretty_message("Found column Unique_Patient_Identifier; we'll assume this is the correct sample ID column.")
  }
  cols_to_check = c(sample_col, ref_col, chr_col, start_col, group_col)
  if (tumor_allele_col != "guess") {
    cols_to_check = c(cols_to_check, tumor_allele_col)
  }
  missing_cols = setdiff(cols_to_check, input_maf_cols)

  if (length(missing_cols) > 0) {
    missing_cols = paste(missing_cols, collapse = ", " )
    msg = "The following MAF columns couldn't be found:"
    msg = paste(msg, missing_cols, sep = "\n\t")
    stop(msg)
  }
  
  # convert all columns to character, then convert Start_Position to numeric
  maf = maf[, names(maf) := lapply(.SD, as.character)]
  maf[[start_col]] = as.numeric(maf[[start_col]])
  
  # drop columns with NAs or empty-ish strings
  nrow_orig = nrow(maf)
  maf = na.omit(maf)
  nrow_curr = nrow(maf)
  if(nrow_curr != nrow_orig) {
    diff = nrow_orig - nrow_curr
    warning(paste0(diff, " rows in the input MAF had NA values in required columns, so they were dropped."))
  }
  looks_empty = apply(maf, 1, function(x) any(grepl("^[.\"\' ]*$", x)))
  num_empty = sum(looks_empty)
  if(num_empty > 0) {
    maf = maf[! looks_empty,]
    warning(paste0(num_empty, " MAF records had fields that looked empty (just whitespace, quotation marks or periods).\n",
                   "These records have been removed from analysis. (One possible cause is indels not being in proper MAF format.)"))
  }
  
  # uppercase bases only
  maf[[ref_col]] = toupper(maf[[ref_col]])
  
  # figure out which column has correct tumor allele data
  if(tumor_allele_col == "guess") {
    tumor_allele_col = "Tumor_Allele"
    if (tumor_allele_col %in% colnames(maf)) {
      pretty_message("Found column Tumor_Allele; we'll assume this is the correct tumor allele column.")
    } else {
      # automated tumor allele determination requires Tumor_Seq_Allele1/Tumor_Seq_Allele2 columns
      allele1_col = "Tumor_Seq_Allele1"
      allele2_col = "Tumor_Seq_Allele2"
      
      if (allele1_col %in% colnames(maf) && ! allele2_col %in% colnames(maf)) {
        pretty_message("Found column Tumor_Seq_Allele1 but not Tumor_Seq_Allele2; we'll assume Tumor_Seq_Allele1 is the correct tumor allele column.")
        maf[[tumor_allele_col]] = toupper(maf[[allele1_col]])
      } else if (allele2_col %in% colnames(maf) && ! allele1_col %in% colnames(maf)) {
        pretty_message("Found column Tumor_Seq_Allele2 but not Tumor_Seq_Allele1;\nwe'll assume Tumor_Seq_Allele2 is the correct tumor allele column.")
        maf[[tumor_allele_col]] = toupper(maf[[allele2_col]])
      } else if (! allele1_col %in% colnames(maf) | ! allele2_col %in% colnames(maf)) {
        stop(paste0("Tumor alleles can't be determined automatically deduced without Tumor_Seq_Allele1 ",
                    "and/or Tumor_Seq_Allele2 columns. Please manually specify with \"tumor_allele_col=...\""))
      } else {
        maf$Tumor_Seq_Allele1 = toupper(maf$Tumor_Seq_Allele1)
        maf$Tumor_Seq_Allele2 = toupper(maf$Tumor_Seq_Allele2)
        
        # take allele 2 as the tumor allele, but when it matches ref, replace with allele 1
        # if that is still equal to ref, record will later be discarded
        tumor_alleles = maf$Tumor_Seq_Allele2
        allele_2_matches_ref = maf$Tumor_Seq_Allele2 == maf[[ref_col]]
        tumor_alleles[allele_2_matches_ref] = maf$Tumor_Seq_Allele1[allele_2_matches_ref]
        maf[[tumor_allele_col]] <- tumor_alleles
      }
    }
  }
  
  # uppercase bases only
  maf[[tumor_allele_col]] = toupper(maf[[tumor_allele_col]])
  
  # drop records where tumor allele is equal to reference allele
  no_variant = maf[[ref_col]] == maf[[tumor_allele_col]]
  num_unvaried = sum(no_variant)
  if(num_unvaried > 0) {
    maf = maf[!no_variant]
    warning(paste0(num_unvaried, " MAF records had tumor alleles identical to reference alleles; these were removed from analysis.\n"))
    # To-do: put BSgenome reference check before this check, and then put failing records from this check into the excluded table
  }
  
  # collect sample group information
  if (! is.null(group_col)) {
    if (is.factor(maf[[group_col]])) {
      warning("You supplied a sample group column as a factor, but it was converted to character.\n",
              "The group ordering will be what you supplied to CESAnalysis() with \"sample_groups\",\n",
              "regardless of factor ordering.")
    }
    sample_groups = as.character(maf[[group_col]])
    if(any(is.na(sample_groups))) {
      stop("Error: There are NA values in your sample groups column.")
    }
  } else {
    sample_groups = cesa@groups[1] # indicates a stageless analysis
  }
  
  # select only the necessary columns and give column names that will stay consistent
  maf = maf[,c(..sample_col, ..chr_col, ..start_col, ..ref_col, ..tumor_allele_col)]
  colnames(maf) = c("Unique_Patient_Identifier", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
  
  new_samples = data.table(Unique_Patient_Identifier = maf$Unique_Patient_Identifier, group = sample_groups)
  new_samples = new_samples[, .(group = unique(group)), by = "Unique_Patient_Identifier"]

  # associate each group name with its order (for analyses with unordered groupings, this may be arbitrary)
  new_samples[, group_index := sapply(group, function(x) which(cesa@groups == x)[1])]
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
  
  # warn the user if some of the declared groups don't appear in the data at all
  if(length(cesa@groups) > 1) {
    missing_groups = cesa@groups[! cesa@groups %in% new_samples[,unique(group)]]
    if (length(missing_groups) > 0) {
      warning(paste0("The following groups were declared in your CESAnalysis, but they weren't present in the MAF data:\n",
                     paste(missing_groups, collapse = ", ")), call. = F)
    }    
  }
  
  # create separate table for excluded records
  excluded = data.table()
  
  # strip chr prefixes from chr column, if present ("NCBI-style")
  maf[, Chromosome := sub('^chr', '', Chromosome)]
  # temporary handling of MT/M
  maf[Chromosome == "M", Chromosome:="MT"]  # changing M to MT
  
  # handle liftOver to the assembly of the CESAnalysis
  if(need_to_liftOver) {
    chain = rtracklayer::import.chain(chain_file)
    names(chain) = sub("^chr", "", names(chain))
    # avoid potential NCBI/UCSC incompatibilities between the user's data and the chain file
    if(any(maf$Chromosome == "MT")) {
      if(! "MT" %in% names(chain) && "M" %in% names(chain)) {
        names(chain) = sub("^M$", "MT", names(chain))
        pretty_message("Assuming that chrM in the chain file refers to mitochondrial DNA (aka chrMT)....")
        pretty_message("If you get reference mismatches on this contig, you may need a different chain file.")
      }
    }
    maf[, rn := 1:.N] # using row number as an identifier to know which intervals fail liftover
    for_liftover = maf[,.(Chromosome, Start_Position, rn)] 
    gr = GenomicRanges::makeGRangesFromDataFrame(df = for_liftover, seqnames.field = "Chromosome", 
                                                 start.field = "Start_Position", end.field = "Start_Position",
                                                  keep.extra.columns = T)
    lifted_over = unlist(rtracklayer::liftOver(gr, chain))
    GenomeInfoDb::seqlevelsStyle(lifted_over) = "NCBI"
    lifted_over = as.data.table(lifted_over)
    
    merged_maf = data.table::merge.data.table(maf, lifted_over, by = "rn")
    merged_maf[, Chromosome := as.character(seqnames)] # seqnames comes back as factor!
    merged_maf[, Start_Position := start]

    failed_liftover = maf[! rn %in% merged_maf$rn,]
    num_failing = nrow(failed_liftover)
    if (num_failing > 0) {
      failed_liftover$Exclusion_Reason = "failed liftOver"
      excluded = rbind(excluded, failed_liftover[, -"rn"])
      pretty_message(paste0("Note: ", num_failing, " records failed liftOver, so they will be set aside."))
    }
    
    maf = merged_maf[, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
    
    # different loci in one genome can get lifted to the same position in the next, due to fixes
    # rarely, mutations at multiple matching sites can get called in a sample, resulting in duplicate records after liftOver
    lifted_to_same = duplicated(maf[,.(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele)])
    maf = maf[! lifted_to_same]
  }
  
  # discard any records with chromosomes not present in reference
  illegal_chroms = maf[! Chromosome %in% GenomeInfoDb::seqlevels(genome_info), unique(Chromosome)]
  if (length(illegal_chroms) > 0) {
    has_bad_chr = maf$Chromosome %in% illegal_chroms
    bad_chr_maf = maf[has_bad_chr]
    bad_chr_maf$Exclusion_Reason = "illegal_chromosome_name"
    maf = maf[! has_bad_chr]
    excluded = rbind(excluded, bad_chr_maf)
    message(paste0("Note: ", length(illegal_chroms), " records excluding for having chromosome names that don't match the genome. "))
  }
  
  # check for multi-nucleotide variants and separate them from MAF data
  # MNVs are only possible in sample/chromosome combinations with more than one MAF record
  poss_mnv = maf[order(Unique_Patient_Identifier, 
                       Chromosome, Start_Position)][, .SD[.N > 1] , by = c("Unique_Patient_Identifier", "Chromosome")]
  
  if (poss_mnv[, .N] > 0) {
    poss_mnv[, dist_to_prev := c(Inf, diff(Start_Position)), by = c("Unique_Patient_Identifier", "Chromosome")]
    poss_mnv[, dist_to_next := c(dist_to_prev[2:.N], Inf), by = c("Unique_Patient_Identifier", "Chromosome")]
    poss_mnv[dist_to_prev < 3 | dist_to_next < 3, is_mnv := T]
    maf[poss_mnv, is_mnv := is_mnv, on = c("Unique_Patient_Identifier", "Chromosome", "Start_Position")]
    mnv_rows = maf[is_mnv == T, which = T]
    num_prefilter = nrow(maf)
    num_mnv = length(mnv_rows)
    
    if (num_mnv > 0) {
      mnv_records = maf[mnv_rows, 1:5]
      mnv_records$Exclusion_Reason = "predicted_MNV"
      maf = maf[! mnv_rows][, is_mnv := NULL]
      excluded = rbind(excluded, mnv_records)
      # To-do: move message to DNP_TNP_remover or otherwise ensure this description remains accurate
      percent = round((num_mnv / num_prefilter) * 100, 1)
      msg = paste0("Note: ", num_mnv, " MAF records (", percent, "%) ",
                   "are within 2 bp of other mutations in the same tumors. These records will not be counted as SNVs ",
                   "since they likely did not arise from independent events (i.e., they're multi-nucleotide variants).")
      pretty_message(msg)
    }
  }
  
  
  
  # remove any MAF records that are not in the coverage, unless generic exome with enforce_generic_exome_coverage = FALSE
  if (covered_regions_name == "genome") {
    num_uncovered = 0
  } else {
    maf_grange = GenomicRanges::makeGRangesFromDataFrame(maf, seqnames.field = "Chromosome", start.field = "Start_Position", 
                                                         end.field = "Start_Position", seqinfo = genome_info)
    
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
        warning(paste0("More than 10% of MAF records are not within the CES genome's generic exome intervals.\n",
                       "Could this be whole-genome data? Or if you know the true covered regions, supply them\n",
                       "with the covered_regions argument."))
      }
      msg = paste0("Note: ", num_uncovered, " MAF records (", percent, 
                   "%) are at loci not covered in the default exome intervals ",
                   "in this CESAnalysis's reference data. The samples in this MAF (along with along other MAFs ",
                   "loaded as generic exome data), will be assumed to all share the same coverage, ",
                   "which we'll call \"exome+\". To instead exclude uncovered records across all generic WES data, ",
                   "create a new CESAnalysis and load data with enforce_generic_exome_coverage = TRUE.")
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
  
  
  # No support for mitochondrial mutations yet
  is_mt <- maf[,Chromosome] == "MT"
  if (any(is_mt)) {
    mt_maf = maf[is_mt,]
    mt_maf$Exclusion_Reason = "mitochondrial"
    excluded = rbind(excluded, mt_maf)
    maf <- maf[! is_mt,]
    msg = paste("Note:", sum(is_mt), "mitochondrial variants have been excluded (sadly, mitochondrial analysis is not yet supported).")
    pretty_message(msg)
  }
  
  
  # Ensure reference alleles of mutations match reference genome (Note: Insertions won't match if their reference allele is "-")
  message("Checking that reference alleles match the reference genome...")
  ref_allele_lengths = nchar(maf[, Reference_Allele])
  ref_alleles_to_test = maf[, Reference_Allele]
  end_pos = maf[, Start_Position] + ref_allele_lengths - 1 # for multi-base deletions, check that all deleted bases match reference
  
  reference_alleles <- as.character(BSgenome::getSeq(bsg, maf[,Chromosome],
                                                     strand="+", start=maf[,Start_Position], end=end_pos))
  num_prefilter = nrow(maf)
  
  # For insertions, can't evaluate if record matches reference since MAF reference will be just "-"
  # Could conceivably verify that the inserted bases don't match reference....
  ref_mismatch_records = maf[Reference_Allele != reference_alleles & Reference_Allele != '-',]
  num_ref_mismatch = nrow(ref_mismatch_records)
  prefilter_total_snv = maf[Reference_Allele %in% c("A", "C", "T", "G"), .N]
  maf = maf[Reference_Allele == reference_alleles | Reference_Allele == '-',] # filter out ref mismatches
  
  # if significant numbers of SNVs don't match reference, don't run
  mismatch_snv = ref_mismatch_records[Reference_Allele %in% c("A", "C", "T", "G"), .N]
  
  bad_snv_frac = mismatch_snv / prefilter_total_snv
  bad_snv_percent = round(bad_snv_frac * 100, 1)
  
  if (num_ref_mismatch > 0) {
    ref_mismatch_records$Exclusion_Reason = "reference_mismatch"
    excluded = rbind(excluded, ref_mismatch_records)
    percent = round((num_ref_mismatch / num_prefilter) * 100, 1)
    
    msg = paste0("Note: ", num_ref_mismatch, " MAF records out of ", num_prefilter, " (", percent, "%, including ", bad_snv_percent,
                          "% of SNV records) are excluded for having reference alleles that do not match the reference genome.")
    pretty_message(msg)
    if(mismatch_snv > 0) {
      warning(mismatch_snv, " SNV records do not match the given reference genome. You should probably figure out why",
                     " and make sure that the rest of your data set is okay to use before continuing.")
    }
  }
   else {
    pretty_message("Reference alleles look good.")
  }
  
  nt = c("A", "T", "C", "G")
  maf[Reference_Allele %in% nt & Tumor_Allele %in% nt, variant_type := "snv"]
  maf[is.na(variant_type), variant_type := "indel"]
  
  # drop any samples that had all mutations excluded
  new_samples = new_samples[Unique_Patient_Identifier %in% maf$Unique_Patient_Identifier]
  cesa@samples = rbind(cesa@samples, new_samples)
  setcolorder(cesa@samples, c("Unique_Patient_Identifier", "coverage", "covered_regions", "group", "group_index"))
  setkey(cesa@samples, "Unique_Patient_Identifier")
  
  if (nrow(excluded) > 0) {
    colnames(excluded) = c(colnames(maf)[1:5], "Exclusion_Reason")
    cesa@excluded = rbind(cesa@excluded, excluded) 
  }
  
  cesa@maf = rbind(cesa@maf, maf, fill = T)
  if(annotate) {
    message("Annotating variants...")
    cesa@advanced$recording = F
    cesa = annotate_variants(cesa)
    cesa@advanced$recording = T
  }

  current_snv_stats = maf[variant_type == "snv", .(num_samples = length(unique(Unique_Patient_Identifier)), num_snv = .N)]
  message(paste0("Loaded ", current_snv_stats$num_snv, " SNVs from ", current_snv_stats$num_samples, " samples into CESAnalysis."))
  
  return(cesa)
}

