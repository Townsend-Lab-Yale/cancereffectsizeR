#' Load MAF-formatted somatic mutation data
#' 
#' Load MAF data from a text file or data.frame/data.table into a CESAnalysis object. If column names
#' don't match MAF format specifications (Chromosome, Start_Position, etc., with Tumor_Sample_Barcode
#' used as the sample ID column), you can supply your own column names. When your CESAnalysis includes
#' chronological tumor progression states (e.g., pre- vs. post-treatment, stages 1-4, primary vs. metastatic),
#' specify "progression_col". By default, data is assumed to be derived from whole-exome sequencing. Whole-genome
#' data and targeted sequencing data are also supported when the "coverage" option is specified. If the data
#' you are loading is from a different genome build than your CESAnalysis, you can use the "chain_file" option 
#' to supply a UCSC-style chain file, and then your MAF coordinates will be automatically converted with liftOver.
#' 
#' 
#' @importFrom IRanges "%within%"
#' @param cesa the CESAnalysis object to load the data into
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.table or data.frame format
#' @param annotate Annotate mutations with gene and other reference information (required for effect size analysis)
#' @param sample_col column name with sample ID data (Tumor_Sample_Barcode or Unique_Patient_Identifier)
#' @param chr_col column name with chromosome data  (Chromosome)           
#' @param start_col column name with start position (Start_Position)
#' @param ref_col column name with reference allele data (Reference_Allele)
#' @param tumor_allele_col column name with alternate allele data; by default,
#'   values from Tumor_Seq_Allele2 and Tumor_Seq_Allele1 columns are used
#' @param progression_col column in MAF with tumor progression state (see docs)
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
#' @return CESAnalysis object with the specified MAF data loaded (in addition to any
#'   previously loaded data)
#' @export
load_maf = function(cesa = NULL, maf = NULL, annotate = TRUE, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", coverage = "exome", covered_regions = NULL,
                    covered_regions_name = NULL, covered_regions_padding = 0, progression_col = NULL, chain_file = NULL, enforce_generic_exome_coverage = FALSE) {
  
  if (is.null(cesa)) {
    stop("You need to supply a CESAnalysis object to load the MAF data into.")
  }
  
  if (! cesa@ref_key %in% ls(.ces_ref_data)) {
    preload_ref_data(cesa@ref_key)
  }
  
  # .ces_ref_data should be populated with reference data
  bsg = .ces_ref_data[[cesa@ref_key]]$genome
  
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
  cesa_anno_status = cesa@advanced$annotated
  if (is.null(cesa_anno_status) || identical(logical(0), cesa_anno_status)) {
    cesa@advanced$annotated = annotate
  } else if (cesa_anno_status == T & annotate == F){
    stop("The CESAnalysis already contains annotated variants, so you need to run with annotate = T.", call. = F)
  } else if (cesa_anno_status == F & annotate == T) {
    stop("The CESAnalysis already contains unannotated records. Either run annotate_variants() to annotate\n",
         "them, or re-run load_maf with annotate = F.", call. = F)
  }

  
  
  if (is.null(maf)) {
    stop("Supply MAF data via maf=[file path or data.table/data.frame].")
  }
  
  if (is.null(progression_col) & length(cesa@progressions) != 1) {
    stop("You must specify progression_col in the MAF since this CESAnalysis incorporates chronological tumor progression states.")
  }
  
  if (! is.null(progression_col) & length(cesa@progressions) == 1) {
    stop(paste0("This CESAnalysis does not incorporate tumor progression states, so you can't use the \"progression_col\" argument.\n",
                "Create a new CESAnalysis with \"progression_order\" specified to include this information."))
  }
  
  # validate coverage
  previous_covered_regions_names = names(cesa@coverage) # may be NULL
  if (! is.character(coverage) || ! coverage %in% c("exome", "genome", "targeted") || length(coverage) > 1) {
    stop("Argument coverage must be \"exome\", \"genome\", or \"targeted\"")
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
    # require letter/number in the name to avoid input mistakes
    if (! grepl("[0-9a-zA-Z]", covered_regions_name)) {
      stop("check your covered_regions_name (looks like just punctuation?)")
    }
    # don't allow user to name covered regions genome, exome, or exome+, since they have reserved use
    if(tolower(covered_regions_name) %in% c("exome", "exome+", "genome")) {
      stop("Sorry, your covered_regions_name is reserved for internal use. Please pick another name.")
    }
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
  
  # exome data can take covered_regions from a "generic exome" if one is present in the genome data directory
  using_generic_exome = FALSE
  if (coverage == "exome" && is.null(covered_regions)) {
    if(! check_for_genome_data(cesa, "generic_exome_gr")) {
      stop("This genome has no generic exome intervals, so to load exome data you must supply covered_regions (see docs)")
    }
    
    if("exome+" %in% names(cesa@coverage)) {
      covered_regions_name = "exome+"
      covered_regions = cesa@coverage[["exome+"]]
    } else {
      covered_regions_name = "exome"
      covered_regions = get_genome_data(cesa, "generic_exome_gr")
    }
    
    using_generic_exome = TRUE
    message(crayon::silver("Assuming this data has generic exome coverage (use the covered_regions argument if this isn't accurate)...."))
  }
  
  # custom covered_regions, when used, must be a GRanges object or a path to a BED-formatted text file
  genome_info = GenomeInfoDb::seqinfo(bsg)
  if(coverage != "genome" && ! using_generic_exome) {
    if (is(covered_regions, "character")) {
      if (length(covered_regions) != 1) {
        stop("covered_regions expected to be either a GRanges object or the path to a BED file.")
      }
      if (! file.exists(covered_regions)) {
        "covered_regions BED file could not be found (check the path)"
      }
      covered_regions = rtracklayer::import(covered_regions, format = "bed")
    } else if (! is(covered_regions, "GRanges")) {
      stop("Argument covered_regions expected to be a GRanges object or the path to a BED-formatted text file.")
    }
    
    
    # set coverage gr to match CESAnalysis genome (if this fails, possibly the genome build does not match)
    GenomeInfoDb::seqlevelsStyle(covered_regions) = "NCBI"
    
    tryCatch({
      msg = paste0("Your covered_regions don't seem compatible with the CESAnalysis reference genome.\n",
                    "Make sure it uses the same genome assembly. It may also help to subset to just the\n",
                    "primary chromosomes, if any obscure contigs are present in your regions.\n",
                   "Original warning/error:")
      GenomeInfoDb::seqlevels(covered_regions) = GenomeInfoDb::seqlevels(genome_info)
      GenomeInfoDb::seqinfo(covered_regions) = genome_info
    }, error = function(e) {
      message(msg)
      stop(conditionMessage(e))
    }, warning = function(w) {
      message(msg)
      stop(conditionMessage(w))
    })

    # drop any metadata
    GenomicRanges::mcols(covered_regions) = NULL
    
    # sort, reduce, unstrand
    covered_regions = GenomicRanges::reduce(GenomicRanges::sort(covered_regions), drop.empty.ranges = T)
    GenomicRanges::strand(covered_regions) = "*"
    
    # require genome name to match the CESAnalysis (too many potential errors if we allow anonymous or mismatched genome)
    expected_genome = GenomeInfoDb::genome(bsg)[1]
    gr_genome = GenomeInfoDb::genome(covered_regions)[1]
    if (expected_genome != gr_genome) {
      stop(paste0("The genome name of the covered_regions GRanges (", gr_genome, ") does not match the CESAnalysis (",
                  expected_genome, ")."))
    }
    
    # Validate and apply covered_regions_padding
    if (! is(covered_regions_padding, "numeric") || length(covered_regions_padding) > 1 || covered_regions_padding < 0 ||
        covered_regions_padding - as.integer(covered_regions_padding) != 0) {
      stop("covered_regions_padding should be 1-length integer", call. = F)
    }
    if (covered_regions_padding > 1000) {
      warning(sprintf("covered_regions_padding=%d is awfully high!", covered_regions_padding), call. = F)
    }
    if (covered_regions_padding > 0) {
      # Suppress the out-of-range warning since we'll trim afterwards
      withCallingHandlers(
      {
        GenomicRanges::start(covered_regions) = GenomicRanges::start(covered_regions) - covered_regions_padding
        GenomicRanges::end(covered_regions) = GenomicRanges::end(covered_regions) + covered_regions_padding
      }, warning = function(w) 
        {
          if (grepl("out-of-bound range", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      )
      covered_regions = GenomicRanges::reduce(GenomicRanges::trim(covered_regions))
    }
    
    
    # if covered_regions_name was already used in a previous load_maf call, make sure granges match exactly
    # otherwise, add the coverage information to the CESAnalysis
    if (covered_regions_name %in% names(cesa@coverage)) {
      if (! identical(covered_regions, cesa@coverage[[covered_regions_name]])) {
        stop("MAF data was previously loaded in using the same covered_regions_name (", covered_regions_name, "),\n",
             "but the covered_regions do not exactly match. Perhaps the input BED files (or GRanges) are from\n",
             "different sources, or different values of covered_regions_padding were used.")
      }
    }
    # make sure it really looks like exome data, if possible
    if (coverage == "exome" & check_for_genome_data(cesa, "generic_exome_gr")) {
      covered_regions_bases_covered = sum(IRanges::width(IRanges::ranges(covered_regions)))
      generic_bases_covered = sum(IRanges::width(IRanges::ranges(get_genome_data(cesa, "generic_exome_gr"))))
      if (covered_regions_bases_covered / generic_bases_covered < .4) {
        warning(paste0("Coverage is set to exome, but your covered_regions are less than 40% of the size of this genome's default exome intervals.\n",
                       "This might make sense if your exome capture array is very lean, but if this is actually targeted sequencing data,\n",
                       "create a new CESAnalysis and run load_maf() with coverage = \"targeted\"."))
      }
    }
  } else {
    # don't allow covered_regions_padding when it's not being used
    if (length(covered_regions_padding) != 1 || covered_regions_padding != 0) {
      stop("You can't use covered_regions_padding without supplying covered_regions.", call. = F)
    }
  }
  
  # genome data always has full coverage; other data has its coverage GRange saved in the @coverage list
  if (coverage != "genome") {
    cesa@coverage[[covered_regions_name]] = covered_regions
  }
  
  
  select_cols = c(sample_col, chr_col, start_col, ref_col, progression_col)
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
    message(crayon::silver("Found column Unique_Patient_Identifier; we'll assume this is the correct sample ID column."))
  }
  cols_to_check = c(sample_col, ref_col, chr_col, start_col, progression_col)
  if (tumor_allele_col != "guess") {
    cols_to_check = c(cols_to_check, tumor_allele_col)
  }
  
  for (col in cols_to_check) {
    if (! col %in% input_maf_cols) {
      missing_cols = c(missing_cols, col)
    }
  }
  
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
    warning(paste0(diff, " records in the input MAF had NA values in required columns. These have been dropped from analysis."))
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
      message(crayon::silver("Found column Tumor_Allele; we'll assume this is the correct tumor allele column."))
    } else {
      # automated tumor allele determination requires Tumor_Seq_Allele1/Tumor_Seq_Allele2 columns
      # if these columns are present, the tumor_allele_adder function will handle capitalization and other validation
      allele1_col = "Tumor_Seq_Allele1"
      allele2_col = "Tumor_Seq_Allele2"
      
      if (allele1_col %in% colnames(maf) && ! allele2_col %in% colnames(maf)) {
        message(crayon::silver("Found column Tumor_Seq_Allele1 but not Tumor_Seq_Allele2;\nwe'll assume Tumor_Seq_Allele1 is the correct tumor allele column."))
        maf[[tumor_allele_col]] = toupper(maf[[allele1_col]])
      } else if (allele2_col %in% colnames(maf) && ! allele1_col %in% colnames(maf)) {
        message(crayon::silver("Found column Tumor_Seq_Allele2 but not Tumor_Seq_Allele1;\nwe'll assume Tumor_Seq_Allele2 is the correct tumor allele column."))
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
    # To-do: put BSGenome reference check before this check, and then put failing records from this check into the excluded table
  }
  
  # collect tumor progression information
  if (! is.null(progression_col)) {
    if (is.factor(maf[[progression_col]])) {
      message("Warning: You supplied tumor progression as a factor, but it will be converted to character.")
      message("The progression ordering will be what you supplied to CESAnalysis() with \"progression_order\",")
      message("regardless of your factor ordering.")
    }
    sample_progressions = as.character(maf[[progression_col]])
    if(any(is.na(sample_progressions))) {
      stop("Error: There are NA values in your sample progressions column.")
    }
  } else {
    sample_progressions = cesa@progressions[1] # indicates a stageless analysis
  }
  
  # select only the necessary columns and give column names that will stay consistent
  maf = maf[,c(..sample_col, ..chr_col, ..start_col, ..ref_col, ..tumor_allele_col)]
  colnames(maf) = c("Unique_Patient_Identifier", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
  
  new_samples = data.table(Unique_Patient_Identifier = maf$Unique_Patient_Identifier, progression_name = sample_progressions)
  new_samples = new_samples[, .(progression_name = unique(progression_name)), by = "Unique_Patient_Identifier"]

  # associate each progression name with its chronological index (first progression stage is 1, next is 2, etc.)
  new_samples[, progression_index := sapply(progression_name, function(x) which(cesa@progressions == x)[1])]
  new_samples[, coverage := coverage]
  new_samples[, covered_regions := covered_regions_name]
  
  # ensure no sample has an illegal progression
  bad_progressions = setdiff(new_samples[, unique(progression_name)], cesa@progressions)
  if(length(bad_progressions) > 0) {
    stop(paste0("The following progressions were not declared in your CESAnalysis, but they were found in your MAF progressions column:\n",
                paste(bad_progressions, collapse = ", ")))
  }
  # see if any sample appears more than once in sample table (happens when one sample has multiple listed progressions)
  repeated_samples = new_samples[duplicated(Unique_Patient_Identifier), unique(Unique_Patient_Identifier)]
  if(length(repeated_samples) > 0) {
    stop(paste0("The following samples are associated with multiple progressions in the input data:\n", paste(repeated_samples, collapse=", ")))
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
  
  # warn the user if some of the declared progressions don't appear in the data at all
  if(length(cesa@progressions) > 1) {
    missing_progressions = cesa@progressions[! cesa@progressions %in% new_samples[,unique(progression_name)]]
    if (length(missing_progressions) > 0) {
      warning(paste0("The following tumor progression states were declared in your CESAnalysis, but they weren't present in the MAF data:\n",
                     paste(missing_progressions, collapse = ", ")))
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
        message(crayon::silver("Assuming that chrM in the chain file refers to mitochondrial DNA (aka chrMT)...."))
        message(crayon::silver("If you get reference mismatches on this contig, you may need a different chain file."))
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
      message(crayon::silver(paste0("Note: ", num_failing, " records failed liftOver, so they will be set aside.")))
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
  
  
  # remove any MAF records that are not in the coverage (except for generic exome data; then, just warn)
  if (coverage == "genome") {
    num_uncovered = 0
  } else {
    maf_grange = GenomicRanges::makeGRangesFromDataFrame(maf, seqnames.field = "Chromosome", start.field = "Start_Position", 
                                                         end.field = "Start_Position", seqinfo = genome_info)
    is_uncovered = ! maf_grange %within% covered_regions
    num_uncovered = sum(is_uncovered)
  }


  if (num_uncovered > 0) {
    total = nrow(maf)
    percent = round((num_uncovered / total) * 100, 1)
    if (using_generic_exome && enforce_generic_exome_coverage == FALSE) {
      # covered_regions is generic exome; take union with the sample data and call that exome+
      # when there is already an "exome+" covered_regions, included those in the intersection
      cesa@coverage[["exome+"]] = GenomicRanges::reduce(GenomicRanges::union(covered_regions, maf_grange[is_uncovered]))
      covered_regions_name = "exome+"
      new_samples[, covered_regions := covered_regions_name]
      
      # mark previous generic exome samples as exome+
      if(cesa@samples[, .N] > 0) {
        cesa@samples[covered_regions_name == "exome", covered_regions := covered_regions_name]
      }
      
      # warn if a lot of records are in uncovered areas; may indicate whole-genome data or low quality exome data
      if (percent > 10) {
        warning(paste0("More than 10% of MAF records are not within the CES genome's generic exome intervals.\n",
                       "Could this be whole-genome data? Or if you know the true covered regions, supply them\n",
                       "with the covered_regions argument."))
      }
      message(paste0("Note: ", num_uncovered, " MAF records out of ", total, " (", percent, 
                     "%) are at loci not covered by the generic exome capture regions\n",
                     "in this genome's CES reference data. It will be assumed that all generic exome samples in this\n",
                     "CESAnalysis, including those from previous or future load_maf calls, have coverage at\n",
                     "all of these sites; we'll call this \"exome+\" coverage.\n",
                     "To instead exclude uncovered records, re-run with enforce_generic_exome_coverage = TRUE."))
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
  
  
  
  # check for potential DNP/TNPs and separate them from data set
  # this is done before any other filtering to ensure catching as many of them as possible
  message(crayon::silver("Searching for possible multinucleotide variants..."))
  num.prefilter = nrow(maf)
  dnp_tnp_results = DNP_TNP_remover(maf)
  maf = dnp_tnp_results$kept[,1:5] # To-do: fix return value of DNP_TNP
  pred.mnv.maf = dnp_tnp_results$removed
  num.pred.mnv = nrow(pred.mnv.maf)
  
  if (num.pred.mnv > 0) {
    pred.mnv.maf$Exclusion_Reason = "predicted_MNV"
    excluded = rbind(excluded, pred.mnv.maf)
    # To-do: move message to DNP_TNP_remover or otherwise ensure this description remains accurate
    percent = round((num.pred.mnv / num.prefilter) * 100, 1)
    message(paste0("Note: ", num.pred.mnv, " mutation records out of ", num.prefilter, " (", percent, "%) ",
                   "are within 2 bp of other mutations in the same tumors. "))
    message("These records will be excluded since they may be multi-nucleotide mutation events rather than independent SNVs.")
  }
  
  # No support for mitochondrial mutations yet
  is_mt <- maf[,Chromosome] == "MT"
  if (any(is_mt)) {
    mt_maf = maf[is_mt,]
    mt_maf$Exclusion_Reason = "mitochondrial"
    excluded = rbind(excluded, mt_maf)
    maf <- maf[! is_mt,]
    message(paste("Note:", sum(is_mt), "mitochondrial mutations have been removed from the data (mitochondrial analysis not yet supported)."))
  }
  
  
  # Ensure reference alleles of mutations match reference genome (Note: Insertions won't match if their reference allele is "-")
  message(crayon::silver("Checking that reference alleles match reference genome..."))
  ref_allele_lengths = nchar(maf[, Reference_Allele])
  ref_alleles_to_test = maf[, Reference_Allele]
  end_pos = maf[, Start_Position] + ref_allele_lengths - 1 # for multi-base deletions, check that all deleted bases match reference
  
  reference_alleles <- as.character(BSgenome::getSeq(bsg, maf[,Chromosome],
                                                     strand="+", start=maf[,Start_Position], end=end_pos))
  num.prefilter = nrow(maf)
  
  # For insertions, can't evaluate if record matches reference since MAF reference will be just "-"
  # Could conceivably verify that the inserted bases don't match reference....
  reference.mismatch.maf = maf[Reference_Allele != reference_alleles & Reference_Allele != '-',]
  num.nonmatching = nrow(reference.mismatch.maf)
  prefilter_total_snv = maf[Reference_Allele %in% c("A", "C", "T", "G"), .N]
  maf = maf[Reference_Allele == reference_alleles | Reference_Allele == '-',] # filter out ref mismatches
  
  # if significant numbers of SNVs don't match reference, don't run
  mismatch_snv = reference.mismatch.maf[Reference_Allele %in% c("A", "C", "T", "G"), .N]
  
  bad_snv_frac = mismatch_snv / prefilter_total_snv
  bad_snv_percent = round(bad_snv_frac * 100, 1)
  
  if (num.nonmatching > 0) {
    reference.mismatch.maf$Exclusion_Reason = "reference_mismatch"
    excluded = rbind(excluded, reference.mismatch.maf)
    percent = round((num.nonmatching / num.prefilter) * 100, 1)
    
    message(crayon::silver(paste0("Note: ", num.nonmatching, " mutation records out of ", num.prefilter, " (", percent, "%, including ", bad_snv_percent,
                          "% of SNV records) have reference alleles that do not actually match the reference genome.")))
    message(crayon::silver("These records will be excluded from effect size analysis."))
    
    if(bad_snv_frac > .01) {
      warning(paste0(bad_snv_percent, "% of SNV records do not match the given reference genome. You should probably figure out why",
                     " and make sure that the rest of your data set is okay to use before continuing."))
    }
  }
   else {
    message(crayon::silver("Reference alleles look good."))
  }
  
  nt = c("A", "T", "C", "G")
  maf[Reference_Allele %in% nt & Tumor_Allele %in% nt, Variant_Type := "SNV"]
  maf[is.na(Variant_Type), Variant_Type := "indel"]
  
  # drop any samples that had all mutations excluded
  new_samples = new_samples[Unique_Patient_Identifier %in% maf$Unique_Patient_Identifier]
  cesa@samples = rbind(cesa@samples, new_samples)
  setcolorder(cesa@samples, c("Unique_Patient_Identifier", "coverage", "covered_regions", "progression_name", "progression_index"))
  setkey(cesa@samples, "Unique_Patient_Identifier")
  
  if (nrow(excluded) > 0) {
    colnames(excluded) = c(colnames(maf)[1:5], "Exclusion_Reason")
    cesa@excluded = rbind(cesa@excluded, excluded) 
  }
  
  
  if(annotate) {
    if (length(cesa@mutations) == 0) {
      message("Annotating variants...")
    } else {
      message("Annotating new variants...")
    }
    
    maf[Variant_Type == "SNV", snv_id := paste0(Chromosome, ':', Start_Position, '_', Reference_Allele, '>', Tumor_Allele)]
    
    # A key point here is that if an AAC is already annotated, all its constituent SNVs are, too; also,
    # when an SNV has been annotated using CES, any associated AACs are always put in the AAC table. 
    # Therefore, only SNVs that are in the new data that aren't already in the SNV table need to be annotated,
    # and there associated AACs are new, too
    new_variants = setdiff(na.omit(maf$snv_id), cesa@mutations$snv$snv_id)
    maf_to_annotate = maf[snv_id %in% new_variants]
    cesa_for_anno = suppressMessages(CESAnalysis(genome = cesa@ref_key))
    cesa_for_anno@maf = maf_to_annotate
    cesa_for_anno@samples = cesa@samples
    cesa_for_anno@coverage = cesa@coverage
    cesa_for_anno = annotate_variants(cesa_for_anno)
    
    # Add the pair of annotation columns from annotate_variants
    maf[cesa_for_anno@maf,  c("genes", "assoc_aac") := .(genes, assoc_aac), on = "snv_id"]
    maf[Variant_Type != "SNV", c("genes", "assoc_aac") := NA] # set to NA instead of 1-item list containing NULL

    # Rarely, variants get thrown out during annotation if their trinuc contexts are non-specific (N's)
    if (cesa_for_anno@excluded[, .N] > 0) {
      snvs_excluded_in_anno = cesa_for_anno@excluded[, paste0(Chromosome, ':', Start_Position, '_', Reference_Allele, '>', Tumor_Allele)]
      maf = maf[! snv_id %in% snvs_excluded_in_anno]
      cesa@excluded = rbind(cesa@excluded, cesa_for_anno@excluded)
    }
    
    # Update coverage of existing mutations, if there any and if the new MAF data uses new covered regions
    if(cesa@maf[, .N] > 0 & ! covered_regions_name %in% previous_covered_regions_names & covered_regions_name != "genome") {
      prev_snv = cesa@mutations$snv
      snv_gr = GenomicRanges::makeGRangesFromDataFrame(prev_snv, seqnames.field = "chr", start.field = "pos", end.field = "pos")
      is_covered = snv_gr %within% cesa@coverage[[covered_regions_name]]
      prev_snv[is_covered, covered_in := list(lapply(covered_in, function(x) c(x, covered_regions_name))), by = "snv_id"]
      cesa@mutations$snv = prev_snv
      
      # Apply to amino acid changes
      aac_coverage = prev_snv[, .(aac_id = unlist(assoc_aac), covered_in), by = "snv_id"]
      aac_coverage = aac_coverage[, .(covered_in = list(sort(unique(unlist(covered_in))))), by = "aac_id"]
      cesa@mutations$amino_acid_change[aac_coverage, covered_in := covered_in, on = "aac_id"]
    }
    cesa@mutations$amino_acid_change = rbind(cesa@mutations$amino_acid_change, cesa_for_anno@mutations$amino_acid_change)
    cesa@mutations$snv = rbind(cesa@mutations$snv, cesa_for_anno@mutations$snv)
    setkey(cesa@mutations$amino_acid_change, 'aac_id')
    setkey(cesa@mutations$snv, 'snv_id')
  }
  
  cesa@maf = rbind(cesa@maf, maf)
  
  
  
  current_snv_stats = maf[Variant_Type == "SNV", .(num_samples = length(unique(Unique_Patient_Identifier)), num_snv = .N)]
  message(paste0("Loaded ", current_snv_stats$num_snv, " SNVs from ", current_snv_stats$num_samples, " samples into CESAnalysis."))
  
  return(cesa)
}

