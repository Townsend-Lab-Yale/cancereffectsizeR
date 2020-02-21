#' load_maf
#' @description Load MAF-formatted mutation data into a CESAnalysis object (the core cancereffectsizeR object)
#' @importFrom IRanges "%within%"
#' @param cesa the CESAnalysis object to load the data into
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.table or data.frame format
#' @param sample_col column name with sample ID data
#' @param chr_col column name with chromosome data             
#' @param start_col column name with start position ( currently only analyze SNVs, so no end position needed)
#' @param ref_col column name with reference allele data
#' @param tumor_allele_col column name with alternate allele data (defaults to "guess", which examines ref_col,
#'                         "Tumor_Seq_Allele1", and "Tumor_Seq_Allele2" and determines from there)
#' @param covered_regions for panel sequencing data, a BED file of covered interals, or a data frame with first three columns chr, start (1-based), end (inclusive)
#' @param progression_col column in MAF with subset data (e.g., column contains data like "Primary" and "Metastatic" in each row)
#' @param chain_file a chain file (text format, name ends in .chain) to convert MAF records to the genome build used in the CESAnalysis
#' @return CESAnalysis object with somewhat cleaned-up MAF data added 
#' @export
load_maf = function(cesa = NULL, maf = NULL, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", covered_regions = "exome",
                    progression_col = NULL, chain_file = NULL) {

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
  if (is.null(cesa)) {
    stop("You need to supply a CESAnalysis object to load the MAF data into.")
  }
  
  if (is.null(maf)) {
    stop("Supply MAF data via maf=[file path or data.table/data.frame].")
  }
  if (is.null(progression_col) & length(cesa@progressions@order) != 1) {
    stop("You must supply a progression_col for your MAF since your analysis covers multiple progression stages")
  }
  
  if (! is.null(progression_col) & length(cesa@progressions@order) == 1) {
    stop(paste0("This CESAnalysis is not stage-specific, so you can't use the \"progression_col\" argument.\n",
                "Create a new CESAnalysis with \"progression_order\" specified to include this information.")) 
  }
  
  if (is.null(covered_regions)) {
    stop("You must supply covered_regions. For whole-exome/whole-genome data, use the default setting \"exome\".")
  }
  
  select_cols = c(sample_col, chr_col, start_col, ref_col, progression_col)
  if (tumor_allele_col == "guess") {
    select_cols = c(select_cols, "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Allele")
  } else {
    select_cols = c(select_cols, tumor_allele_col)
  }
  
  bad_maf_msg = "Input MAF is expected to be a data.frame or the filename of an MAF-formatted tab-delimited text file."
  if ("character" %in% class(maf)) {
    if(length(maf) > 1) {
      stop(bad_maf_msg)
    }
    if(file.exists(maf)) {
      # read in file and suppress warnings about missing columns since this function handles the validation
      withCallingHandlers(
      {
        maf = data.table::fread(maf, sep = "\t", quote ="", blank.lines.skip = T, select = select_cols)
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
  
  # convert all columns to character except Start_Position
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
  
  
  maf[[ref_col]] = toupper(maf[[ref_col]])
  # figure out which column has correct tumor allele data
  if(tumor_allele_col == "guess") {
    tumor_allele_col = "Tumor_Allele"
    if (tumor_allele_col %in% colnames(maf)) {
      message("Found data column \"Tumor_Allele\"; will assume this is the correct column for tumor allele.")
      maf[[tumor_allele_col]] = toupper(maf[[tumor_allele_col]])
    } else {
      # automated tumor allele determination requires Tumor_Seq_Allele1/Tumor_Seq_Allele2 columns
      # if these columns are present, the tumor_allele_adder function will handle capitalization and other validation
      allele1_col = "Tumor_Seq_Allele1"
      allele2_col = "Tumor_Seq_Allele2"
      if (! allele1_col %in% colnames(maf) | ! allele2_col %in% colnames(maf)) {
        stop(paste0("Tumor alleles can't be deduced automatically deduced without Tumor_Seq_Allele1 ",
                    "and Tumor_Seq_Allele2 columns. Please manually specify with \"tumor_allele_col=...\""))
      }
      maf$Tumor_Seq_Allele1 = toupper(maf$Tumor_Seq_Allele1)
      maf$Tumor_Seq_Allele2 = toupper(maf$Tumor_Seq_Allele2)
      message("Determining tumor alleles...")
      
      # take allele 2 as the tumor allele, but when it matches ref, replace with allele 1
      # if that is still equal to ref, record will later be discarded
      tumor_alleles = maf$Tumor_Seq_Allele2
      allele_2_matches_ref = maf$Tumor_Seq_Allele2 == maf[[ref_col]]
      tumor_alleles[allele_2_matches_ref] = maf$Tumor_Seq_Allele1[allele_2_matches_ref]
      maf[[tumor_allele_col]] <- tumor_alleles
    }
  } else {
    maf[[tumor_allele_col]] = toupper(maf[[tumor_allele_col]])
  }
  
  # drop records where tumor allele is equal to reference allele
  no_variant = maf[[ref_col]] == maf[[tumor_allele_col]]
  num_unvaried = sum(no_variant)
  if(num_unvaried > 0) {
    maf = maf[!no_variant]
    warning(paste0(num_unvaried, " MAF records had tumor alleles identical to reference alleles; these were removed from analysis.\n"))
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
    sample_progressions = rep("1", nrow(maf)) # if analysis ignores tumor progression levels, all tumors get assigned level 1
  }
  
  # add samples CESProgressions object (don't actually assign to the CESAnalysis until the end of this function)
  progressions_update = cancereffectsizeR:::add_samples_to_CESProgressions(progressions = cesa@progressions, samples = maf[[sample_col]], sample_progressions = sample_progressions)
  
  
  # select only the necessary columns and give column names that will stay consistent
  maf = maf[,c(..sample_col, ..chr_col, ..start_col, ..ref_col, ..tumor_allele_col)]
  sample_col = "Unique_Patient_Identifier"
  chr_col = "Chromosome"
  start_col = "Start_Position"
  ref_col = "Reference_Allele"
  tumor_allele_col = "Tumor_Allele"
  colnames(maf) =  c(sample_col, chr_col, start_col, ref_col, tumor_allele_col)
  
  # if CESAnalysis already contains samples, make sure the new data has no repeat samples
  if (nrow(cesa@maf) > 0) {
    previous_samples = cesa@maf$Unique_Patient_Identifier
    new_samples = maf$Unique_Patient_Identifier
    if (length(intersect(previous_samples, new_samples)) > 0) {
      stop(paste0("Error: Sample identifiers in new data have overlap with those already in the CESAnalysis.\n",
                  "If this is intentional, merge the data manually before loading it."))
    }
  }

  # create MAF df to hold excluded records
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
        message(silver("Assuming that chrM in the chain file refers to mitochondrial DNA (aka chrMT)...."))
        mesage(silver("If you get reference mismatches on this contig, you may need a different chain file."))
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
      message(silver(paste0("Note: ", num_failing, " records failed liftOver, so they will be set aside.")))
    }
    maf = merged_maf[, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
  }
  

  # get coverage info for new samples and incorporate with any existing info
  coverage_update = cesa@coverage
  if (class(covered_regions) == "character" && covered_regions[1] == "exome" && length(covered_regions) == 1) {
    if ("exome" %in% names(coverage_update$samples_by_coverage)) {
      coverage_update$samples_by_coverage[["exome"]] = c(coverage_update$samples_by_coverage[["exome"]], unique(maf[[sample_col]]))
    } else {
      coverage_update$samples_by_coverage[["exome"]] = unique(maf[[sample_col]])
    }
  } else {
    # create GRanges for panel sequencing coverage
    genome_info = GenomeInfoDb::seqinfo(cesa@genome)
    if ("data.frame" %in% class(covered_regions)) {
      grange_cols = colnames(covered_regions)[1:3]
      coverage = GenomicRanges::makeGRangesFromDataFrame(covered_regions, seqnames.field = grange_cols[1], start.field = grange_cols[2], end.field = grange_cols[3], 
                                                         seqinfo = genome_info)
    } else if (class(covered_regions) == "character") {
      coverage = read.table(file = covered_regions, comment.char = '#', sep = '\t', quote = '', blank.lines.skip = T, stringsAsFactors = F)
      coverage = GenomicRanges::makeGRangesFromDataFrame(coverage, starts.in.df.are.0based = T, seqnames.field = "V1", start.field = "V2", end.field = "V3", seqinfo = genome_info)
    } else {
      stop("Error: covered_regions should be a path to a BED-formatted file, or a data-frame-like object with chr, start (1-based), end (inclusive) in the first three columns.")
    }
    
    # combine overlapping intervals
    coverage = GenomicRanges::reduce(coverage) 
    
    # name panel coverage groups "panel_1", "panel_2", etc.
    num_groups = length(coverage_update$samples_by_coverage)
    new_coverage_group_num = ifelse("exome" %in% names(coverage_update$samples_by_coverage), num_groups, num_groups + 1)
    new_group_name = paste0("panel_", new_coverage_group_num)
    
    # update sample_by_coverage and collection of granges
    coverage_update$samples_by_coverage[[new_group_name]] = unique(maf[,Unique_Patient_Identifier])
    
    # need to put the granges into an environment because list structure couldn't handle it in testing
    if (! "granges" %in% names(coverage_update)) {
      coverage_update$granges = new.env()
    }
    coverage_update$granges[[new_group_name]] = coverage
    
    # remove any MAF records that aren't in the provided coverage
    maf_grange = GenomicRanges::makeGRangesFromDataFrame(maf, seqnames.field = chr_col, start.field = start_col, 
                                          end.field = start_col, seqinfo = genome_info)
    is_uncovered = ! maf_grange %within% coverage
    num_uncovered = sum(is_uncovered)
    
    if (num_uncovered > 0) {
      uncovered.maf = maf[is_uncovered,]
      total = nrow(maf)
      percent = round((num_uncovered / total) * 100, 1)
      uncovered.maf$Exclusion_Reason = paste0("uncovered_in_", new_group_name)
      maf = maf[!is_uncovered,]
      excluded = rbind(excluded, uncovered.maf)
      message(paste0("Warning: ", num_uncovered, " MAF records out of ", total, " (", percent, 
                     "%) are at loci not covered in the input covered_regions. ",
                     "\nThese mutations will be excluded from analysis."))
    }
  }
  
  
  
  # check for potential DNP/TNPs and separate them from data set
  # this is done before any other filtering to ensure catching as many of them as possible
  message("Searching for possible multinucleotide variants...")
  num.prefilter = nrow(maf)
  dnp_tnp_results = cancereffectsizeR::DNP_TNP_remover(maf)
  maf = dnp_tnp_results$kept[,1:5] # To-do: fix return value of DNP_TNP
  pred.mnv.maf = dnp_tnp_results$removed
  num.pred.mnv = nrow(pred.mnv.maf)
  
  if (num.pred.mnv > 0) {
    pred.mnv.maf$Exclusion_Reason = "predicted_MNV"
    excluded = rbind(excluded, pred.mnv.maf)
    # To-do: move message to DNP_TNP_remover or otherwise ensure this description remains accurate
    percent = round((num.pred.mnv / num.prefilter) * 100, 1)
    message(paste0("Note: ", num.pred.mnv, " mutation records out of ", num.prefilter, " (", percent, "%) ",
                   "are within 2 bp of other mutations in the same tumors."))
    message("These records will be excluded from effect size analysis.")
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
  message("Checking that reference alleles match reference genome...")
  ref_allele_lengths = nchar(maf[, Reference_Allele])
  ref_alleles_to_test = maf[, Reference_Allele]
  end_pos = maf[, Start_Position] + ref_allele_lengths - 1 # for multi-base deletions, check that all deleted bases match reference
  
  reference_alleles <- as.character(BSgenome::getSeq(cesa@genome, maf[,Chromosome],
                                                     strand="+", start=maf[,Start_Position], end=end_pos))
  num.prefilter = nrow(maf)
  
  # For insertions, can't evaluate if record matches reference since MAF reference will be just "-"
  # Could conceivably verify that the inserted bases don't match reference....
  reference.mismatch.maf = maf[Reference_Allele != reference_alleles & Reference_Allele != '-',]
  num.nonmatching = nrow(reference.mismatch.maf)
  maf = maf[Reference_Allele == reference_alleles | Reference_Allele == '-',] # filter out ref mismatches
  
  # if significant numbers of SNVs don't match reference, don't run
  mismatch_snv = reference.mismatch.maf[Reference_Allele %in% c("A", "C", "T", "G"), .N]
  total_snv = maf[Reference_Allele %in% c("A", "C", "T", "G"), .N]
  bad_snv_frac = mismatch_snv / total_snv
  bad_snv_percent = round(bad_snv_frac * 100, 1)
  if(bad_snv_frac > .05) {
    stop(paste0(bad_snv_percent, "% of SNV records do not match the given reference genome. Remove or fix these records before running."))
  }
  
  if (num.nonmatching > 0) {
    reference.mismatch.maf$Exclusion_Reason = "reference_mismatch"
    excluded = rbind(excluded, reference.mismatch.maf)
    percent = round((num.nonmatching / num.prefilter) * 100, 1)
    message(paste0("Note: ", num.nonmatching, " mutation records out of ", num.prefilter, " (", percent, "%, including ", bad_snv_percent,
                    "% of SNV records) have reference alleles that do not actually match the reference genome."))
    message("These records will be excluded from effect size analysis.")
  } else {
    message("Reference alleles look good.")
  }
  
  # Count SNVs that will be included in effect size analysis
  message("Collecting all SNVs...")
  num.total = nrow(maf) # indels won't be filtered, but still another filtering step later
  bases = c("A","T","G","C")
  snv.maf = maf[Reference_Allele %in% bases & Tumor_Allele %in% bases,]
  num.non.snv = num.total - nrow(snv.maf)
  if (num.non.snv > 0) {
    percent = round((num.non.snv / num.total) * 100, 1)
    message(paste0("Note: ", num.non.snv, " mutation records out of ", num.total, " (", percent, "%) are not SNVs.\n",
                   "(That is, ref or tumor allele were something other than a single A/C/T/G.)"))
    message("Indels will be run through dNdScv for your convenience, but they will not be included in SNV selection analysis.")
  }
  
  num.good.snv = nrow(snv.maf)
  
  if (num.good.snv == 0) {
    stop("Error: No SNVs are left to analyze!")
  }
  num.samples = length(unique(snv.maf[,Unique_Patient_Identifier]))
  message(paste0("Loaded ", num.good.snv, " SNVs from ", num.samples, " samples into CESAnalysis."))
  
  cesa@coverage = coverage_update
  cesa@progressions = progressions_update
  
  cesa@maf = rbind(cesa@maf, maf)
  nt = c("A", "T", "C", "G")
  snv_stats = cesa@maf[Reference_Allele %in% nt & Tumor_Allele %in% nt, .(num_samples = length(unique(Unique_Patient_Identifier)), num_snv = .N)]
  cesa@status[["MAF data"]] = paste0(snv_stats$num_snv, " SNV records from ", snv_stats$num_samples, " samples (view with maf())")
  if (nrow(excluded) > 0) {
    colnames(excluded) = c(colnames(maf), "Exclusion_Reason")
    cesa@excluded = rbind(cesa@excluded, excluded) 
  }
  return(cesa)
}

