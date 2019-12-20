#' load_maf
#' @description Load MAF-formatted mutation data into a CESAnalysis object (the core cancereffectsizeR object)
#' @importFrom IRanges "%within%"
#' @param cesa the CESAnalysis object to load the data into
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.frame format
#' @param sample_col column name with sample ID data
#' @param chr_col column name with chromosome data             
#' @param start_col column name with start position ( currently only analyze SNVs, so no end position needed)
#' @param ref_col column name with reference allele data
#' @param tumor_allele_col column name with alternate allele data (defaults to "guess", which examines ref_col,
#'                         "Tumor_Seq_Allele1", and "Tumor_Seq_Allele2" and determines from there)
#' @param covered_regions for panel sequencing data, a BED file of covered interals, or a data frame with first three columns chr, start (1-based), end (inclusive)
#' @param progression_col column in MAF with subset data (e.g., column contains data like "Primary" and "Metastatic" in each row)
#' @return CESAnalysis object with somewhat cleaned-up MAF data added 
#' @export
load_maf = function(cesa = NULL, maf = NULL, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", covered_regions = "exome",
                    progression_col = NULL) {
  if (is.null(cesa)) {
    stop("You need to supply a CESAnalysis object to load the MAF data into.")
  }
  
  if (is.null(maf)) {
    stop("Supply MAF data via maf=[file path or data frame].")
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
  
  bad_maf_msg = "Input MAF is expected to be a data.frame or the filename of an MAF-formatted tab-delimited text file."
  if ("character" %in% class(maf)) {
    if(length(maf) > 1) {
      stop(bad_maf_msg)
    }
    if(file.exists(maf)) {
      maf = tryCatch(
        {
          # To-do: Switch to fread to make reading large files faster,
          # but need to handle optional presence of Tumor_Seq_Allele_1/2 columns
          #data.table::fread(maf, sep = "\t", header=T, stringsAsFactors = F, quote="", select = data.c("hey"))
          data.table::fread(maf, sep = "\t", quote ="", blank.lines.skip = T)
        },
        error = function(e) {
          message("Unable to read specified MAF file:")
          message(e)
        }
      )
    } else {
      stop("Specified MAF not found.")
    }
  }
  if(! "data.frame" %in% class(maf)) {
    stop(bad_maf_msg)
  } else if(nrow(maf) == 0) {
    stop("Input MAF data set is empty.")
  }
  
  missing_cols = character()
  input_maf_cols = colnames(maf)
  
  
  cols_to_check = c(sample_col, ref_col, chr_col, start_col)
  if (tumor_allele_col != "guess") {
    cols_to_check = c(cols_to_check, tumor_allele_col)
  }
  if (! is.null(progression_col)) {
    cols_to_check = c(cols_to_check, progression_col)
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
  
  
  # figure out which column has correct tumor allele data
  if(tumor_allele_col == "guess") {
    if ("Tumor_Allele" %in% colnames(maf)) {
      message("Found data column \"Tumor_Allele\"; will assume this is the correct column for tumor allele.")
      tumor_allele_col = "Tumor_Allele"
    } else {
      # automated tumor allele determination requires Tumor_Seq_Allele1/Tumor_Seq_Allele2 columns
      # if these columns are present, the tumor_allele_adder function will handle capitalization and other validation
      allele1_col = "Tumor_Seq_Allele1"
      allele2_col = "Tumor_Seq_Allele2"
      if (! allele1_col %in% colnames(maf) | ! allele2_col %in% colnames(maf)) {
        stop(paste0("Tumor alleles can't be deduced automatically deduced without Tumor_Seq_Allele1",
                    "and Tumor_Seq_Allele2 columns. Please manually specify with \"tumor_allele_col=...\""))
      }
      message("Determining tumor alleles...")
      tumor_allele_adder_returns <- cancereffectsizeR::tumor_allele_adder(Reference_Allele = maf[,ref_col],
                                                                          Tumor_Seq_Allele1 = maf[,allele1_col],
                                                                          Tumor_Seq_Allele2 =  maf[,allele2_col])
      maf <- maf[tumor_allele_adder_returns$original_index,]
      maf[,tumor_allele_col] <- tumor_allele_adder_returns$Tumor_allele
    }
  }
  maf[,ref_col] = toupper(maf[,ref_col])
  maf[,tumor_allele_col] = toupper(maf[,tumor_allele_col])
  
  # collect tumor progression information
  if (! is.null(progression_col)) {
    if (is.factor(maf[,progression_col])) {
      message("Warning: You supplied tumor progression as a factor, but it will be converted to character.")
      message("The progression ordering will be what you supplied to CESAnalysis() with \"progression_order\",")
      message("regardless of any ordering in your input MAF data frame.")
    }
    sample_progressions = as.character(maf[,progression_col])
    if(any(is.na(sample_progressions))) {
      stop("Error: There are NA values in your sample progressions column.")
    }
  } else {
    sample_progressions = rep("1", nrow(maf)) # if analysis ignores tumor progression levels, all tumors get assigned level 1
  }
  
  
  # select only the necessary columns and give column names matching MAFdf format
  maf = maf[,c(sample_col, chr_col, start_col, ref_col, tumor_allele_col)]
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
  excluded = data.frame()
  
  
  # add samples CESProgressions object (don't actually assign to the CESAnalysis until the end of this function)
  progressions_update = cancereffectsizeR:::add_samples_to_CESProgressions(progressions = cesa@progressions, samples = maf[,sample_col], sample_progressions = sample_progressions)
  
  
  # strip chr prefixes from chr column, if present ("NCBI-style")
  maf[,chr_col] = gsub('^chr', '', maf[,chr_col])
  
  # temporary handling of MT/M
  maf[maf[,chr_col] == "M", chr_col] <- "MT" # changing M to MT
  
  # get coverage info for new samples and incorporate with any existing info
  coverage_update = cesa@coverage
  if (class(covered_regions) == "character" && covered_regions[1] == "exome" && length(covered_regions) == 1) {
    if ("exome" %in% names(coverage_update$samples_by_coverage)) {
      coverage_update$samples_by_coverage[["exome"]] = c(coverage_update$samples_by_coverage[["exome"]], unique(maf[,sample_col]))
    } else {
      coverage_update$samples_by_coverage[["exome"]] = unique(maf[,sample_col])
    }
  } else {
    # create GRanges for panel sequencing coverage
    genome_info = GenomeInfoDb::Seqinfo(cesa@genome)
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
    coverage_update$samples_by_coverage[[new_group_name]] = unique(maf[, sample_col])
    
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
  is_mt <- maf[,chr_col] == "MT"
  if (any(is_mt)) {
    mt_maf = maf[is_mt,]
    mt_maf$Exclusion_Reason = "mitochondrial"
    excluded = rbind(excluded, mt_maf)
    maf <- maf[! is_mt,]
    message(paste("Note:", sum(is_mt), "mitochondrial mutations have been removed from the data (mitochondrial analysis not yet supported)."))
  }
  
  
  # Ensure reference alleles of mutations match reference genome (Note: Insertions won't match if their reference allele is "-")
  message("Checking that reference alleles match reference genome...")
  ref_allele_lengths = nchar(maf[,ref_col])
  ref_alleles_to_test = maf[,ref_col]
  end_pos = maf[,start_col] + ref_allele_lengths - 1 # for multi-base deletions, check that all deleted bases match reference
  
  reference_alleles <- as.character(BSgenome::getSeq(cesa@genome, maf[,chr_col],
                                                     strand="+", start=maf[,start_col], end=end_pos))
  num.prefilter = nrow(maf)
  
  # For insertions, can't evaluate if record matches reference since MAF reference will be just "-"
  # Could conceivably verify that the inserted bases don't match reference....
  reference.mismatch.maf = maf[maf[,ref_col] != reference_alleles & maf[,ref_col] != '-',]
  num.nonmatching = nrow(reference.mismatch.maf)
  maf = maf[maf[,ref_col] == reference_alleles | maf[,ref_col] == '-',] # filter out ref mismatches
  
  if (num.nonmatching > 0) {
    reference.mismatch.maf$Exclusion_Reason = "reference_mismatch"
    excluded = rbind(excluded, reference.mismatch.maf)
    percent = round((num.nonmatching / num.prefilter) * 100, 1)
    message(paste0("Note: ", num.nonmatching, " mutation records out of ", num.prefilter, " (", percent, "%) ",
                   "have reference alleles that do not actually match the reference genome."))
    message("These records will be excluded from effect size analysis.")
  } else {
    message("Reference alleles look good.")
  }
  
  # Count SNVs that will be included in effect size analysis
  message("Collecting all SNVs...")
  num.total = nrow(maf) # indels won't be filtered, but still another filtering step later
  bases = c("A","T","G","C")
  snv.maf = maf[maf[,ref_col] %in% bases  & maf[,tumor_allele_col] %in% bases,]
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
  num.samples = length(unique(snv.maf[, sample_col]))
  message(paste0("Loaded ", num.good.snv, " SNVs from ", num.samples, " samples into CESAnalysis."))
  
  cesa@coverage = coverage_update
  cesa@progressions = progressions_update 
  cesa@maf = MAFdf(rbind(cesa@maf, maf))
  if (nrow(excluded) > 0) {
    colnames(excluded) = c(colnames(maf), "Exclusion_Reason")
    cesa@excluded = MAFdf(rbind(cesa@excluded, excluded))  
  }
  return(cesa)
}

