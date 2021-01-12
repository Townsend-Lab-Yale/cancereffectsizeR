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
#' @param refset_env a refset data environment
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.table or data.frame format
#' @param sample_col column name with sample ID data (Tumor_Sample_Barcode or Unique_Patient_Identifier)
#' @param chr_col column name with chromosome data  (Chromosome)           
#' @param start_col column name with start position (Start_Position)
#' @param ref_col column name with reference allele data (Reference_Allele)
#' @param tumor_allele_col column name with alternate allele data; by default,
#'   values from Tumor_Seq_Allele2 and Tumor_Seq_Allele1 columns are used.
#' @param chain_file a LiftOver chain file (text format, name ends in .chain) to convert MAF
#'   records to the genome build used in the CESAnalysis.
#' @return data.table with core MAF columns, any other requested columns, and a "problem" column
#' @keywords internal
read_in_maf = function(maf, refset_env, chr_col = "Chromosome", start_col = "Start_Position", 
                       ref_col = "Reference_Allele", tumor_allele_col = "guess", sample_col = "Tumor_Sample_Barcode", 
                       more_cols = NULL, chain_file = NULL) {
  
  select_cols = c(sample_col, chr_col, start_col, ref_col)
  if (tumor_allele_col == "guess") {
    select_cols = c(select_cols, "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Allele")
  } else {
    select_cols = c(select_cols, tumor_allele_col)
  }
  
  # when sample_col has default value, will also check for Unique_Patient_Identifier if the default isn't found (since CESAnalysis creates this column)
  if (sample_col == "Tumor_Sample_Barcode") {
    select_cols = c(select_cols, "Unique_Patient_Identifier")
  }
  
  
  if (! is.null(more_cols)) {
    if(! is.character(more_cols)) {
      stop("more_cols should be type character (\"all\" to read all columns) or left NULL.")
    }
    
    if(length(more_cols) == 1 && more_cols == "all") {
      select_cols = NULL
    } else {
      select_cols = unique(c(select_cols, more_cols))
    }
  }
  
  

  bad_maf_msg = "Input MAF is expected to be a data frame or the filename of an MAF-formatted tab-delimited text file."
  if (is.character(maf)) {
    if(length(maf) > 1) {
      stop(bad_maf_msg)
    }
    if(file.exists(maf)) {
      # read in file and suppress warnings about missing columns since this function handles the validation
      withCallingHandlers(
        {
          if (is.null(select_cols)) {
            maf = fread(maf, sep = "\t", quote ="", blank.lines.skip = T)
          } else {
            maf = fread(maf, sep = "\t", quote ="", blank.lines.skip = T, select = select_cols)
          }
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
    maf = as.data.table(maf)
    if (! is.null(select_cols)) {
      cols_to_keep = names(maf)[which(names(maf) %in% select_cols)]
      maf = maf[, ..cols_to_keep]
    }
  }
  missing_cols = character()
  input_maf_cols = colnames(maf)
  
  if (sample_col == "Tumor_Sample_Barcode" & ! sample_col %in% input_maf_cols & "Unique_Patient_Identifier" %in% input_maf_cols) {
    sample_col = "Unique_Patient_Identifier"
    pretty_message("Found column Unique_Patient_Identifier; we'll assume this is the correct sample ID column.")
  }
  cols_to_check = c(sample_col, ref_col, chr_col, start_col)
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
  
  
  setnames(maf, c(chr_col, start_col, ref_col, tumor_allele_col, sample_col),
           c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele",
             "Unique_Patient_Identifier"), skip_absent = TRUE) # tumor_allele_col will be absent when set to "guess"
  maf[, c("Chromosome", "Unique_Patient_Identifier") := list(as.character(Chromosome),
                                                             as.character(Unique_Patient_Identifier))]
  maf[, Start_Position := as.numeric(Start_Position)]
  maf[, Reference_Allele := toupper(Reference_Allele)]
  
  # figure out which column has correct tumor allele data
  maf_cols = colnames(maf)
  if(tumor_allele_col == "guess") {
    tumor_allele_col = "Tumor_Allele"
    if (! tumor_allele_col %in% maf_cols) {
      # automated tumor allele determination requires Tumor_Seq_Allele1/Tumor_Seq_Allele2 columns
      allele1_col = "Tumor_Seq_Allele1"
      allele2_col = "Tumor_Seq_Allele2"
      
      if ("Tumor_Seq_Allele1" %in% maf_cols && ! "Tumor_Seq_Allele2" %in% maf_cols) {
        pretty_message("Found column Tumor_Seq_Allele1 but not Tumor_Seq_Allele2; we'll assume Tumor_Seq_Allele1 is the correct tumor allele column.")
        setnames(maf, "Tumor_Seq_Allele1", "Tumor_Allele")
        maf[, Tumor_Allele := toupper(Tumor_Allele)]
      } else if ("Tumor_Seq_Allele2" %in% maf_cols && ! "Tumor_Seq_Allele1" %in% maf_cols) {
        pretty_message("Found column Tumor_Seq_Allele2 but not Tumor_Seq_Allele1;\nwe'll assume Tumor_Seq_Allele2 is the correct tumor allele column.")
        setnames(maf, "Tumor_Seq_Allele2", "Tumor_Allele")
        maf[, Tumor_Allele := toupper(Tumor_Allele)]
      } else if (! any(c("Tumor_Seq_Allele1", "Tumor_Seq_Allele2") %in% maf_cols)) {
        stop("Tumor alleles column not found. Please specify with \"tumor_allele_col=...\"")
      } else {
        # take allele 2 as the tumor allele, but when it matches ref, replace with allele 1
        # if that is still equal to ref, record will later be discarded
        maf[, Tumor_Allele := toupper(Tumor_Seq_Allele2)]
        maf[Tumor_Allele == Reference_Allele, Tumor_Allele := toupper(Tumor_Seq_Allele1)]
        maf[, c("Tumor_Seq_Allele1", "Tumor_Seq_Allele2") := NULL]
      }
    }
  } else {
    maf[, Tumor_Allele := toupper(Tumor_Allele)]
  }
  
  maf_cols = c("Unique_Patient_Identifier", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
  maf[, problem := NA_character_]
  looks_empty = apply(maf[, ..maf_cols], 1, function(x) any(is.na(x)) | any(grepl("^[.\"\' ]*$", x)))
  num_bad = sum(looks_empty)
  maf[looks_empty, problem := 'missing_values']
  
  # Problem if tumor allele matches reference allele
  no_variant = maf[Reference_Allele == Tumor_Allele, which = T]
  maf[no_variant, problem := 'not_variant']
  
  duplicate_records = duplicated(maf[,.(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele)])
  maf[duplicate_records, problem := 'duplicate_record']
  
  # run liftOver if chain file supplied
  if(! is.null(chain_file)) {
    message("Preparing and running liftOver...")
    chain = rtracklayer::import.chain(chain_file)
    names(chain) = sub("^chr", "", names(chain))
    prelift_cols = copy(names(maf))
    maf[, rn := 1:.N] # using row number as an identifier to know which intervals fail liftover
    for_liftover = maf[is.na(problem), .(Chromosome, Start_Position, rn)]
    for_liftover[, Chromosome := sub('^chr', '', Chromosome)]
    gr = GenomicRanges::makeGRangesFromDataFrame(df = for_liftover, seqnames.field = "Chromosome", 
                                                 start.field = "Start_Position", end.field = "Start_Position",
                                                 keep.extra.columns = T)
    lifted_over = unlist(rtracklayer::liftOver(gr, chain))
    GenomeInfoDb::seqlevelsStyle(lifted_over) = "NCBI"
    lifted_over = as.data.table(lifted_over)
    
    merged_maf = merge.data.table(maf, lifted_over, by = "rn")
    merged_maf[, Chromosome := as.character(seqnames)] # seqnames comes back as factor!
    merged_maf[, Start_Position := start]
    
    failed_liftover = maf[is.na(problem) & ! rn %in% merged_maf$rn, -"rn"]
    maf = merged_maf[, ..prelift_cols]
    if (failed_liftover[, .N] > 0) {
      failed_liftover[, problem := 'failed_liftOver']
      maf = rbind(maf, failed_liftover)
    }
    
    # Different loci in one genome build can get lifted to the same position in the
    # another, due to fixes. Therefore, liftOver sometimes reveals duplicate records
    # when two mutations in same sample get lifted to the same locus in a newer build.
    lifted_to_same = duplicated(maf[,.(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele)])
    maf[lifted_to_same, problem := 'duplicate_record_after_liftOver']
  }
  
  # strip chr prefixes from chr column, if present ("NCBI-style"), and flag unsupported chromosomes
  supported_chr = refset_env$supported_chr
  maf[, supported := FALSE][Chromosome %in% supported_chr, supported := T]
  maf[supported == FALSE, stripped_chr := sub('^chr', '', Chromosome) ]
  maf[supported == FALSE & stripped_chr %in% supported_chr, c("Chromosome", "supported") := list(stripped_chr, TRUE)]
  
  maf[supported == FALSE, problem := 'unsupported_chr']
  maf[, c("stripped_chr", "supported") := NULL]
  
  
  # Ensure reference alleles of mutations match reference genome (Note: Insertions won't match if their reference allele is "-")
  message("Checking that reference alleles match the reference genome...")
  ref_allele_lengths = nchar(maf[is.na(problem), Reference_Allele])
  ref_alleles_to_test = maf[is.na(problem), Reference_Allele]
  end_pos = maf[is.na(problem), Start_Position] + ref_allele_lengths - 1 # for multi-base deletions, check that all deleted bases match reference
  reference_alleles <- as.character(BSgenome::getSeq(refset_env$genome, maf[is.na(problem), Chromosome],
                                                     strand="+", start=maf[is.na(problem), Start_Position], end=end_pos))
  
  # For insertions, can't evaluate if record matches reference since MAF reference will be just "-"
  # Could conceivably verify that the inserted bases don't match reference....
  maf[is.na(problem), actual_ref := reference_alleles]
  maf[is.na(problem) & Reference_Allele != '-' & Reference_Allele != actual_ref, problem := 'reference_mismatch']
  maf[, actual_ref := NULL]
  return(maf)
}
  