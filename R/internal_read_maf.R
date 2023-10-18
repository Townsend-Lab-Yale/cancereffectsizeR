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
#' @param separate_old_problems When TRUE (as used by load_maf), respect old problems that
#'   look like they came from cancereffectsizeR (typically from preload_maf). These get
#'   separated as "old_problem", and the records won't be checked. chain_file must be
#'   NULL.
#' @return data.table with core MAF columns, any other requested columns, and a "problem" column
#' @keywords internal
read_in_maf = function(maf, refset_env, chr_col = "Chromosome", start_col = "Start_Position", 
                       ref_col = "Reference_Allele", tumor_allele_col = "guess", sample_col = "Unique_Patient_Identifier", 
                       more_cols = NULL, chain_file = NULL, separate_old_problems = FALSE) {
  
  select_cols = c(sample_col, chr_col, start_col, ref_col, 'problem')
  if (tumor_allele_col == "guess") {
    select_cols = c(select_cols, "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Allele")
  } else {
    select_cols = c(select_cols, tumor_allele_col)
  }
  
  # when sample_col has default value, will also check for Tumor_Sample_Barcode if the default isn't found
  if (sample_col == "Unique_Patient_Identifier") {
    select_cols = c(select_cols, "Tumor_Sample_Barcode")
  }
  
  if (! is.null(more_cols)) {
    if(! is.character(more_cols)) {
      stop("more_cols specified incorrectly (internally)")
    }
    
    if(length(more_cols) == 1 && more_cols == "all") {
      select_cols = NULL
    } else {
      select_cols = unique(c(select_cols, more_cols))
    }
  }
  
  if (! is.null(chain_file)) {
    if(!is(chain_file, "character") || length(chain_file) != 1 || !(endsWith(chain_file, ".chain"))) {
      stop("Argument chain_file expected to be the filename/path of a text file ending in .chain")
    }
    if(! file.exists(chain_file)) {
      stop("The chain_file specified could not be found; check the file path.")
    }
  }
  
  # If it looks like this function's liftOver functionality has previously been used
  # on this data, we'll keep the columns that it generated.
  if (is.null(chain_file) && ! is.null(select_cols)) {
    if (! is.null(select_cols))
    select_cols = unique(c(select_cols, "prelift_chr", "prelift_start", "liftover_strand_flip", "prelift_variant_id"))
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
          # readChar(maf, 2) == '##' could indicate a VCF file
          # fread(maf, skip = '#CHROM') would skip to end of VCF headers
          if (is.null(select_cols)) {
            maf = fread(maf, quote ="", blank.lines.skip = T)
          } else {
            maf = fread(maf, quote ="", blank.lines.skip = T, select = select_cols)
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
  
  if (sample_col == "Unique_Patient_Identifier" && ! sample_col %in% input_maf_cols &&
      "Tumor_Sample_Barcode" %in% input_maf_cols) {
    sample_col = "Tumor_Sample_Barcode"
    pretty_message("Found column Tumor_Sample_Barcode; we'll assume that this column identifies patients.")
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
  suppressWarnings(maf[, c("Tumor_Seq_Allele1", "Tumor_Seq_Allele2") := NULL])
  
  maf_cols = c("Unique_Patient_Identifier", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
  
  old_problems = data.table()
  if('problem' %in% names(maf)) {
    if(separate_old_problems == TRUE) {
      if(length(setdiff(unique(maf$problem), c(preload_problems, 'NA', NA))) == 0) {
        old_problems = maf[! is.na(problem)]
        maf = maf[is.na(problem)]
        
        # Some excluded records should retain variant IDs
        identify_maf_variants(old_problems)
        problems_that_retain_id = c('duplicate_record', 'duplicate_from_TCGA_sample_merge',
                                    'duplicate_record_after_liftOver',
                                    'merged_into_dbs_variant', 'merged_with_nearby_variant', 
                                    'merged_into_other_variant')
        old_problems[! problem %in% problems_that_retain_id, 
                     c('variant_id', 'variant_type') := list(NA_character_, NA_character_)]
        if(! is.null(chain_file)) {
          pretty_message("Note: Records with problems identified in previous runs of preload_maf() won't be run through liftOver.")
          old_problems[problem %in% problems_that_retain_id, prelift_variant_id := variant_id]
          old_problems
        }
        setnames(old_problems, 'problem', 'old_problem')
        old_problems[, problem := NA_character_]
      } else {
        msg = paste0("Found a \"problem\" column in MAF, but if the table came from cancereffectsizeR, ",
                     "it's been subsequently altered, so it can't be used and will be overwritten.")
      }
    }
    maf[, problem := NULL]
  } 
  maf[, problem := NA_character_]
  check = maf[, lapply(.SD, function(x) is.na(x) | grepl("^[.\"\' ]*$", x)), .SDcols = maf_cols]
  looks_empty = apply(check, 1, any)
  num_bad = sum(looks_empty)
  maf[looks_empty, problem := 'missing_values']
  
  duplicate_records = duplicated(maf[,.(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)])
  maf[duplicate_records, problem := 'duplicate_record']
  
  # Change duplicate record annotation for TCGA patients with multiple (essentially replicate) samples
  if('Tumor_Sample_Barcode' %in% names(maf)) {
    if(sample_col == 'Unique_Patient_Identifier') {
      maf[, within_sample_dup := FALSE]
      sample_barcode_dups = duplicated(maf[,.(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele)])
      maf[sample_barcode_dups, within_sample_dup := TRUE]
      maf[, is_tcga_patient := Unique_Patient_Identifier %like% '^TCGA-.{7}$']
      maf[within_sample_dup == FALSE & problem == 'duplicate_record' & is_tcga_patient == TRUE, problem := 'duplicate_from_TCGA_sample_merge']
      maf[, c('within_sample_dup', 'is_tcga_patient') := NULL]
    }
    maf[, Tumor_Sample_Barcode := NULL]
  }
  
  # If MAF has unsupported chromosome names, try stripping chr prefixes. We'll run this again after liftOver.
  supported_chr = refset_env$supported_chr
  deal_with_chr_prefix = function(maf) {
    maf[, supported := Chromosome %in% supported_chr]
    maf[supported == FALSE, stripped_chr := sub('^chr', '', Chromosome) ]
    maf[supported == FALSE & stripped_chr %in% supported_chr, c("Chromosome", "supported") := list(stripped_chr, TRUE)]
    maf[supported == FALSE & is.na(problem), problem := 'unsupported_chr']
    maf[, c("stripped_chr", "supported") := NULL]
  }
  deal_with_chr_prefix(maf)
  identify_maf_variants(maf)
  
  # run liftOver if chain file supplied
  if(! is.null(chain_file)) {
    message("Preparing and running liftOver...")
    chain = rtracklayer::import.chain(chain_file)
    
    # We'll actually use UCSC style for liftOver since we don't know what style the input
    # MAF is, and it appears more consistent than NCBI.
    suppressWarnings({seqlevelsStyle(chain) = 'UCSC'})
    names(chain) = seqnames(seqinfo(chain)) # names should change with seqlevelsStyle, but they don't
    maf[, rn := 1:.N] # using row number as an identifier to know which intervals fail liftover
    for_liftover = maf[is.na(problem), .(Chromosome, Start_Position, rn)]
    
    # Assume positive strand (if there are ever serious MAF files with minus strand records, may need to change)
    gr = GenomicRanges::makeGRangesFromDataFrame(df = for_liftover, seqnames.field = "Chromosome", 
                                                 start.field = "Start_Position", end.field = "Start_Position",
                                                 keep.extra.columns = T)
    
    GenomeInfoDb::seqlevelsStyle(gr) = 'UCSC'
    GenomicRanges::strand(gr) = "+"
    lifted_over = unlist(rtracklayer::liftOver(gr, chain))
    genome(lifted_over) = genome(refset_env$genome)[1]
    suppressWarnings({GenomeInfoDb::seqlevelsStyle(lifted_over) = seqlevelsStyle(refset_env$genome)})
    lifted_over = as.data.table(lifted_over)
    
    merged_maf = merge.data.table(maf, lifted_over, by = "rn")
    merged_maf[, liftover_strand_flip := FALSE][strand == "-", liftover_strand_flip := TRUE]
    merged_maf[, prelift_variant_id := variant_id]
    setnames(merged_maf, c("Chromosome", "Start_Position"), c("prelift_chr", "prelift_start"))
    merged_maf[, Chromosome := as.character(seqnames)] # seqnames comes back as factor!
    merged_maf[, Start_Position := start]
    merged_maf[, c("start", "seqnames", "width", "strand", "end") := NULL] 
    
    # A few records (e.g., ~0.3% when lifting TCGA UCEC from hg38 to hg19) will need
    # positions/alleles corrected because of the orientation of the local contig being
    # reversed across genome builds. All of these basically get the same treatment: flip
    # start position to opposite end of variant, which means (start - length + 1), and
    # reverse complement the reference and alternate alleles. Indels have a '-' in either
    # ref or alt columns, which don't get changed. Since reverseComplement returns
    # identity on '-', we don't actually have to worry about this.
    merged_maf[liftover_strand_flip == TRUE, 
               c("Start_Position", "Reference_Allele", "Tumor_Allele") := 
                 .(Start_Position - nchar(Reference_Allele) + 1,
                   as.character(Biostrings::reverseComplement(DNAStringSet(Reference_Allele))),
                   as.character(Biostrings::reverseComplement(DNAStringSet(Tumor_Allele))))]
    deal_with_chr_prefix(merged_maf) # in case liftOver brought seqlevelsstyle out of alignment with refset
    identify_maf_variants(merged_maf) # update variant_id

    failed_liftover = maf[is.na(problem) & ! rn %in% merged_maf$rn, -"rn"]
    maf[, rn := NULL]
    setcolorder(merged_maf, names(maf))
    maf = merged_maf[, -"rn"]
    if (failed_liftover[, .N] > 0) {
      failed_liftover[, problem := 'failed_liftOver']
      failed_liftover[, c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele") := .(NA_character_, NA_real_, NA_character_, NA_character_)]
      failed_liftover[, c("prelift_chr", "prelift_start") := list(Chromosome, Start_Position)]
      failed_liftover[, prelift_variant_id := variant_id]
      failed_liftover[, variant_id := NA_character_]
      maf = rbind(maf, failed_liftover, fill = T)
    }
    
    # Different loci in one genome build can get lifted to the same position in the
    # another, due to fixes. Therefore, liftOver sometimes reveals duplicate records
    # when two mutations in same sample get lifted to the same locus in a newer build.
    lifted_to_same = duplicated(maf[,.(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele)])
    maf[lifted_to_same, problem := 'duplicate_record_after_liftOver']
  } else {
    # If liftOver is not being used (which would catch problems earlier), make sure all records are in-bounds.
    max_lengths = GenomeInfoDb::seqlengths(refset_env$genome)
    maf[Start_Position > max_lengths[Chromosome], problem := 'out_of_bounds']
  }
  
  # Ensure reference alleles of mutations match reference genome (Note: Insertions won't match if their reference allele is "-")
  # Can safely skip the reference allele check if the table has previously been checked and still has the same ref alleles
  skip_ref_check = FALSE
  prev_md5 = attr(maf, 'ref_md5')
  if(! is.null(prev_md5)) {
    prev_md5_noproblem = attr(maf, 'ref_md5_noproblem')
    current_md5 = digest::digest(maf[, .(variant_id, Reference_Allele)])
    if (current_md5  %in% c(prev_md5, prev_md5_noproblem)) {
      skip_ref_check = TRUE
    }
  }

  if(! skip_ref_check) {
    message("Checking that reference alleles match the reference genome...")
    ref_allele_lengths = nchar(maf[is.na(problem), Reference_Allele])
    ref_alleles_to_test = maf[is.na(problem), Reference_Allele]
    end_pos = maf[is.na(problem), Start_Position] + ref_allele_lengths - 1 # for multi-base deletions, check that all deleted bases match reference
    reference_alleles = as.character(BSgenome::getSeq(refset_env$genome, maf[is.na(problem), Chromosome],
                                                      strand="+", start=maf[is.na(problem), Start_Position], end=end_pos))
    
    # For insertions, can't evaluate if record matches reference since MAF reference will be just "-".
    # "Other" variants will also not be evaluated.
    maf[is.na(problem), actual_ref := reference_alleles]
    maf[is.na(problem) & variant_type %in% c('del', 'sbs', 'dbs') & Reference_Allele != actual_ref, 
        c("problem", "variant_type", "variant_id") := list('reference_mismatch', NA_character_, NA_character_)]
    maf[, actual_ref := NULL]
  }
  
  # Problem if tumor allele matches reference allele
  no_variant = maf[is.na(problem) & Reference_Allele == Tumor_Allele, which = T]
  maf[no_variant, c('problem', 'variant_id', 'variant_type') := list('not_variant', NA_character_, NA_character_)]
  
  maf = rbind(maf, old_problems, fill = T)
  setcolorder(maf, maf_cols)
  return(maf)
}

#' Annotate MAF data with variant types and variant IDs
#' 
#' @param maf a validated data.table with MAF-like columns
#' @return input table with variant_type and variant_id columns
#' @keywords internal
identify_maf_variants = function(maf) {
  suppressWarnings(maf[, c("variant_id", "variant_type") := NULL])
  # Format Start_Position carefully, to avoid issues like position 100000 getting rendered as "1e5"
  maf[, start_char := format(Start_Position, scientific = F, justify = 'none', trim = T)]
  maf[Reference_Allele %ilike% '^[ACTG]$' & Tumor_Allele %ilike% '^[ACTG]$', 
      c("variant_type", "variant_id") := .("sbs", paste0(Chromosome, ':', start_char, '_', toupper(Reference_Allele), '>', toupper(Tumor_Allele)))]
  maf[Reference_Allele %ilike% '^[ACTG]{2}$' & Tumor_Allele %ilike% '^[ACTG]{2}$', 
      c("variant_type", "variant_id") := .('dbs', paste0(Chromosome, ':', start_char, '_', toupper(Reference_Allele), '>', toupper(Tumor_Allele)))]
  maf[Reference_Allele == '-' & Tumor_Allele %ilike% '^[ACTG]+$', 
      c("variant_type", "variant_id") := .('ins', paste0(Chromosome, ':', start_char, '_ins_', toupper(Tumor_Allele)))]
  maf[Tumor_Allele == '-' & Reference_Allele %like% '^[ACTG]+$', 
      c("variant_type", "variant_id") := .('del', paste0(Chromosome, ':', start_char, '_del_', nchar(Reference_Allele)))]
  maf[is.na(variant_id), c("variant_type", "variant_id") := .('other', paste0(Chromosome, ':', start_char, '_', toupper(Reference_Allele), 
                                                                              '>', toupper(Tumor_Allele)))]
  
  # "complex" variants result from MNV detection and will get validated separately
  complex_variants_ind = maf[variant_type == 'other' & Reference_Allele %like% ',' & Tumor_Allele %like% ',', which = T]
  maf[, is_complex := FALSE][complex_variants_ind, is_complex := TRUE]
  
  # Don't allow commas at start or end, and require matched number in ref/tumor.
  maf[is_complex == TRUE & (Reference_Allele %like% '(^,)|(,$)' | Tumor_Allele %like% '(^,)|(,$)' | 
                              nchar(gsub('[^,]', '', Reference_Allele)) != 
                              nchar(gsub('[^,]', '', Tumor_Allele))),
                            c("is_complex", "variant_type", "variant_id") :=
                              .(FALSE, 'illegal', NA_character_)]
                            
  if(length(complex_variants_ind) > 0) {
    complex_variants = maf[complex_variants_ind]
    complex_variants[, rn := complex_variants_ind]
    maf[is_complex == TRUE, rn := complex_variants_ind]
    
    tmp = copy(complex_variants)
    to_identify = complex_variants[, .(Chromosome, Start_Position, Reference_Allele, Tumor_Allele,
                                       Unique_Patient_Identifier, rn)][0]
    
    # iteratively process comma-delimited variant IDs from left to right
    while(tmp[, .N] > 0) {
      left_variant = tmp[, .(Chromosome,
                             Start_Position,
                             Reference_Allele = gsub(',.*', '', Reference_Allele),
                             Tumor_Allele = gsub(',.*', '', Tumor_Allele),
                             Unique_Patient_Identifier, rn)]
      
      left_variant[, to_add := 0]
      left_variant[Reference_Allele %like% '\\(\\+\\d+\\)',  to_add := as.numeric(gsub('.*\\(\\+(\\d+).*', '\\1', Reference_Allele))]
      left_variant[, Reference_Allele := gsub('.*\\)', '', Reference_Allele)]
      left_variant[, Start_Position := Start_Position + to_add]
      left_variant[, to_add := NULL]
      
      to_identify = rbind(to_identify, left_variant)
      
      tmp = tmp[Reference_Allele %like% ',', 
                .(Chromosome,
                  Start_Position,
                  Reference_Allele = gsub('.*?,', '', Reference_Allele),
                  Tumor_Allele = gsub('.*?,', '', Tumor_Allele),
                  Unique_Patient_Identifier, rn)]
    }
    
    to_identify = identify_maf_variants(to_identify)
    bad_complex_rn = to_identify[variant_type == 'illegal', unique(rn)]
    maf[bad_complex_rn, c("variant_type", "variant_id") := .('illegal', NA_character_)]
    good_ids_by_rn = to_identify[! bad_complex_rn, .(variant_id = paste0(variant_id, collapse = ',')), by = 'rn']
    maf[good_ids_by_rn, variant_id := i.variant_id, on = 'rn']
    maf[, rn := NULL]
  }
  
  
  # Variants that don't match sbs, DBS, indel formats are either MNV substitutions,
  # or length-changing substitutions (e.g., CT->AAG), or should only contain ACTG.
  maf[variant_type == 'other' & ! is_complex & (! Reference_Allele %ilike% '^[ACTG]+$' | 
                                                 ! Tumor_Allele %ilike% '^[ACTG]+$'),
      c("variant_type", "variant_id") := .('illegal', NA_character_)]
  
  maf[, is_complex := NULL]
  maf[, start_char := NULL]
  return(maf)
}

  