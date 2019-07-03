#' Reads in an MAF file or data frame and creates a CESAnalysis object (the data structure used by cancereffectsizeR)
#' @param maf	Path of tab-delimited text file in MAF format, or an MAF in data.frame format
#' @param sample_col column name with sample ID data
#' @param ref_col column name with reference allele data
#' @param tumor_allele_col column name with alternate allele data (defaults to "guess", which examines ref_col,
#'                "Tumor_Seq_Allele1", and "Tumor_Seq_Allele2" and determines from there)
#' @param start_col column name with start position ( currently only analyze SNVs, so no end position needed)
#' @param chr_col column name with chromosome data
#' @param progression_col column in MAF with subset data (e.g., column contains data like "Primary" and "Metastatic" in each row)
#' @param progression_order evolutionary order of events in progression_col. (e.g. c("Primary", "Metastatic")
#' @export


CESAnalysis = function(maf = NULL, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                   ref_col = "Reference_Allele", tumor_allele_col = "guess",
                   progression_col = NULL, progression_order = NULL) {

	if (xor(is.null(progression_col), is.null(progression_order))) {
	  stop("progression_col and progression_order must be defined together, or left null for analyses that don't distinguish tumor progression")
	}
	if (is.null(maf)) {
	  stop("maf must be defined")
	}

  
	bad_maf_msg = "Input MAF is expected to be a data.frame or the filename of an MAF-formatted tab-delimited text file."
	if (class(maf) == "character") {
		if(length(maf) > 1) {
			stop(bad_maf_msg)
		}
		if(file.exists(maf)) {
		    maf = tryCatch(
		      {
		        # To-do: Switch to fread to make reading large files faster, 
		        # but need to handle optional presence of Tumor_Seq_Allele_1/2 columns
		        #data.table::fread(maf, sep = "\t", header=T, stringsAsFactors = F, quote="", select = data.c("hey"))
		        read.table(maf, sep = "\t", header = T, stringsAsFactors = F, quote ="", comment.char = '#')
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
	if(class(maf) != "data.frame") {
		stop(bad_maf_msg)
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
	  msg = "The following required MAF columns couldn't be found:"
	  msg = paste(msg, missing_cols, sep = "\n\t")
	  stop(msg)
	}
	
	
	# To-do: move this to a locus/build validation function
	maf[,chr_col] = gsub('MT', 'M', maf[,chr_col])

	if(tumor_allele_col == "guess") {
	  # automated tumor allele determination requires Tumor_Seq_Allele1/Tumor_Seq_Allele2 columns
	  # if these columns are present, the tumor_allele_adder function will handle capitalization and other validation
    allele1_col = "Tumor_Seq_Allele1" 
    allele2_col = "Tumor_Seq_Allele2"
	  if (! allele1_col %in% colnames(maf) | ! allele2_col %in% colnames(maf)) {
	    stop(paste0("Tumor alleles can't be deduced automatically deduced without Tumor_Seq_Allele1",
	                "and Tumor_Seq_Allele2 columns. Please manually specify with \"tumor_allele_col=...\""))
	  }
    message("Determining tumor alleles...")
    maf[,tumor_allele_col] = cancereffectsizeR::tumor_allele_adder(maf[,ref_col], maf[,allele1_col], maf[,allele2_col])
	} else {
	  maf[,ref_col] = toupper(maf[,ref_col])
	  maf[,tumor_allele_col] = toupper(maf[,tumor_allele_col])
	}

	# collect tumor progression information
	if (! is.null(progression_col)) {
		if (is.factor(maf[,progression_col])) {
			message("Warning: You supplied tumor progression as a factor, but it will be converted to character.")
			message("The progression ordering will be what you supplied with \"progression_order\",")
			message("regardless of any ordering in your input MAF data frame.")
		}
		progressions = as.character(maf[,progression_col])
	} else {
		progressions = rep("1", nrow(maf)) # if analysis ignores tumor progression levels, all tumors get assigned level 1
		progression_order = c("1")
	}
	progressions = CESProgressions(samples = maf[,sample_col], sample_progressions = progressions, order = progression_order)
		

	# select only the necessary columns and give column names with consistency enforced by CESAnalysis class
	maf = maf[,c(sample_col, chr_col, start_col, ref_col, tumor_allele_col)]
	sample_col = "Unique_Patient_Identifier"
	chr_col = "Chromosome"
	start_col = "Start_Position"
	ref_col = "Reference_Allele"
	tumor_allele_col = "Tumor_Allele"
	colnames(maf) =  c(sample_col, chr_col, start_col, ref_col, tumor_allele_col)
  
  
	# strip chr prefixes from chr column, if present
	maf[,chr_col] = gsub('^chr', '', maf[,chr_col])
	

	# check for potential DNP/TNPs and separate them from data set
	# this is done before any other filtering to ensure catching as many of them as possible
	message("Searching for hidden multinucleotide variants...")
	num.prefilter = nrow(maf)
	dnp_tnp_results = cancereffectsizeR::DNP_TNP_remover(maf)
	maf = dnp_tnp_results$kept[,1:5] # To-do: fix return value of DNP_TNP
	pred.mnv.maf = dnp_tnp_results$removed
	num.pred.mnv = nrow(pred.mnv.maf)
	
	if (num.pred.mnv > 0) {
	  # To-do: move message to DNP_TNP_remover or otherwise ensure this description remains accurate
	  percent = round((num.pred.mnv / num.prefilter) * 100, 1)
	  message(paste0("Note: ", num.pred.mnv, " mutation records out of ", num.prefilter, " (", percent, "%) ",
	                 "are within 2 bp of other mutations in the same tumors."))
	  message("These records (saved as @excluded$pred.mnv.maf) will be excluded from effect size analysis.")
	}
	
	# Ensure reference alleles of mutations match reference genome (Note: Insertions won't match if their reference allele is "-")
	message("Checking that reference alleles match reference genome...")
	ref_chr = paste0('chr', maf[,chr_col])
	reference_alleles <- as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, ref_chr,
	                                                   strand="+", start=maf[,start_col], end=maf[,start_col]))
	num.prefilter = nrow(maf)
	reference.mismatch.maf = maf[maf[,ref_col] != reference_alleles,]
	num.nonmatching = nrow(reference.mismatch.maf)
	maf = maf[maf[,ref_col] == reference_alleles,]

	if (num.nonmatching > 0) {
	  percent = round((num.nonmatching / num.prefilter) * 100, 1)
	  message(paste0("Note: ", num.nonmatching, " mutation records out of ", num.prefilter, " (", percent, "%) ",
	                 "have reference alleles that do not actually match the reference genome."))
	  message("These records (saved as @excluded$reference.mismatch.maf) will be excluded from effect size analysis.")
	} else {
	  message("Reference alleles look good.")
	}

	# Count SNVs that will be included in effect size analysis
	message("Collecting all SNVs...")
	num.total = nrow(maf) # indels won't be filtered, but still another filtering steff later 
	bases = c("A","T","G","C")
	snv.maf = maf[maf[,ref_col] %in% bases  & maf[,tumor_allele_col] %in% bases,]
	num.non.snv = num.total - nrow(snv.maf)
	if (num.non.snv > 0) {
	  percent = round((num.non.snv / num.total) * 100, 1)
	  message(paste0("Note: ", num.non.snv, " mutation records out of ", num.total, " (", percent, "%) are not SNVs.\n",
	                 "(That is, ref or tumor allele were something other than a single A/C/T/G.)"))
	  message("Indels will be included in gene mutation rate calculation, but they will not be included in effect size analysis.")
	}
	
	num.good.snv = nrow(snv.maf)
	num.samples = length(unique(snv.maf[, sample_col]))
	message(paste0(num.good.snv, " SNVs from ", num.samples, " samples will be included in effect size analysis."))
  
  # declare environment to hold excluded MAFs
  excluded = new.env()
  excluded[["reference.mismatch.maf"]] = MAFdf(reference.mismatch.maf)
  excluded[["pred.mnv.maf"]] = MAFdf(pred.mnv.maf)


	# declare CESAnalysis object
	x = new("CESAnalysis", main.maf = MAFdf(maf), excluded = excluded, progressions = progressions)
	return(x)
}

