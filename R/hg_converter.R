#' hg_converter
#'
#' Convert MAF files from one genomic coordinate system to another.
#' Developed for use with MAF files found on the
#' NCI Genome Data Commons. Adapted from the liftOver vignettes page at \url{https://master.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html}
#'
#' @param chain A directory on your system with the "*over.chain" file to be used with liftOver. For instance, the chain to convert from hg38 to hg19 can to be downloaded from UCSC \url{http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/} (hg38ToHg19.over.chain.gz).
#' @param maf_to_convert MAF file, with relevant fields consisten with the GDC MAF specification (\url{https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/}), used in the conversion.
#' @param chrom_col_name Name of the column header containing chromosome information.
#' @param new_build_name This will be the string in the "NCBI_Build" column header of the new MAF file.
#' @return The returned output will be a MAF file with all the same data fields as the original but with "Start_Position" and "End_Position" converted to the build specified by the 'chain' input.
#' @export


#-##' @examples
#-##' hg_converter(chain = "input_data/hg38Tohg19.chain", maf_to_convert = MAF_of_interest)

hg_converter <- function(chain,maf_to_convert,chrom_col_name = "Chromosome",new_build_name = "Converted_from_GRCh38_to_hg19") {

  # flips the nucleotide to the other strand

  flip_function <- function(nucleotide){
    if(nucleotide=="A"){
      return("T")
    }
    else if(nucleotide=="T"){
      return("A")
    }
    else if(nucleotide=="C"){
      return("G")
    }
    else if(nucleotide=="G"){
      return("C")
    }
    else if(nucleotide=="N"){
      return("N")
    }
  }

  ch = rtracklayer::import.chain(chain)

  message("Loading in specified MAF...")

  # Keeping the convention of the original hg38 --> hg19 conversion
  # although this has now been updated for any conversion.

  hg38_Grange <- GenomicRanges::makeGRangesFromDataFrame(df = maf_to_convert,
                                          keep.extra.columns = T,
                                          ignore.strand = F,
                                          seqinfo = NULL,
                                          seqnames.field = chrom_col_name,
                                          start.field = "Start_Position",
                                          end.field = "End_Position",
                                          strand.field = "Strand",
                                          starts.in.df.are.0based = F)
  GenomeInfoDb::genome(hg38_Grange) <- "GRCh38"

  GenomeInfoDb::seqlevelsStyle(hg38_Grange) = "UCSC"

  hg19_Grange <- rtracklayer::liftOver(hg38_Grange,ch)

  hg19_Grange <- unlist(hg19_Grange)
  GenomeInfoDb::genome(hg19_Grange) <- "hg19"

  message(paste("Number of rows in the MAF that failed to convert: ", length(hg38_Grange) - length(hg19_Grange)))

  hg19_df <- as.data.frame(hg19_Grange)

  colnames(hg19_df)[which(colnames(hg19_df)=="seqnames")] <- "Chromosome"
  colnames(hg19_df)[which(colnames(hg19_df)=="start")] <- "Start_Position"
  colnames(hg19_df)[which(colnames(hg19_df)=="end")] <- "End_Position"

  hg19_df$NCBI_Build <- new_build_name

  #reduce chromosome column to just numbers/letters
  hg19_df[,"Chromosome"] <- gsub("chr", "", hg19_df[,"Chromosome"])

  #To keep things consistent with always being on the "+" strand, we need to flip the genes that were just flipped.
  to_flip <- which(hg19_df$strand=="-" & hg19_df$Variant_Type=="SNP")
  hg19_df$Reference_Allele[to_flip] <- unlist(lapply(X = hg19_df$Reference_Allele[to_flip],FUN = flip_function))
  hg19_df$Tumor_Seq_Allele1[to_flip] <- unlist(lapply(X = hg19_df$Tumor_Seq_Allele1[to_flip],FUN = flip_function))
  hg19_df$Tumor_Seq_Allele2[to_flip] <- unlist(lapply(X = hg19_df$Tumor_Seq_Allele2[to_flip],FUN = flip_function))

  return(hg19_df)
}

