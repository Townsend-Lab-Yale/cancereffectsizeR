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
#' @import rtracklayer
#' @import GenomicRanges
#' @export


#-##' @examples
#-##' hg_converter(chain = "input_data/hg38Tohg19.chain", maf_to_convert = MAF_of_interest)






# library(rtracklayer)
# library(GenomicRanges)


hg_converter <- function(chain,maf_to_convert,chrom_col_name = "Chromosome",new_build_name = "Converted_from_GRCh38_to_hg19"){

  # flips the nucleotide to the other strand

  flip.function <- function(nucleotide){
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

  hg38.Grange <- makeGRangesFromDataFrame(df = maf_to_convert,
                                          keep.extra.columns = T,
                                          ignore.strand = F,
                                          seqinfo = NULL,
                                          seqnames.field = chrom_col_name,
                                          start.field = "Start_Position",
                                          end.field = "End_Position",
                                          strand.field = "Strand",
                                          starts.in.df.are.0based = F)
  GenomeInfoDb::genome(hg38.Grange) <- "GRCh38"

  GenomeInfoDb::seqlevelsStyle(hg38.Grange) = "UCSC"

  hg19.Grange <- rtracklayer::liftOver(hg38.Grange,ch)

  hg19.Grange <- unlist(hg19.Grange)
  GenomeInfoDb::genome(hg19.Grange) <- "hg19"

  message(paste("Number of rows in the MAF that failed to convert: ", length(hg38.Grange) - length(hg19.Grange)))

  hg19.df <- as.data.frame(hg19.Grange)

  colnames(hg19.df)[which(colnames(hg19.df)=="seqnames")] <- "Chromosome"
  colnames(hg19.df)[which(colnames(hg19.df)=="start")] <- "Start_Position"
  colnames(hg19.df)[which(colnames(hg19.df)=="end")] <- "End_Position"

  hg19.df$NCBI_Build <- new_build_name

  #reduce chromosome column to just numbers/letters
  hg19.df[,"Chromosome"] <- unlist(lapply(X = as.character(hg19.df[,"Chromosome"]),function(X) substr(X, start=4, stop=nchar(X))))

  #To keep things consistent with always being on the "+" strand, we need to flip the genes that were just flipped.
  to.flip <- which(hg19.df$strand=="-" & hg19.df$Variant_Type=="SNP")
  hg19.df$Reference_Allele[to.flip] <- unlist(lapply(X = hg19.df$Reference_Allele[to.flip],FUN = flip.function))
  hg19.df$Tumor_Seq_Allele1[to.flip] <- unlist(lapply(X = hg19.df$Tumor_Seq_Allele1[to.flip],FUN = flip.function))
  hg19.df$Tumor_Seq_Allele2[to.flip] <- unlist(lapply(X = hg19.df$Tumor_Seq_Allele2[to.flip],FUN = flip.function))

  return(hg19.df)
}








