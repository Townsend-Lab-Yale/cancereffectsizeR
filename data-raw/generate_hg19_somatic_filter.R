## 
library(data.table)
library(GenomicRanges)





## Step 1: Population variant frequency data source and generation
# Annovar has publicly available gnomAD data in an easy-to-use table format (large files)
# Need these columns from exome file: #Chr, Start, End, non_cancer_AF_popmax (columns 1, 2, 3, 21)
# awk -F "\t" 'BEGIN {OFS = "\t"} $21 > .01 {print $1, $2, $3, $21}' hg19_gnomad211_exome.txt > common_germline_from_annovar_gnomad211_exome.txt
# 
# For genome file, there is no "non-cancer" cohort, so use AF_popmax (field 7) instead
# awk -F "\t" 'BEGIN {OFS = "\t"} $21 > .01 {print $1, $2, $3, $21}' hg19_gnomad211_exome.txt > common_germline_from_annovar_gnomad211_exome.txt

# Read in data files
common_exome = fread("common_germline_from_annovar_gnomad211_exome.txt")

# verify all numeric, and already confirmed that they are MAF-like intervals (1-based, closed)
common_exome[, all(is.numeric(non_cancer_AF_popmax) & ! is.na(non_cancer_AF_popmax))]
common_exome = common_exome[, .(chr = `#Chr`, start = Start, end = End)]
exome_gr = makeGRangesFromDataFrame(common_exome, starts.in.df.are.0based = F)


common_genome = fread("common_germline_from_annovar_gnomad211_genome.txt")
common_genome[, all(is.numeric(AF_popmax) & ! is.na(AF_popmax))]
common_genome = common_genome[, .(chr = `#Chr`, start = Start, end = End)]
genome_gr = makeGRangesFromDataFrame(common_genome, starts.in.df.are.0based = F)
combined = sort(reduce(c(exome_gr, genome_gr)))

# .xz compression seems to be most efficient on GRanges (worth it for large GRanges)
saveRDS(combined, 'gnomad_common_variation_granges.rds', compress = 'xz')



## Step 2: RepeatMasker
# Download RepeatMasker hg19 annotations from UCSC genome browser
rmsk = fread("~/Documents/reference/rmsk.hg19.bed") # will respect BED coordinate format
rmsk = rmsk[, .(genoName, genoStart, genoEnd)]
rmsk = makeGRangesFromDataFrame(rmsk, seqnames.field = "genoName", start.field = "genoStart", end.field = "genoEnd",
                                keep.extra.columns = F, ignore.strand = T, starts.in.df.are.0based = T)
seqlevelsStyle(rmsk) = "NCBI" # remove chr prefixes to match curr_maf data
rmsk = rmsk[seqnames(rmsk) %in% c(1:22, 'X', 'Y')] # restrict to primary contigs
rmsk = reduce(rmsk)

saveRDS(rmsk, "repeat_masked_hg19_granges.rds", compress = 'xz')


