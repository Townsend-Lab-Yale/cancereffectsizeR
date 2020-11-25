
annovar_gnomad_exome = "~/progs/annovar/humandb/hg19_gnomad211_exome.txt"
annovar_gnomad_genome = "~/progs/annovar/humandb/hg19_gnomad211_genome.txt"

# genome file is large: to save space, throw out repeat-masked ranges
rmsk_file = "~/reference/rmsk.hg19.bed"

## Start with exome file
gnomad_exome = fread(annovar_gnomad_exome, nrows = 5)
cols_to_keep = which(colnames(gnomad_exome) %in% c("#Chr", "Start", "End", "non_cancer_AF_popmax")) # 1, 2, 3, 21

# 
# awk -F "\t" 'BEGIN {OFS = "\t"} $21 > .01 {print $1, $2, $3, $21}' hg19_gnomad211_exome.txt > common_germline_from_annovar_gnomad211_exome.txt
# 

common_exome = fread("~/progs/annovar/humandb/common_germline_from_annovar_gnomad211_exome.txt")

# verify all numeric, and already confirmed that they are MAF-like intervals (1-based, closed)
common_exome[, all(is.numeric(non_cancer_AF_popmax) & ! is.na(non_cancer_AF_popmax))]
common_exome = common_exome[, .(chr = `#Chr`, start = Start, end = End)]
exome_gr = makeGRangesFromDataFrame(common_exome, starts.in.df.are.0based = F)

## Repeat for genome file
gnomad_genome = fread(annovar_gnomad_genome, nrows = 5)

# Use AF_popmax because there is "non-cancer" cohort with this release's genome data
cols_to_keep = which(colnames(gnomad_genome) %in% c("#Chr", "Start", "End", "AF_popmax")) # 1, 2, 3, 7

# 
# awk -F "\t" 'BEGIN {OFS = "\t"} $7 > .01 {print $1, $2, $3, $7}' hg19_gnomad211_genome.txt > common_germline_from_annovar_gnomad211_genome.txt
# 
common_genome = fread("~/progs/annovar/humandb/common_germline_from_annovar_gnomad211_genome.txt")
common_genome[, all(is.numeric(AF_popmax) & ! is.na(AF_popmax))]
common_genome = common_genome[, .(chr = `#Chr`, start = Start, end = End)]
genome_gr = makeGRangesFromDataFrame(common_genome, starts.in.df.are.0based = F)

combined = sort(reduce(c(exome_gr, genome_gr)))

# Read in rmsk file and reduce ranges
rmsk = fread(rmsk_file)
rmsk = rmsk[, .(genoName, genoStart, genoEnd)]
rmsk = makeGRangesFromDataFrame(rmsk, seqnames.field = "genoName", start.field = "genoStart", end.field = "genoEnd",
                                keep.extra.columns = F, ignore.strand = T, starts.in.df.are.0based = T)
seqlevelsStyle(rmsk) = "NCBI" # remove chr prefixes to match curr_maf data
rmsk = rmsk[seqnames(rmsk) %in% c(1:22, 'X', 'Y')] # restrict to primary contigs

# This covers almost half the genome, almost entirely due to rmsk. However, there so many
# single-base ranges from common genome/exome that most of the ranges, and therefore most
# of the memory used, is for a small proportion of total coverage.
combined = sort(reduce(c(combined, rmsk)))


