library(GenomicRanges)

# name for the reference data set
ref_set_version = "ces_hg19_v1"

# path for BED file defining default exome coverage (optional, but see below)
exome_bed = 'xgen-exome-research-panel-targets.bed'

# where to save output reference data 
ref_set_dir = system.file(paste0("ref_sets/", ref_set_version), package = "cancereffectsizeR")


# save an environment (dictionary) with some basic genome info
genome_info = new.env()

# for displaying to user
genome_info[['build_name']] = "hg19"
genome_info[['species']] = "human"

# name to use when loading genome with BSgenome::getBSgenome()
genome_info[['BSgenome']] = "hg19"

# stick to primary contigs (mitochondrial variants would be fun to include, but there
# are several technical hurdles, unfortunately)
genome_info[['supported_chr']] = c(as.character(1:22), 'X', 'Y')
saveRDS(genome_info , paste0(ref_set_dir, "/genome_build_info.rds"))

# Compute and save genome-wide trinucleotide contexts
# (For ces_hg19_v1, this exactly reproduces the existing hg19 counts packaged with deconstructSigs.)
bsg = BSgenome::getBSgenome(genome_info[['BSgenome']])
seqlevelsStyle(bsg) = "NCBI"

# deconstructSigs_trinuc_string is internal to cancereffectsizeR (find with getAnywhere if necessary).
# Here, we need unique trinucleotide contexts (without mutations) in deconstructSigs ordering,
# which is why we produce this by converting from their column headings.
context_names = unique(sub("\\[([ACTG]).*\\]", "\\1", deconstructSigs_trinuc_string))
context_names = sort(context_names)
reverse_complement_names = unique(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(context_names))))

wg = BSgenome::getSeq(bsg)
genome_tri = Biostrings::trinucleotideFrequency(wg[good_chr])
genome_tri = colSums(genome_tri)

genome_counts = genome_tri[context_names] + genome_tri[reverse_complement_names]
genome_counts = data.frame(x = genome_counts)
saveRDS(genome_counts, paste0(ref_set_dir, "/tri.counts.genome.rds"))

# Optional: load and save default exome intervals ()
# If you don't set a default exome, the user must always supply coverage intervals.
default_exome = rtracklayer::import(exome_bed, format = "bed")
default_exome = reduce(sort(default_exome), drop.empty.ranges = T)
seqlevelsStyle(default_exome) = "NCBI"
good_chr = genome_info[['supported_chr']]
default_exome = default_exome[seqnames(default_exome) %in% good_chr]
default_exome = unstrand(default_exome)
# need to reorder seqlevels (1, 2, ..., X, Y) in cases where the input bed file wasn't sorted
seqlevels(default_exome)[seqlevels(default_exome) %in% good_chr] = good_chr

# Add interval padding of 100 nt
start(default_exome) = start(default_exome) - 100
end(default_exome) = end(default_exome) + 100
default_exome = reduce(default_exome, drop.empty.ranges = T)

seqinfo(default_exome) = seqinfo(bsg)
keepSeqlevels(default_exome, good_chr) # drop extra seqlevels
saveRDS(default_exome, paste0(ref_set_dir, "/default_exome_gr.rds"))

# get generic captured exome sequence and tabulate trinuc contexts
exome_seq = getSeq(bsg, default_exome)
exome_tri_contexts = Biostrings::trinucleotideFrequency(exome_seq)
exome_tri_contexts = colSums(exome_tri_contexts)

# reorder the counts as desired, then save as a data frame since that's what deconstructSigs wants
exome_counts = exome_tri_contexts[context_names] + exome_tri_contexts[reverse_complement_names]
exome_counts = data.frame(x = exome_counts)
saveRDS(exome_counts, paste0(ref_set_dir, "/tri.counts.exome.rds"))

