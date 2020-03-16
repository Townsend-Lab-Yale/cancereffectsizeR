# create list of all the trinuc numbers in every gene

# Load the pre-existing RefCDS object
hg19_dir = get_genome_data_directory("hg19")
RefCDS = readRDS(paste0(hg19_dir, "/RefCDS.rds"))

# Also using deconstructSigs_trinuc_string, part of CES sysdata
# Contains  92 context-specific SNV mutations, in the order used by deconstructSigs

# Convert stuff like "A[C>G]A" format to "ACA" (just the reference trinucleotide context)
context_names = sub("\\[([ACTG]).*\\]", "\\1", deconstructSigs_trinuc_string)

# deconstructSigs only includes C/T as central nucleotides; we need reverse complement for A/G in center
reverse_complement_names = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(context_names)))

gene_trinuc_comp  = new.env()

# go through all the transcripts in the RefCDS object...
for(i in 1:length(RefCDS)){
  if (i%%1000 == 0) {
    message(paste("Finished with", i, "RefCDS transcripts."))
  }

  total_counts = integer(96) # number of distinct deconstructSigs SNVs

  # for each transcripts, need to consider each exon and 1 base up/downstream for trinucleotide context
  intervals = RefCDS[[i]]$intervals_cds # two-column matrix with start/stop coordinates of each cds
  cds_lengths = intervals[,2] - intervals[,1] 

  # if transcript is on negative strand, flip exon order
  if (RefCDS[[i]]$strand == -1) {
    cds_lengths = rev(cds_lengths)
  }
  start = 1
  for (cds_length in cds_lengths) {
    end = start + cds_length
    # xscat and subseq are much more efficient than plain concatenation and subsetting
    seq = Biostrings::xscat(Biostrings::subseq(RefCDS[[i]]$seq_cds1up, start = start, width = 1), 
            Biostrings::subseq(RefCDS[[i]]$seq_cds, start = start, end = end), 
            Biostrings::subseq(RefCDS[[i]]$seq_cds1down, start = end, width = 1))
    
    # get trinucleotide counts for cds
    # this function returns full 64-cell table, including counts of 0 (otherwise subsetting below wouldn't work)
    cds_counts = Biostrings::trinucleotideFrequency(seq)

    # order the counts in deconstructSigs order and add them to total counts for this transcript
    total_counts = cds_counts[context_names] + cds_counts[reverse_complement_names] + total_counts

    # next cds sequence starts with next base in the seq_cds sequence
    start = end + 1
  }

  # add trinuc counts to environment (dropping names since they don't match deconstructSigs format)
  gene_trinuc_comp[[RefCDS[[i]]$gene_name]] = unname(total_counts)
}


output_file = paste0(hg19_dir, "/gene_trinuc_comp.rds")
saveRDS(gene_trinuc_comp,output_file)




