#' Create a custom refset
#'
#' Use this function to create and save a directory of custom reference data that can be
#' used with cancereffectsizeR instead of supplied refsets like \code{ces.refset.hg19}. All
#' arguments are required except default_exome/exome_interval_padding, which are recommended.
#' 
#' To run this function, you'll need to have output from \code{build_RefCDS()}.
#' 
#' @param output_dir Name/path of an existing, writeable output directory where all data
#'   will be saved. The name of this directory will serve as the name of the custom refset.
#' @param refcds_output Transcript information in the two-item list (consisting of RefCDS
#'   and gr_genes) that is output by \code{build_RefCDS}.
#' @param species_name Name of the species, primarily for display (e.g., "human").
#' @param genome_build_name Name of the genome build, primarily for display (e.g., "hg19").
#' @param BSgenome_name The name of the BSgenome package to use (e.g., "hg19"); will used by
#'   cancereffectsizeR to load the reference genome via BSgenome::getBSgenome().
#' @param supported_chr Character vector of supported chromosomes. Note that
#'   cancereffectsizeR uses NCBI-style chromosome names, which means no chr prefixes ("X",
#'   not "chrX"). Mitochondrial contigs shouldn't be included since they would require
#'   special handling that hasn't been implemented.
#' @param default_exome_bed A BED file, typically acquired from some exome capture kit
#'   documentation, that defines default exome capture intervals for the genome. When users run
#'   cancereffectsizeR without specifying the genomic intervals that are covered by their
#'   sequencing data, these intervals inform some calculations. They can also be used to
#'   exclude non-covered variant calls. Including a default exome is optional but
#'   recommended, and users can always choose to user their own coverage definitions
#'   instead.
#' @param exome_interval_padding Number of bases to pad start/end of each covered
#'   interval, to allow for some variants to be called just outside of targeted regions,
#'   where there still may be pretty good sequencing coverage.
#' @export
create_refset = function(output_dir, refcds_output, species_name, genome_build_name, 
                         BSgenome_name, supported_chr = c(1:22, 'X', 'Y'), default_exome_bed = NULL,
                         exome_interval_padding = 100) {
  
  if (! is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir should be 1-length character (path for a new directory)")
  }
  
  if(! dir.exists(output_dir)) {
    stop("The supplied output_dir doesn't exist.")
  }
  
  if(file.access(output_dir, mode = 2) != 0) {
    stop("output_dir doesn't appear to be writeable.")
  }
  
  if (! is(refcds_output, "list") || length(refcds_output) != 2) {
    stop("refcds_output doesn't look valid; it should be a two-item list from build_RefCDS()")
  }
  
  if (! is.character(species_name) || length(species_name) != 1) {
    stop("species_name should be 1-length character")
  }
  if (! is.character(genome_build_name) || length(genome_build_name) != 1) {
    stop("genome_build_name should be 1-length character")
  }
  
  if (! is.character(BSgenome_name) || length(BSgenome_name) != 1) {
    stop("BSgenome_name should be 1-length character")
  }
  
  if (! is.character(supported_chr) || length(supported_chr) == 0) {
    stop("supported_chr should be type character")
  }
  
  if (! is.null(default_exome_bed)) {
    if(! is.character(default_exome_bed) || length(default_exome_bed) != 1) {
      stop("default_exome_bed should be 1-length character (file path).")
    }
    if (! file.exists(default_exome_bed)) {
      stop("File specified in default_exome_bed doesn't exist.")
    }
  } else {
    warning("No default exome supplied (okay, but recommended to include if available).")
  }
  
  if (! is.numeric(exome_interval_padding) || length(exome_interval_padding) != 1) {
    stop("exome_interval_padding should be 1-length numeric")
  }
  if (exome_interval_padding > 1000) {
    warning("You chose a high exome_interval_padding!")
  }
  
  # save an environment (dictionary) with some basic genome info
  genome_info = new.env(parent = emptyenv()) # emptyenv or when you call saveRDS, parents will get sucked in
  
  # for displaying to user
  genome_info[['build_name']] = genome_build_name
  genome_info[['species']] = species_name
  
  # name to use when loading genome with BSgenome::getBSgenome()
  genome_info[['BSgenome']] = BSgenome_name
  genome_info[['supported_chr']] = supported_chr
  
  exome_bed = default_exome_bed
  
  
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
  
  message("Loading reference genome and counting trinucleotide contexts...")
  wg = BSgenome::getSeq(bsg)
  genome_tri = Biostrings::trinucleotideFrequency(wg[supported_chr])
  genome_tri = colSums(genome_tri)
  
  genome_counts = genome_tri[context_names] + genome_tri[reverse_complement_names]
  genome_counts = data.frame(x = genome_counts)
  

  output_dir = sub('/+$', '', output_dir)
  if (! is.null(exome_bed)) {
    message("Loading default exome and counting trinucleotide contexts...")
    default_exome = rtracklayer::import(exome_bed, format = "bed")
    
    seqlevelsStyle(default_exome) = "NCBI"
    default_exome = default_exome[as.character(seqnames(default_exome)) %in% supported_chr]
    default_exome = unstrand(default_exome)
    # need to reorder seqlevels (1, 2, ..., X, Y) in cases where the input bed file wasn't sorted
    seqlevels(default_exome) = supported_chr
    default_exome = reduce(sort(default_exome), drop.empty.ranges = T)
    
    # Add interval padding of 100 nt
    start(default_exome) = start(default_exome) - exome_interval_padding
    end(default_exome) = end(default_exome) + exome_interval_padding
    default_exome = reduce(default_exome, drop.empty.ranges = T)
    
    seqinfo(default_exome) = seqinfo(bsg)
    default_exome = keepSeqlevels(default_exome, supported_chr) # drop extra seqlevels
    
    # get generic captured exome sequence and tabulate trinuc contexts
    exome_seq = getSeq(bsg, default_exome)
    exome_tri_contexts = Biostrings::trinucleotideFrequency(exome_seq)
    exome_tri_contexts = colSums(exome_tri_contexts)
    
    # reorder the counts as desired, then save as a data frame since that's what deconstructSigs wants
    exome_counts = exome_tri_contexts[context_names] + exome_tri_contexts[reverse_complement_names]
    exome_counts = data.frame(x = exome_counts)
    
    saveRDS(default_exome, paste0(output_dir, "/default_exome_gr.rds"))
    saveRDS(exome_counts, paste0(output_dir, "/tri.counts.exome.rds"))
  }
  
  message("Saving data files...")
  saveRDS(genome_info , paste0(output_dir, "/genome_build_info.rds"))
  saveRDS(genome_counts, paste0(output_dir, "/tri.counts.genome.rds"))
  
  refcds = refcds_output[[1]]
  gr_genes = refcds_output[[2]]
  saveRDS(refcds, paste0(output_dir, "/RefCDS.rds"))
  saveRDS(gr_genes, paste0(output_dir, "/gr_genes.rds"))
  saveRDS(unique(gr_genes$names), paste0(output_dir, "/gene_names.rds"))
  
}












# Optional: load and save default exome intervals ()
# If you don't set a default exome, the user must always supply coverage intervals.


