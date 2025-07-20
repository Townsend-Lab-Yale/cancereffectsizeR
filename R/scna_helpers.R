#' Determine even integer ploidy (2x, 4x, 8x)
#' 
#' Applies the approach used in Steele et al. to call sample ploidy.
#' We exclude chrX in samples not specifically annotated as female.
#'
#' @param cna_calls calls
#' @return A data.table of ploidy calls (2x, 4x, 8x).
calculate_ploidy = function(cna_calls) {
  calls = copy(cna_calls)
  calls[, size := end - start + 1]
  
  if(! 'sex' %in% names(calls)) {
    stop('Prepare the SCNA calls with prep_ASCAT3_segments() before calculating ploidy.')
  }
  
  if(anyNA(calls$sex)) {
    message('FYI, excluding chrX segments from ploidy calculations in samples with NA (or male) sex.')
  }
  
  # From Steele:
  # "Ploidy for each copy number profile was calculated as the relative length weighted sum of TCN
  # across a sample. The proportions of the genome that displayed LOH (pLOH) were also calculated.
  # Samples with a ploidy above −3/2 × pLOH + 3, meaning an LOH-adjusted ploidy of 3 or greater, were
  # deemed to be genome-doubled samples. By contrast, samples with a ploidy above −5/2 × pLOH + 5,
  # meaning an LOH-adjusted ploidy of 5 or greater, were deemed to be twice genome-doubled samples.
  # All other samples were considered as non-genome-doubled samples."
  ploidy_calc = calls[chr != 'X' | sex == 'F', .(our_prop_loh = .SD[nMinor == 0, sum(size)]/sum(size),
                                                 our_ploidy = .SD[, sum(total_copy * size)]/sum(size)),
                      by = 'sample']
  
  ploidy_calc[, our_wgd_call := as.numeric(our_ploidy > 3 - 1.5 * our_prop_loh)]
  ploidy_calc[, our_wgt_call := as.numeric(our_ploidy > 5 - 2.5 * our_prop_loh)]
  ploidy_calls = ploidy_calc[, .(sample, prop_LOH = our_prop_loh, mean_total_copy = our_ploidy, 
                                 is_wgd = our_wgd_call, is_wgt = our_wgt_call)]
  ploidy_calls[is_wgd == 0 & is_wgt == 0, simple_ploidy := 2]
  ploidy_calls[is_wgd == 1 & is_wgt == 0, simple_ploidy := 4]
  ploidy_calls[is_wgt == 1, simple_ploidy := 8]
  return(ploidy_calls)
}

#' Categorize copy changes relative to ploidy
#' 
#' @param cna_calls calls
#' @return A copy of the input table, with annotation columns added.
#' @keywords internal
annotate_copy_change_consequence = function(cna_calls) {
  # Including is.na(copy_change_detail) so that a logic error will reveal itself by leaving something unclassified
  calls = copy(cna_calls)
  calls[nMajor > exp_major_copy & nMinor > exp_minor_copy, 
        let(copy_change_detail = 'gain_gain', copy_change_class = 'gain')]
  calls[nMajor > exp_major_copy & nMinor == exp_minor_copy & is.na(copy_change_detail),
        let(copy_change_detail = 'gain_neutral', copy_change_class = 'gain')]
  
  # This case is possible for haploid chrX: nMinor has increased and nMajor hasn't.
  calls[nMajor == exp_major_copy & nMinor > exp_minor_copy & is.na(copy_change_detail),
        let(copy_change_detail = 'gain_neutral', copy_change_class = 'gain')]
  
  # Another haploid chrX case: nMajor has decreased and nMinor = exp_minor_copy.
  # (This can't happen otherwise, because nMinor is always less than nMajor.)
  # Within this case, when nMajor = 0 in diploid samples (that is, exp_minor_copy = 0), the assigned
  # category will be loss instead.
  calls[nMajor < exp_major_copy & nMinor == exp_minor_copy & nMajor > 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'neutral_decrease', copy_change_class = 'decrease')]
  
  # Last haploid chrX case: nMajor has decreased from expectation and nMinor has increased. No overall change.
  # Since our exp_major_copy is capped at 4 (for WGDx2 samples; in which case exp_minor_copy is 2), 
  # we can't have a case where major decreases, minor increases, and total copy changes (will always be 3/3).
  # We're calling this mixed rather than neutral.
  calls[nMajor < exp_major_copy & nMinor > exp_minor_copy & is.na(copy_change_detail),
        let(copy_change_detail = 'gain_decrease', copy_change_class = 'mixed')]
  
  
  calls[nMajor > exp_major_copy & nMinor < exp_minor_copy & nMinor > 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'gain_decrease', copy_change_class = 'mixed')]
  
  # Want exp_minor_copy > 0 due to haploid chrX.
  calls[nMajor > exp_major_copy & nMinor == 0 & exp_minor_copy > 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'gain_loss', copy_change_class = 'loh')]
  
  calls[nMajor == exp_major_copy & nMinor == exp_minor_copy & is.na(copy_change_detail),
        let(copy_change_detail = 'neutral_neutral', copy_change_class = 'neutral')]
  calls[nMajor == exp_major_copy & nMinor < exp_minor_copy & nMinor > 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'neutral_decrease', copy_change_class = 'decrease')]
  calls[nMajor == exp_major_copy & nMinor == 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'neutral_loss', copy_change_class = 'loh')]
  
  calls[nMajor < exp_major_copy & nMinor < exp_minor_copy & nMinor > 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'decrease_decrease', copy_change_class = 'decrease')]
  calls[nMajor < exp_major_copy & nMajor > 0 & nMinor == 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'decrease_loss', copy_change_class = 'loh')]
  
  calls[nMajor == 0 & is.na(copy_change_detail),
        let(copy_change_detail = 'loss_loss', copy_change_class = 'loss')]
  
  stopifnot(! anyNA(calls$copy_change_detail))
  return(calls[])
}

#' Make a BISCUT-style chromosome info data.frame
#' @param arm_coordinates data.table that defines the p/q arms of each chromosome.
#' @param acrocentromeric_chr character vector with the names of chromosomes that are acrocentromeric
#' @return description
make_biscut_chr_info = function(arm_coordinates, acrocentromeric_chr = c('13', '14', '15', '21', '22')) {
  
  if (! require("BISCUT")) {
    stop("BISCUT must be installed. Ask for help, because you need a custom version of it.")
  }
  coordinates = arm_coordinates
  acro_chr = acrocentromeric_chr
  required_cols = c('chr', 'p_start_effective', 'p_end', 'q_start', 'q_end_effective', 'cen_start', 'cen_end')
  missing_cols = setdiff(required_cols, names(coordinates))
  if(length(missing_cols) > 0) {
    stop('Required columns missing from coordinates table:\n',
         paste(missing_cols, collapse = ', '), '.')
  }
  # We substitute in p effect start and q effective end for p_start and q_end so that BISCUT can recognize
  # telomere-anchored calls. Omitting chrY.
  biscut_style_chr_info = coordinates[chr %in% c(1:22, 'X'), .(chromosome_info = chr, p_start = p_start_effective, 
                                                               p_end, q_start, q_end = q_end_effective, cen_start, cen_end)]
  
  # If q_end is not defined, it means that the chromosome (generally chrX) has no segments in the data set
  # in which q_end_effective coordinates were defined.
  biscut_style_chr_info = biscut_style_chr_info[! is.na(q_end)]
  
  # By BISCUT convention, p_start set to 1 for acrocentromeric chromosomes, but no events will be called on these arms.
  biscut_style_chr_info[chromosome_info %in% acro_chr, p_start := 1]
  
  # For BISCUT, we're going to re-define p_end and q_start to exclude centromeric regions
  biscut_style_chr_info[, p_end := cen_start - 1]
  biscut_style_chr_info[, q_start := cen_end + 1]
  
  # Only used for BISCUT visualizations; must supply values or BISCUT crashes (or used to)
  biscut_style_chr_info[, let(centromere = cen_start, size = q_end)] 
  stopifnot(! anyNA(biscut_style_chr_info))
  
  biscut_style_chr_info[chromosome_info == 'X', chromosome_info := '23']
  biscut_style_chr_info[, chromosome_info := as.numeric(chromosome_info)]
  biscut_style_chr_info = biscut_style_chr_info[order(chromosome_info)]
  biscut_style_chr_info = as.data.frame(biscut_style_chr_info) # validation function crashes if data.table
  biscut_style_chr_info = BISCUT::validate_chr_coordinates(biscut_style_chr_info)
  return(biscut_style_chr_info)
}

# BISCUT README:
# "If you have your data in absolute allelic copy number format, to obtain Segment_Mean for each
# segment, simply add the allelic copy numbers and then divide that by the ploidy of the sample and
# take log2."
# (README also says to set Num_Probes to 10 for non-probe data.)
#
# We're not going to do ploidy correction. It causes neutral segments to get called as deletions/amplifications!
prep_for_biscut_preprocess = function(calls) {
  calls = calls[, .(Sample = sample, Chromosome = chr, Start = start, End = end,
                    Num_Probes = 10, Segment_Mean = log2(total_copy/exp_total_copy))]
  
  
  #calls[, Segment_Mean := log2(total_copy/ploidy)]
  #calls[, c('total_copy', 'ploidy') := NULL]
  calls[Chromosome == 'X', Chromosome := '23']
  calls[, Chromosome := as.numeric(Chromosome)]
  return(calls[])
}

# Assumes we want telcent_thres = .001; assumes existence of coordinates table (created above)
read_in_biscut_preprocess = function(file_prefix) {
  bp_files = list.files(paste0(file_prefix, '/'), full.names = T)
  bp = rbindlist(lapply(bp_files, fread), idcol = 'file')
  setnames(bp, 'Sample', 'sample') # BISCUT uses Sample; this script uses sample
  bp[, file := basename(bp_files[file])]
  bp[, file := sub('\\.txt', '', file)]
  bp[, c('arm', 'ampdel', 'telcent') := tstrsplit(sub(paste0(basename(file_prefix), '_'), '', file), '_')]
  bp[, chr := sub('[pq]', '', arm)]
  bp[, pq := gsub('[^pq]', '', arm)]
  bp[nchar(pq) == 0, pq := 'q'] # deal with BISCUT acrocentromeric style
  
  # BISCUT recommends filtering out the very largest events as likely whole-arm. For now, we'll
  # also get rid of the very smallest events, as BISCUT suggests they could represent "noise."
  bp = bp[percent >= .001]
  which_whole_arm = bp[percent > .999, which = TRUE]
  
  bp_whole_arm = bp[which_whole_arm]
  bp = bp[! which_whole_arm]
  return(bp[])
}

find_consensus_changes = function(calls, region_col, region_frac_col, threshold = .99) {
  region_name = region_col
  calls = copy(calls)
  calls[, curr_region := calls[[region_name]]]
  calls[, curr_region_frac := calls[[region_frac_col]]]
  adjustments = unique(calls[, .(sample, curr_region)])
  
  for(allele_change in c('nMinor_change', 'nMajor_change')) {
    copy_change = calls[, .(prop = sum(curr_region_frac)), by = c('sample', 'curr_region', allele_change)]
    setnames(copy_change, allele_change, 'curr_change')
    # rank from greatest to least absolute copy number change in order to find the largest change that meets threshold
    amp_change = copy_change[curr_change > 0][order(sample, curr_region, -curr_change)]
    amp_change[, cum_prop := cumsum(prop), by = c('sample', 'curr_region')]
    del_change = copy_change[curr_change < 0][order(sample, curr_region, curr_change)]
    del_change[, cum_prop := cumsum(prop), by = c('sample', 'curr_region')]
    common_amp = amp_change[cum_prop > threshold, .SD[1], by = c('sample', 'curr_region')]
    common_del = del_change[cum_prop > threshold, .SD[1], by = c('sample', 'curr_region')]
    
    no_change = copy_change[curr_change == 0][, .(cum_prop = sum(prop)), by = c('sample', 'curr_region')]
    common_neutral = no_change[cum_prop > threshold]
    
    # Eventually need to deal with the issue that BISCUT peaks may have overlap here if threshold is <.5.
    stopifnot(merge.data.table(common_amp, common_del, by = c('sample', 'curr_region'), all = F)[, .N == 0])
    
    consensus_change_col = paste0('consensus_', allele_change)
    adjustments[common_amp, (consensus_change_col) := curr_change, on = c('sample', 'curr_region')]
    adjustments[common_del, (consensus_change_col) := curr_change, on = c('sample', 'curr_region')]
    adjustments[common_neutral, (consensus_change_col) := 0, on = c('sample', 'curr_region')]
  }
  
  if('accounted_nMajor_change' %in% names(calls)) {
    nonzero_upstream_nMajor = calls[accounted_nMajor_change != 0, .(accounted_nMajor_change), by = c('sample', 'curr_region')]
    nonzero_upstream_nMinor = calls[accounted_nMinor_change != 0, .(accounted_nMinor_change), by = c('sample', 'curr_region')]
    
    # We subtract the cumulative parent region changes from the current region changes. 
    adjustments[nonzero_upstream_nMajor, remaining_nMajor_change := consensus_nMajor_change - accounted_nMajor_change, on = c('sample', 'curr_region')]
    adjustments[nonzero_upstream_nMinor, remaining_nMinor_change := consensus_nMinor_change - accounted_nMinor_change, on = c('sample', 'curr_region')]
    
    adjustments[is.na(remaining_nMajor_change), remaining_nMajor_change := consensus_nMajor_change]
    adjustments[is.na(remaining_nMinor_change), remaining_nMinor_change := consensus_nMinor_change]
    
    # Neutral regions have consensus changes of zero (and no consensus changes in parent regions); NA changes are not consistent with neutrality
    adjustments[, is_neutral := complete.cases(.SD) & consensus_nMinor_change == 0 & consensus_nMajor_change == 0 & remaining_nMajor_change == 0 & remaining_nMinor_change == 0]
    
    
    adjustments[, let(consensus_nMajor_change = remaining_nMajor_change,
                      consensus_nMinor_change = remaining_nMinor_change)]
    adjustments[, c('remaining_nMajor_change', 'remaining_nMinor_change') := NULL]
  } else {
    calls$accounted_nMajor_change = 0
    calls$accounted_nMinor_change = 0
    adjustments[, is_neutral := complete.cases(.SD) & consensus_nMajor_change == 0 & consensus_nMinor_change == 0]
  }
  
  # Having recorded predominantly neutral regions, set NA changes to 0
  adjustments[is.na(adjustments)] = 0
  
  calls[adjustments, let(curr_nMinor_change = consensus_nMinor_change, 
                         curr_nMajor_change = consensus_nMajor_change), 
        on = c('sample', 'curr_region')]
  
  # Segments that can't absorb the consensus change (i.e., the change would go past neutrality) are exceptions
  calls[, is_exception := FALSE]
  
  # If applying the currently and previously accounted nMajor/nMinor changes takes an allele past neutrailty, the segment is an exception.
  # Examples: cumulative nMajor_change is 5, but nMajor <  5 + exp_major_copy;
  ###         cumulative nMinor_change is -2, but nMinor > exp_minor_copy - 2
  calls[(curr_nMinor_change > 0 & (nMinor < exp_minor_copy + curr_nMinor_change + accounted_nMinor_change)) |
          (curr_nMinor_change < 0 & (nMinor > exp_minor_copy + curr_nMinor_change + accounted_nMinor_change)) | 
          (curr_nMajor_change > 0 & (nMajor < exp_major_copy + curr_nMajor_change + accounted_nMajor_change)) |
          (curr_nMajor_change < 0 & (nMajor > exp_major_copy + curr_nMajor_change + accounted_nMajor_change)), is_exception := TRUE]
  
  calls[, region_is_neutral := FALSE]
  calls[adjustments[is_neutral == TRUE, .(sample, curr_region)], region_is_neutral := TRUE, on = c('sample', 'curr_region')]
  calls[region_is_neutral == TRUE & (nMajor != exp_major_copy | nMinor != exp_minor_copy), is_exception := TRUE]
  
  # Identify segments that are part of consensus regional changes
  calls[, region_has_call := region_is_neutral | curr_nMinor_change != 0 | curr_nMajor_change != 0]
  calls[is_exception == TRUE, let(curr_nMajor_change = 0, curr_nMinor_change = 0)]
  
  # Define region boundaries the same way that p_start_effective/q_end_effective were defined for chromosome coordinates
  # Non-neutral events get annotated with start/end denoting the endpoints of the affected part of the region.
  region_boundaries = calls[, .(region_start = min(start), region_end = max(end)), by = 'curr_region']
  region_events = unique(calls[region_has_call == TRUE & is_exception == FALSE,  
                               .(nMajor_change = curr_nMajor_change, nMinor_change = curr_nMinor_change, start = min(start), end = max(end)), 
                               by = c('sample', 'curr_region')])
  region_events = merge.data.table(region_events, region_boundaries, by = 'curr_region')
  region_events[nMajor_change == 0 & nMinor_change == 0, let(start = region_start, end = region_end)]
  region_events[, let(region_start = NULL, region_end = NULL)]
  
  
  region_change_detail = calls[region_has_call == TRUE & is_exception == FALSE, .(sum(curr_region_frac) > threshold), 
                               by = c('sample', 'curr_region', 'copy_change_detail')][V1 == TRUE, -"V1"]
  region_change_class = calls[region_has_call == TRUE & is_exception == FALSE, .(sum(curr_region_frac) > threshold), 
                              by = c('sample', 'curr_region', 'copy_change_class')][V1 == TRUE, -"V1"]
  
  region_events[region_change_detail, copy_change_detail := copy_change_detail, on = c('sample', 'curr_region')]
  region_events[region_change_class, copy_change_class := copy_change_class, on = c('sample', 'curr_region')]
  setnames(region_events, 'curr_region', region_name)
  
  # Current changes are 0 for exception segments and non-call regions
  calls[, let(accounted_nMajor_change = accounted_nMajor_change + curr_nMajor_change,
              accounted_nMinor_change = accounted_nMinor_change + curr_nMinor_change)]
  
  calls[, c('curr_region', 'curr_region_frac', 'region_is_neutral', 'curr_nMajor_change', 'curr_nMinor_change') := NULL]
  
  setnames(calls, 'is_exception', paste0('is_', region_name, '_exception'))
  setnames(calls, 'region_has_call', paste0('on_', region_name, '_call'))
  
  return(list(calls[], region_events[]))
}

#' @keywords internal
.prep_ascat_segments = function(calls, chr_coordinates, genome_build = 'hg38') {
  # Pretty wrapped user messages
  pmsg = function(...) {
    msg = paste0(strwrap(paste0(...)), collapse = "\n")
    message(msg)
  }
  if(! identical(genome_build, 'hg38')) {
    stop('Currently, genome_build must be "hg38".')
  }
  
  calls = copy(calls)
  if(! is.data.frame(calls)) {
    stop('calls must be a data.table.')
  } else {
    calls = as.data.table(calls)
  }
  required_cols = c('sample', 'chr', 'start', 'end', 'nMajor', 'nMinor')
  missing_cols = setdiff(required_cols, names(calls))
  if(length(missing_cols) > 0) {
    stop('Missing required columns in SNCA segment table:\n', paste0(missing_cols, collapse = ', '), '.')
  }
  
  calls[, chr := as.character(chr)]
  calls[, names(.SD) := lapply(.SD, as.numeric), .SDcols = c('start', 'end', 'nMajor', 'nMinor')]
  
  if(! 'sex' %in% names(calls)) {
    message('No \"sex\" column in calls. Therefore, chrX will be excluded from sample ploidy calculation.')
    calls$sex = NA_character_
  } else {
    calls$sex = as.character(calls$sex)
  }
  
  num_invalid_sex = calls[! is.na(sex) & ! sex %in% c('F', 'M'), .N]
  if(num_invalid_sex > 0) {
    stop("Invalid values in sex column. Acceptable values are F, M, <NA>.")
  }
  if(uniqueN(calls[, .(sample, sex)]) != uniqueN(calls$sample)) {
    stop('One or more samples has an inconsistent value in the sex column. A sample must have uniformly F, M or <NA>.')
  }
  calls[, total_copy := nMajor + nMinor]
  
  
  if(calls[, any(nMinor > nMajor)]) {
    stop('One or more segments have nMinor > nMajor')
  }
  
  if(calls[, any(end < start)]) {
    stop('There are segments with end < start!')
  }
  calls[, next_start := shift(start, -1), by = c('sample', 'chr')]
  num_overlapping_segments = calls[next_start == end, .N]
  
  # ASCAT output is 1-based (according to https://github.com/VanLoo-lab/ascat/issues/144#issuecomment-1540021594).
  if(num_overlapping_segments > 0) {
    pmsg('There are ', format(num_overlapping_segments, big.mark = ','), 
         " sample segments with an end coordinate matching the next segment's start coordinate.",
         " We will reduce the end coordinate by 1 for each of these.")
    calls[next_start == end & start != end, end := end - 1]
  }
  calls[, next_start := NULL]
  
  
  coordinates = copy(chr_coordinates)
  
  # If coordinates table omits chr prefixes, remove them from calls.
  if(! any(coordinates$chr %like% '^chr')) {
    calls[, chr := sub('^chr', '', chr)]
  }
  missing_chr = setdiff(calls$chr, coordinates$chr)
  if(length(missing_chr) > 0) {
    stop('Some chromosomes in calls are not present in chr_coordinates:\n', paste(missing_chr, collapse = ', '), '.')
  }
  
  # Get lowest/highest covered positions for each chromosome; verifying that they are consistent across samples.
  lowest_by_chr_and_sample = calls[, .(lowest = min(start)), by = c('chr', 'sample')]
  if (! uniqueN(lowest_by_chr_and_sample[, .(chr, lowest)]) == uniqueN(lowest_by_chr_and_sample$chr)) {
    stop("Expected the first segment of every chromosome to have the same start position for all samples.")
  }
  highest_by_chr_and_sample = calls[, .(highest = max(end)), by = c('chr', 'sample')]
  if (! uniqueN(highest_by_chr_and_sample[, .(chr, highest)]) == uniqueN(highest_by_chr_and_sample$chr)) {
    stop("Expected the last segment of every chromosome to have the same end position for all samples.")
  }
  
  coordinates[calls[, .(p_start_effective = min(start), 
                        q_end_effective = max(end)), by = 'chr'], 
              let(p_start_effective = p_start_effective,
                  q_end_effective = q_end_effective), on = 'chr']
  
  # In TCGA data, acrocentromeric chr13, 14, 15, 22 (but not acro chr 21) have no p segments in the
  # data set, so no arm calls will be made for these.
  coordinates[p_start_effective > p_end, p_start_effective := NA]
  
  ploidy_calls = calculate_ploidy(calls)
  
  wgd_samples = ploidy_calls[is_wgd == 1 & is_wgt == 0, sample]
  wgt_samples = ploidy_calls[is_wgt == 1, sample]
  diploid_samples = ploidy_calls[is_wgd == 0, sample]
  
  # Expected copy number is based on ploidy
  calls[wgt_samples, let(exp_total_copy = 8, exp_minor_copy = 4, exp_major_copy = 4), on = 'sample']
  calls[wgd_samples,let(exp_total_copy = 4, exp_minor_copy = 2, exp_major_copy = 2), on = 'sample']
  calls[is.na(exp_total_copy), let(exp_total_copy = 2, exp_minor_copy = 1, exp_major_copy = 1)]
  
  
  # Currently other genome builds aren't supported.
  if(genome_build == 'hg38' && 'X' %in% calls$chr) {
    # Confirmed via getSeq that these are correct hg38 PAR intervals:
    # chrX:10001-2781479 (PAR1), chrX:155701383-156030895 (PAR2)
    par1_start = 10001
    par1_end = 2781479
    par2_start = 155701383
    # Since PAR1 start is the start of resolved chrX (all Ns before), and PAR2 end is the end of chrX,
    # shouldn't need to worry about segments that are before PAR1 or after PAR2. We'll check PAR1 anyway.
    if(coordinates[chr == 'X', p_start_effective] < par1_start) {
      stop('Unexpected situation: p_start_effective for chrX precedes the start of PAR1.')
    }
    
    ## In males, split all segments that span PAR1 or PAR2 boundaries
    calls[, to_split_par1 := sex == 'M' & chr == 'X' & start <= par1_end & end > par1_end]
    par1_split_left = calls[to_split_par1 == TRUE]
    par1_split_right = copy(par1_split_left)
    par1_split_left[, end := par1_end]
    par1_split_right[, start := par1_end + 1]
    calls = rbind(calls[to_split_par1 == FALSE], par1_split_left, par1_split_right)
    
    calls[, to_split_par2 := sex == 'M' & chr == 'X' & start < par2_start & end >= par2_start]
    par2_split_left = calls[to_split_par2 == TRUE]
    par2_split_right = copy(par2_split_left)
    par2_split_left[, end := par2_start - 1]
    par2_split_right[, start := par2_start]
    calls = rbind(calls[to_split_par2 == FALSE], par2_split_left, par2_split_right)
    calls[, c('to_split_par1', 'to_split_par2') := NULL]
    
    # The non-PAR region of chrX has exp_minor_copy of 0 for males, regardless of sample ploidy.
    calls[sex == 'M' & chr == 'X' & start > par1_end & end < par2_start, exp_minor_copy := 0]
    calls[, exp_total_copy := exp_minor_copy + exp_major_copy]
  }
  
  # Annotate calls with gain/loss/LOH/etc. states
  calls = annotate_copy_change_consequence(calls)
  calls[, let(nMajor_change = nMajor - exp_major_copy, 
              nMinor_change = nMinor - exp_minor_copy)]
  ploidy_calls[calls, sex := sex, on = 'sample']
  
  calls[coordinates, starting_arm := ifelse(start < q_start, 'p', 'q'), on = 'chr']
  calls[coordinates, ending_arm := ifelse(end < q_start, 'p', 'q'), on = 'chr']
  
  # Calculate effective chr sizes
  chr_in_use = unique(calls$chr)
  coordinates[chr %in% chr_in_use, covered_chr_size :=  fcase(is.na(p_start_effective), q_end - q_start,
                                                      default = q_end_effective - q_start + p_end - p_start_effective)]
  
  covered_genome_size = sum(na.omit(coordinates$covered_chr_size))
  
  calls = assign_cosmic_CN_class(calls)
  calls[, c('CN', 'het', 'size_range') := tstrsplit(CNclass, ':')]
  
  # COSMIC uses consistent size bins except for homozygous deletions, which instead cap out at 1Mb. We
  # will recode to match the other events. The boundaries match those used in
  # assign_cosmic_CN_class(). Large homdels are rare.
  calls[het == 'homdel' & size_range == '>1Mb', size_range := fcase(end - start + 1 <= 1e7, '1Mb-10Mb',
                                                                 end - start + 1 <= 4e7, '10Mb-40Mb',
                                                                 default = '>40Mb')]
  
  prop_inc = calls[total_copy > exp_total_copy, .(prop_genome = sum(end - start)/covered_genome_size), 
                by = c('sample', 'size_range')]
  prop_inc[, prop_of_inc := prop_genome/sum(prop_genome), by = 'sample']
  
  prop_dec = calls[total_copy < exp_total_copy, .(prop_genome = sum(end - start)/covered_genome_size), 
                   by = c('sample', 'size_range')]
  prop_dec[, prop_of_dec := prop_genome/sum(prop_genome), by = 'sample']
  segment_proportions = list(increase = prop_inc, decrease = prop_dec)
  
  
  return(list(calls = calls, 
              ploidy_calls = ploidy_calls, 
              effective_coordinates = coordinates,
              segment_proportions = segment_proportions,
              genome_build = genome_build))
}

call_large_events = function(prepped_calls, arm_chr_threshold = .99,
                                ignore_centromeres_for_chr = FALSE,
                                account_biscut_regions = TRUE, biscut_dir = NULL,
                                run_biscut = FALSE, cores = parallel::detectCores()) {
  threshold = arm_chr_threshold
  
  # Fields for output events table, in a carefully chosen order
  final_event_cols = c("sample", "region_type", "region", "start", "end", "exp_total_copy", "nMajor_change", "nMinor_change", 
                       "cumulative_nMajor_change", "cumulative_nMinor_change", "consistent_impact", "consistent_impact_allele", 
                       "is_cumulative_neutral", "chr", "arm", "pq", "crosses_centromere", "ampdel", "telcent", "prop_other", 
                       "on_biscut_transcent", "biscut_percent", "cumulative_nMajor_change_p", "cumulative_nMinor_change_p", 
                       "cumulative_nMajor_change_q", "cumulative_nMinor_change_q")
  
  if(! is(prepped_calls, 'list') ||
     uniqueN(names(prepped_calls)) != length(prepped_calls) ||
     ! all(c('calls', 'effective_coordinates', 'genome_build') %in% names(prepped_calls))) {
    stop('prepped_calls should be a list as outputted by prep_ascat_segments().')
  }
  
  if(! rlang::is_scalar_integerish(cores) || cores < 1) {
    stop('cores should be 1-length positive integer.')
  }
  
  
  if(is.null(biscut_dir) && account_biscut_regions == TRUE) {
    stop('You must specify a BISCUT directory (for BISCUT-related output files) unless account_biscut_regions = FALSE.')
  }
  
  if(! rlang::is_bool(account_biscut_regions)) {
    stop('account_biscut_regions must be TRUE/FALSE.')
  }
  if(! rlang::is_bool(run_biscut)) {
    stop('run_biscut must be TRUE/FALSE.')
  }
  
  if(account_biscut_regions == FALSE && run_biscut == TRUE) {
    stop('Incompatible arguments: When run_biscut = TRUE, must have account_biscut_regions = TRUE.')
  }
  
  if(! is.null(biscut_dir) && account_biscut_regions == FALSE) {
    stop('A biscut_dir was specified even though account_biscut_regions was set FALSE')
  }
  
  if(! is.null(biscut_dir)) {
    if(! rlang::is_scalar_character(biscut_dir) || nchar(biscut_dir) < 1) {
      stop('biscut_dir should be type character (directory name/path)')
    }
    if(dir.exists(biscut_dir)) {
      stop('Specified directory for BISCUT files already exists.')
    }
  }
  
  if(! rlang::is_bool(ignore_centromeres_for_chr)) {
    stop('ignore_centromeres_for_chr should be TRUE/FALSE.')
  }
  
  calls = copy(prepped_calls$calls)
  coordinates = copy(prepped_calls$effective_coordinates)
  
  
  if(ignore_centromeres_for_chr) {
    coordinates[, effective_size := (cen_start - p_start_effective + 1) + (q_end_effective - cen_end + 1)]
  } else {
    # Effective size will just be the q size on some.
    coordinates[, effective_size := ifelse(is.na(p_start_effective),
                                           q_end_effective - q_start + 1,
                                           q_end_effective - p_start_effective + 1)]
  }
  
  calls = merge.data.table(calls, coordinates, by = 'chr')
  calls[, chr_frac := (end - start + 1) / effective_size]
  calls[, crosses_centromere := start < cen_start & end > cen_end]
  
  # For arm calls, we ignore centromeric regions (arguably should for chr calls, too...)
  calls[starting_arm == 'p', p_width := end - start + 1]
  calls[starting_arm == 'q', q_width := end - start + 1]
  
  # Calls that start or end within centromere will get negative p_width/q_width. 
  # We set p_frac/q_frac to zero for these segments.
  calls[starting_arm == 'p' & (ending_arm == 'q' | end > cen_start), 
        let(p_width = cen_start - start + 1,
            q_width = end - cen_end + 1)]
  
  calls[, p_size := cen_start - p_start_effective + 1]
  calls[, q_size := q_end_effective - cen_end + 1]
  calls[, p_frac := p_width / p_size]
  calls[, q_frac := q_width / q_size]
  
  calls[p_frac <= 0, p_frac := 0]
  calls[q_frac <= 0, q_frac := 0]
  
  
  chr_processed = find_consensus_changes(calls, region_col = 'chr', region_frac_col = 'chr_frac', 
                                         threshold = threshold)
  chr_accounted = chr_processed[[1]]
  chr_events = chr_processed[[2]]
  
  p_segments = chr_accounted[starting_arm == 'p' | ending_arm == 'p'][, arm := paste0(chr, 'p')]
  p_processed = find_consensus_changes(p_segments, region_col = 'arm', region_frac_col = 'p_frac', threshold = threshold)
  
  q_segments = chr_accounted[starting_arm == 'q' | ending_arm == 'q'][, arm := paste0(chr, 'q')]
  q_processed = find_consensus_changes(q_segments, region_col = 'arm', region_frac_col = 'q_frac', threshold = threshold)
  
  accounted_p = p_processed[[1]]
  accounted_q = q_processed[[1]]
  
  # Segments that cross centromere will appear in both accounted_p and accounted_q.
  num_unique_segments = uniqueN(rbind(accounted_p, accounted_q)[, .(sample, chr, start, end)])
  
  p_overlap = accounted_p[ending_arm == 'q']
  q_overlap = accounted_q[starting_arm == 'p']
  accounted_arm_no_overlap = rbind(accounted_p[starting_arm == ending_arm],
                                   accounted_q[starting_arm == ending_arm])
  
  # For segments that span p-arm and q-arm events, we use the nMinor/nMajor accounting that more fully explains nMinor change. (If the same,
  # we repeat with nMajor.) We prioritize nMinor so that LOH arms can be properly credited to the spanning segment.
  pq_overlap = rbindlist(list(p = p_overlap, q = q_overlap), idcol = 'pq')
  pq_overlap[, let(abs_nMinor_accounted = abs(accounted_nMinor_change), abs_nMajor_accounted = abs(accounted_nMajor_change))]
  
  # Have two events (p and q) to consider for each overlap
  stopifnot(pq_overlap[, .N, by = c('sample', 'chr', 'start', 'end')][, all(N == 2)])
  pq_resolved = pq_overlap[order(sample, chr, start, end, -on_arm_call, is_arm_exception, -abs_nMinor_accounted, -abs_nMajor_accounted), .SD[1], by = c('sample', 'chr', 'start', 'end')]
  
  accounted_arms = rbind(pq_resolved[, -c("pq", "abs_nMinor_accounted", "abs_nMajor_accounted")], accounted_arm_no_overlap)
  stopifnot(num_unique_segments == accounted_arms[, .N])
  
  # Rename accounted calls and get rid of chr/arm-related fields.
  accounted_calls = copy(accounted_arms)
  accounted_calls[, c('p_frac', 'q_frac', 'chr_frac', 'arm') := NULL]
  
  
  arm_events = rbindlist(list(p = p_processed[[2]], q = q_processed[[2]]), idcol = 'pq')
  arm_events[, chr := sub('[pq]$', '', arm)]
  
  # Neutral events get arm limits as boundaries; we're not calling them trans-centromeric
  arm_events = merge.data.table(arm_events, coordinates[, .(chr, p_start_effective, p_end, q_end_effective, cen_start, cen_end)], by = 'chr')
  arm_events[nMinor_change == 0 & nMajor_change == 0 & pq == 'p', 
             let(crosses_centromere = FALSE, start = p_start_effective, end = p_end)]
  arm_events[nMinor_change == 0 & nMajor_change == 0 & pq == 'q',
             let(crosses_centromere = FALSE, start = p_end + 1, end = q_end_effective)]
  arm_events[, crosses_centromere := (pq == 'p' & end > cen_end) | (pq == 'q' & start < cen_start)]
  arm_events[crosses_centromere == TRUE, prop_other := ifelse(pq == 'p', (end - p_end)/(q_end_effective - p_end),
                                                              (p_end - start + 1)/(p_end - p_start_effective + 1))]
  
  # The differences between these distributions may be worth some more thought.
  # arm_events[crosses_centromere == T & pq == 'p', MASS::truehist(prop_other)]
  # arm_events[crosses_centromere == T & pq == 'q', MASS::truehist(prop_other)]
  
  arm_events[, c('p_start_effective', 'p_end', 'q_end_effective', 'cen_start', 'cen_end') := NULL]
  
  chr_events[, let(region_type = 'chr', region = chr)]
  arm_events[, let(region_type = 'arm', region = arm)]
  
  # Propagate net changes
  arm_events[chr_events, let(cumulative_nMajor_change = i.nMajor_change + nMajor_change,
                             cumulative_nMinor_change = i.nMinor_change + nMinor_change), on = c('sample', 'chr')]
  arm_events[is.na(cumulative_nMajor_change), let(cumulative_nMajor_change = nMajor_change, cumulative_nMinor_change = nMinor_change)]
  chr_events[, let(cumulative_nMajor_change = nMajor_change, cumulative_nMinor_change = nMinor_change)]
  
  
  # Return early if not doing any BISCUT work.
  if(! account_biscut_regions) {
    chr_arm_events = rbind(chr_events, arm_events, fill = T)
    chr_arm_events[prepped_calls$ploidy_calls, let(exp_total_copy = simple_ploidy, sex = sex), on = 'sample']
    chr_arm_events[sex == 'M' & chr == 'X', exp_total_copy := exp_total_copy / 2]
    setnames(chr_arm_events, c('copy_change_class', 'copy_change_detail'), c('consistent_impact', 'consistent_impact_allele'))
    cols_to_use = intersect(final_event_cols, names(chr_arm_events))
    chr_arm_events = chr_arm_events[, .SD, .SDcols = cols_to_use]
    return(list(calls = calls, events = chr_arm_events))
  }
  
  
  biscut_style_chr_info = make_biscut_chr_info(arm_coordinates = coordinates)
  
  # Not currently necessary; would be if BISCUT accounting were run first (e.g., if spun off into its own function.)
  calls = accounted_calls
  if(! 'accounted_nMajor_change' %in% names(calls)) {
    calls[, let(accounted_nMajor_change = 0,
                accounted_nMinor_change = 0)]
  }
  
  # FYI, currently BISCUT preprocessing filters to autosomes.
  # We neutralize the segments with chr/arm changes so that BISCUT regions are not identified where arm/chr changes 
  # explain the would-be BISCUT regions.
  calls_for_biscut = calls[, .(sample, chr, start, end, nMajor = nMajor - accounted_nMajor_change,
                               nMinor = nMinor - accounted_nMinor_change, exp_total_copy)]
  calls_for_biscut[, let(total_copy = nMajor + nMinor, nMinor = NULL, nMajor = NULL)]
  for_biscut_input = prep_for_biscut_preprocess(calls_for_biscut)
  
  if(! dir.create(biscut_dir)) {
    stop('Unable to create BISCUT directory.')
  }
  biscut_input = paste0(biscut_dir, '/', 'BISCUT_preprocess_input.txt')
  fwrite(for_biscut_input, biscut_input, sep = "\t")
  biscut_breakpoint_dir = paste0(biscut_dir, '/', 'BISCUT_breakpoint_dir')
  message('Calling BISCUT telomeric and centromeric regions...')
  suppressMessages(BISCUT::make_breakpoint_files(segment_file = biscut_input, output_dir = biscut_breakpoint_dir,
                                                 chromosome_coordinates = biscut_style_chr_info, cores = cores))
  biscut_bp = read_in_biscut_preprocess(biscut_breakpoint_dir)
  biscut_bp = biscut_bp[, .SD, .SDcols = c('chr', setdiff(names(biscut_bp), names(coordinates)))]
  
  biscut_bp[, sample_region_id := paste(sample, ampdel, telcent, arm, sep = '.')]
  stopifnot(uniqueN(biscut_bp$sample_region_id) == biscut_bp[, .N])
  setnames(biscut_bp, c('startpos', 'endpos'), c('biscut_bp_start', 'biscut_bp_end'))
  
  # In our calls tables, start/end refer to ASCAT segments. Remove the name-conflicting columns
  # from BISCUT regions. (These BISCUT start/end are in fractional arm units and aren't needed.)
  biscut_bp[, c('start', 'end') := NULL] 
  
  # BISCUT omits q on acrocentromeric chromosomes
  biscut_bp[! arm %like% '[pq]$', arm := paste0(arm, 'q')]
  
  # Find and subtract consensus BISCUT changes within accounted_calls
  # ?foverlap: "Let [a,b] and [c,d] be intervals in x and y [...]. For type='within', the intervals
  # overlap iff a>=c and b<=d."
  setkey(biscut_bp, 'sample', 'chr', 'biscut_bp_start', 'biscut_bp_end')
  setkey(calls, 'sample', 'chr', 'start', 'end')
  
  calls_within_biscut = foverlaps(calls, biscut_bp, type = 'within', nomatch = NULL)
  calls_overlapping_biscut = foverlaps(calls, biscut_bp, type = 'any', nomatch = NULL)
  overlapping_not_within = calls_overlapping_biscut[!calls_within_biscut]
  
  # The copy segments that overlap BISCUT calls but are not entirely within the BISCUT calls are all
  # centromeric, with the BISCUT region cut off at cen_start or cen_end. These are all
  # valid centromere-anchored BISCUT regions. (BISCUT splits centromere-spanning segments in two,
  # often resulting in adjacent BISCUT regions.)
  stopifnot(overlapping_not_within[, all(biscut_bp_start == cen_end + 1 | biscut_bp_end == cen_start - 1)])
  
  # All copies segments that overlap BISCUT calls are eligible to have their BISCUT copy change assessed.
  # We'll get some duplicate segments, which will be handled.
  biscut_bp[, biscut_segment := paste(arm, ampdel, telcent, sep = '.')]
  overlapping_biscut = foverlaps(accounted_calls, biscut_bp[, -"arm"], type = 'any', nomatch = NULL)
  
  # Every BISCUT region is represented
  stopifnot(uniqueN(overlapping_biscut$sample_region_id) == uniqueN(biscut_bp$sample_region_id))
  not_overlapping_biscut = accounted_calls[! overlapping_biscut, on = c('sample', 'chr', 'start', 'end')]
  
  # We take width as max(end) - min(start). If some segments are absent in a sample, then it's possible that this summed
  # width is large enough to make it so that no copy changes meet biscut_adjust_threshold.
  overlapping_biscut[, biscut_width := max(end) - min(start) + 1, by = 'sample_region_id']
  overlapping_biscut[, biscut_segment_frac := (end - start + 1) / biscut_width]
  
  biscut_adjust_threshold = .5
  biscut_processed = find_consensus_changes(calls = overlapping_biscut, region_col = 'biscut_segment', 
                                            region_frac_col = 'biscut_segment_frac', threshold = biscut_adjust_threshold)
  
  # Remove BISCUT non-events.
  which_no_biscut_event = biscut_processed[[2]][nMajor_change == 0 & nMinor_change == 0, which = T]
  no_biscut_event = biscut_processed[[2]][which_no_biscut_event]
  biscut_processed[[2]] = biscut_processed[[2]][! which_no_biscut_event]
  biscut_accounted = biscut_processed[[1]]
  biscut_accounted[no_biscut_event, 
                   let(is_biscut_segment_exception = FALSE, on_biscut_segment_call = FALSE),
                   on = c('sample', 'biscut_segment')]
  
  # Some segments will be covered by more than one biscut region, due to overlapping tel/cent changes (sometimes on different alleles)
  # First, find which segments have more than one corresponding nMinor/nMajor change.
  
  # Note that the BISCUT tool avoids chrX
  #distinct_segment_changes = unique(biscut_events[, .(sample, chr, start, end, nMajor_change, nMinor_change)])
  distinct_segment_changes = unique(biscut_accounted[, .(sample, chr, start, end, accounted_nMajor_change, accounted_nMinor_change, 
                                                         on_biscut_segment_call, is_biscut_segment_exception, biscut_segment)])
  which_have_conflicts = distinct_segment_changes[, .N, by = c('sample', 'chr', 'start', 'end')][N > 1, -'N']
  conflict_segments = distinct_segment_changes[which_have_conflicts, on = names(which_have_conflicts)]
  conflict_segments[, let(abs_nMinor_accounted = abs(accounted_nMinor_change), 
                          abs_nMajor_accounted = abs(accounted_nMajor_change))]
  
  # As with chr/arm events, when multiple BISCUT changes overlap the same segment, we take the change that
  # brings the segment closest to expected minor copy; when tied, we use the major copy.
  biscut_resolved = conflict_segments[order(sample, chr, start, end, -on_biscut_segment_call, is_biscut_segment_exception, 
                                            -abs_nMinor_accounted, -abs_nMajor_accounted), .SD[1], by = c('sample', 'chr', 'start', 'end')]
  
  # Pull in the accepted nMinor/nMajor accounting
  biscut_accounted[biscut_resolved, let(accounted_nMajor_change = i.accounted_nMajor_change,
                                        accounted_nMinor_change = i.accounted_nMinor_change,
                                        on_biscut_segment_call = i.on_biscut_segment_call,
                                        is_biscut_segment_exception = i.is_biscut_segment_exception), 
                   on = c('sample', 'chr', 'start', 'end')]
  
  stopifnot(uniqueN(biscut_accounted[, .(sample, chr, start, end, 
                                         accounted_nMajor_change, accounted_nMinor_change, 
                                         is_biscut_segment_exception, on_biscut_segment_call)]) ==
              uniqueN(overlapping_biscut[, .(sample, chr, start, end)]))
  biscut_accounted = unique(biscut_accounted, by = c('sample', 'chr', 'start', 'end'))
  
  biscut_events = biscut_processed[[2]]
  biscut_events[, chr := sub('[pq]\\.(amp|del).*$', '', biscut_segment)]
  biscut_events[coordinates, let(crosses_centromere = start < cen_start & end > cen_end), on = 'chr']
  biscut_events[, let(ampdel = ifelse(biscut_segment %like% '\\.amp\\.', 'amp', 'del'),
                      pq = ifelse(biscut_segment %like% 'p\\.(amp|del)', 'p', 'q'),
                      telcent = ifelse(biscut_segment %like% '\\.tel$', 'tel', 'cent'))]
  biscut_events[, arm := paste0(chr, pq)]
  
  # Sometimes, there are centromere-bounded BISCUT calls on both p and q arms. For those that are
  # consistent with each other (as in, both amp or both del), we will note the possibility that they
  # are the same event.
  biscut_transcent = merge.data.table(biscut_events[pq == 'p' & crosses_centromere == T],
                                      biscut_events[pq == 'q' & crosses_centromere == T], 
                                      by = c('sample', 'chr'),
                                      all = F, suffixes = c('.p', '.q'))[end.p > start.q]
  
  # It may not be impossible to have inconsistent overlapping BISCUT calls, hasn't happened yet.
  stopifnot(biscut_transcent[, all(ampdel.p == ampdel.q)])
  
  # Usually the copy changes match (e.g., 98.9% in TCGA)
  # biscut_transcent[, mean(nMinor_change.p == nMinor_change.q & nMajor_change.p == nMajor_change.q)]
  
  # Take consensus changes
  biscut_transcent[, consensus_nMinor_change := ifelse(abs(nMinor_change.p) < abs(nMinor_change.q), nMinor_change.p, nMinor_change.q)]
  biscut_transcent[, consensus_nMajor_change := ifelse(abs(nMajor_change.p) < abs(nMajor_change.q), nMajor_change.p, nMajor_change.q)]
  
  transcent_events = biscut_transcent[, .(sample, chr, ampdel = ampdel.p, nMajor_change = consensus_nMajor_change, 
                                          nMinor_change = consensus_nMinor_change, crosses_centromere = TRUE, start = start.p, end = end.q)]
  biscut_events[, on_biscut_transcent := FALSE]
  biscut_events[transcent_events, on_biscut_transcent := TRUE, on = c('sample', 'chr', 'ampdel', 'crosses_centromere')]
  stopifnot(biscut_events[on_biscut_transcent == T, .N] / 2 == transcent_events[, .N])
  
  # biscut_events records BISCUT regions with consensus copy number changes.
  # Unlike arm/chr events, there is no BISCUT neutrality event category.
  stopifnot(biscut_events[nMinor_change == 0 & nMajor_change == 0, .N] == 0)
  
  biscut_accounted[biscut_events, on_biscut_transcent := on_biscut_transcent, on = c('sample', 'biscut_segment')]
  setnames(biscut_accounted, c('is_biscut_segment_exception', 'on_biscut_segment_call'),
           c('is_biscut_exception', 'on_biscut_call'))
  
  not_overlapping_biscut[, let(on_biscut_call = FALSE, is_biscut_exception = FALSE, on_biscut_transcent = FALSE)]
  biscut_accounted = biscut_accounted[, .SD, .SDcols = names(not_overlapping_biscut)]
  biscut_accounted = rbind(biscut_accounted, not_overlapping_biscut)
  
  # Error-check: No duplicate segments in biscut_accounted, and the exact sample segments from
  # accounted_calls are still present.
  stopifnot(uniqueN(biscut_accounted[, .(sample, chr, start, end)]) == biscut_accounted[, .N],
            biscut_accounted[, .N] == accounted_calls[, .N],
            accounted_calls[biscut_accounted, .N, on = c('chr', 'sample', 'start', 'end'), 
                            nomatch = NULL] == accounted_calls[, .N])
  
  biscut_accounted = biscut_accounted[, .(sample, chr, start, end, nMajor, nMinor, total_copy, exp_total_copy, 
                                          CNclass, het, size_range,
                                          starting_arm, ending_arm,
                                          crosses_centromere, accounted_nMajor_change, accounted_nMinor_change,
                                          cna_impact_allele = copy_change_detail, cna_impact = copy_change_class, 
                                          on_chr_call, on_arm_call, on_biscut_call, on_biscut_transcent, 
                                          is_chr_exception, is_arm_exception, is_biscut_exception)]
  
  biscut_events[biscut_bp, biscut_percent := percent, on = c('biscut_segment', 'sample')]
  biscut_events[, let(region_type = 'biscut', region = biscut_segment, biscut_segment = NULL)]
  transcent_events[, let(region_type = 'transcent', region = paste(chr, ampdel, 'transcent', sep = '.'))]
  
  biscut_events[arm_events, let(cumulative_nMajor_change = i.cumulative_nMajor_change + nMajor_change,
                                cumulative_nMinor_change = i.cumulative_nMinor_change + nMinor_change), on = c('sample', 'arm')]
  biscut_events[is.na(cumulative_nMajor_change), let(cumulative_nMajor_change = nMajor_change, cumulative_nMinor_change = nMinor_change)]
  
  transcent_events[arm_events[pq == 'p'], let(cumulative_nMajor_change_p = cumulative_nMajor_change + nMajor_change,
                                              cumulative_nMinor_change_p = cumulative_nMinor_change + nMinor_change), on = c('sample', 'chr')]
  transcent_events[arm_events[pq == 'q'], let(cumulative_nMajor_change_q = cumulative_nMajor_change + nMajor_change,
                                              cumulative_nMinor_change_q = cumulative_nMinor_change + nMinor_change), on = c('sample', 'chr')]
  transcent_events[is.na(cumulative_nMajor_change_p), cumulative_nMajor_change_p := nMajor_change]
  transcent_events[is.na(cumulative_nMajor_change_q), cumulative_nMajor_change_q := nMajor_change]
  transcent_events[is.na(cumulative_nMinor_change_p), cumulative_nMinor_change_p := nMinor_change]
  transcent_events[is.na(cumulative_nMinor_change_q), cumulative_nMinor_change_q := nMinor_change]
  
  # Since totally different arm events may cover the p/q sides of the transcent event, we only define cumulative changes
  # when they happen to be equal across arms (~80% of the time, here)
  # transcent_events[, mean(cumulative_nMajor_change_p == cumulative_nMajor_change_q &
  #                           cumulative_nMinor_change_p == cumulative_nMinor_change_q)]
  transcent_events[cumulative_nMajor_change_p == cumulative_nMajor_change_q, cumulative_nMajor_change := cumulative_nMajor_change_p]
  transcent_events[cumulative_nMinor_change_p == cumulative_nMinor_change_q, cumulative_nMinor_change := cumulative_nMinor_change_p]
  
  large_events = rbindlist(list(chr_events, arm_events, biscut_events, transcent_events), fill = TRUE)
  
  # Copy change annotations are potentially confusing in the events table. Better to get overall copy impact
  # from the annotated CNA calls.
  setnames(large_events, c('copy_change_class', 'copy_change_detail'), c('consistent_impact', 'consistent_impact_allele'))
  
  # biscut events are never neutral, and arms don't get called as neutral when their chr have non-neutral events.
  large_events[, is_cumulative_neutral := nMinor_change == 0 & nMajor_change == 0]
  large_events[prepped_calls$ploidy_calls, let(exp_total_copy = simple_ploidy, sex = sex), on = 'sample']
  large_events[sex == 'M' & chr == 'X', exp_total_copy := exp_total_copy / 2]
  large_events = large_events[, .SD, .SDcols = final_event_cols]
  
  fwrite(biscut_style_chr_info, paste0(biscut_dir, '/BISCUT_style_chr_info_', prepped_calls$genome_build, '.txt'), sep = "\t")
  output = list(calls = biscut_accounted, events = large_events, biscut_regions = biscut_bp)
  if(run_biscut) {
    biscut_out = do_biscut(breakpoint_file_dir = biscut_breakpoint_dir, 
                           results_dir = paste0(biscut_dir, '/biscut-out'), cores = cores, 
                           chromosome_coordinates = biscut_style_chr_info)
    output = c(output, list(biscut_output = biscut_out))
  }
  
  # Merge previous info into the new output
  output = c(output, prepped_calls[! names(prepped_calls) == 'calls'])
  return(output)
}

assign_cosmic_CN_class = function(cna_calls) {
  if(! 'total_copy' %in% names(cna_calls)) {
    cna_calls[, total_copy := nMinor + nMajor]
  }
  cna_calls[total_copy == 0, copy_group := '0']
  cna_calls[total_copy == 1, copy_group := '1']
  cna_calls[total_copy == 2, copy_group := '2']
  cna_calls[total_copy > 2 & total_copy < 5, copy_group := '3-4']
  cna_calls[total_copy > 4 & total_copy < 9, copy_group := '5-8']
  cna_calls[total_copy > 8, copy_group := '9+']
  cna_calls[nMajor == 0, het_status := 'homdel']
  cna_calls[nMinor == 0 & nMajor > 0, het_status := 'LOH']
  cna_calls[nMinor > 0, het_status := 'het']
  cna_calls[, width := end - start + 1]
  
  # FYI, it appears that Steele made the same choices for open/closed sides of intervals
  cna_calls[width <= 1e5, size_group := '0-100kb']
  cna_calls[width > 1e5 & width <= 1e6, size_group := '100kb-1Mb']
  cna_calls[width > 1e6 & width <= 1e7, size_group := '1Mb-10Mb']
  cna_calls[width > 1e7 & width <= 4e7, size_group := '10Mb-40Mb']
  cna_calls[width > 4e7, size_group := '>40Mb']
  
  # homdel group has fewer size bins (maxes out at 1Mb)
  cna_calls[nMajor == 0 & width > 1e6, size_group := '>1Mb']
  
  cna_calls[, CNclass := paste(copy_group, het_status, size_group, sep = ':')]
  cna_calls[, let(copy_group = NULL, width = NULL, size_group = NULL)]
  return(cna_calls[])
}


#' Tabulate copy changes by genomic window.
#'
#' Across a chromosome, counts the number of samples with gains, decreases, LOH, non-LOH decrease,
#' etc. All changes are with respect to expected ploidy.
#' @param chr Chromosome name as 1-length character
#' @param chr_ranges From get_genomic_windows()
#' @param account When TRUE (default), When account = TRUE, arm-level and chr-level changes are 
#' accounted for (i.e., ignored) when counting changes.
#' @export
get_chr_change_counts = function(calls, chr, chr_ranges, account = TRUE) {
  if(! 'exp_total_copy' %in% names(calls)) {
    stop('Need exp_total_copy in calls to identify diploid/nondiploid sample counts.')
  }
  ## Arguments no longer used (now doing raw burden, rather than per sample)
  # if(is.null(num_samples)) {
  #   stop('Must supply num_samples')
  # }
  # 
  # if(is.null(num_nondiploid)) {
  #   stop('Must supply num_nondiploid')
  # }
  curr_chr = chr
  calls = calls[chr == curr_chr]
  chr_ranges = chr_ranges[chr == curr_chr]
  
  # if(identical(account, TRUE)) {
  #   account = c('chr', 'arm', 'biscut')
  # }
  # if(! identical(account, FALSE)) {
  #   if(! is.character(account)) {
  #     stop('account must be TRUE/FALSE or a subset of c("biscut", "arm", "chr")')
  #   }
  #   account = tolower(unique(account))
  #   if(!(all(account %in% c('biscut', 'arm', 'chr')))) {
  #     stop('account must be TRUE/FALSE or a subset of c("biscut", "arm", "chr")')
  #   }
  #   if('biscut' %in% account) {
  #     
  #   }
  # }
  if(account == TRUE) {
    calls[on_biscut_call == F, nMajor := nMajor - accounted_nMajor_change]
    calls[on_biscut_call == F, nMinor := nMinor - accounted_nMinor_change] 
    
    # accounting may flip nMajor/nMinor
    calls[nMajor < nMinor, let(nMajor = nMinor, nMinor = nMajor)]
  }
  
  setkey(chr_ranges, starting, ending)
  setkey(calls, start, end)
  
  ol = foverlaps(chr_ranges, calls[, .(sample, chr, start, end, nMajor, nMinor, exp_total_copy)], nomatch = NULL)
  ol[, mp := floor((starting + ending) / 2) ]
  
  counts = ol[, .(loh = .SD[nMinor == 0, uniqueN(sample)],
                  decrease = .SD[nMajor + nMinor < exp_total_copy,
                                 uniqueN(sample)],
                  nonLOHdecrease = .SD[nMajor + nMinor < exp_total_copy & nMinor > 0, uniqueN(sample)],
                  # prop_decrease = .SD[nMajor + nMinor < exp_total_copy & nMinor > 0 & is_decrease_eligible == T,
                  #                     uniqueN(sample)]/num_nondiploid,
                  gain = .SD[nMajor + nMinor > exp_total_copy, uniqueN(sample)],
                  single_gain = .SD[nMajor + nMinor - exp_total_copy == 1, uniqueN(sample)],
                  technical_gain = .SD[nMajor + nMinor > 2, uniqueN(sample)],
                  covering = uniqueN(sample)), by = 'mp']
  counts$chr = curr_chr
  return(counts)
}

# should account be T/F?
get_local_rates = function(calls, chr_ranges, chr, account = FALSE, 
                               project_sample_counts = NULL, gene_coord) {
  if(! 'project' %in% names(calls)) {
    stop('Missing column in calls: project')
  }
  if(is.null(project_sample_counts)) {
    stop("You need to supply project_sample_counts (for now).")
  }
  change_counts = rbindlist(lapply(split(calls, calls$project), get_chr_change_counts, chr = chr, 
                                   chr_ranges = chr_ranges, account = account), idcol = 'project')
  
  # Use WGD samples for decrease
  change_counts[project_sample_counts, let(gain_rate = gain/N, decrease_rate = nonLOHdecrease/N_nondiploid), on = 'project'] #decrease_rate = decrease/N), on = 'project']
  rates = change_counts[, .(project, mp, gain_rate, decrease_rate, covering)]
  rates[, let(norm_gain = (gain_rate - min(gain_rate))/(max(gain_rate) - min(gain_rate)),
              norm_decrease = (decrease_rate - min(decrease_rate))/(max(decrease_rate) - min(decrease_rate))),
        by = 'project']
  
  
  gain_pam = rates[! is.nan(norm_gain), .(pamout = .(cluster::pam(norm_gain, k = 3)), med = median(norm_gain)), by = 'mp']
  gain_pam[, c('low', 'mid', 'high') := transpose(setDT(lapply(gain_pam$pamout, \(x) sort(x$medoids))))]
  
  decrease_pam =  rates[! is.nan(norm_decrease), .(pamout = .(cluster::pam(norm_decrease, k = 3)), med = median(norm_decrease)), by = 'mp']
  decrease_pam[, c('low', 'mid', 'high') := transpose(setDT(lapply(decrease_pam$pamout, \(x) sort(x$medoids))))]
  
  rates_by_mp = rbind(gain_pam[, .(mp, low, mid, high, med, type = 'gain')],
                    decrease_pam[, .(mp, low, mid, high, med, type = 'decrease')])
  
  max_changes = change_counts[, .(max_gain = max(gain_rate), max_decrease = max(decrease_rate)), by = c('project')]
  change_counts[max_changes, let(max_gain = max_gain, max_decrease = max_decrease), on = 'project']
  change_counts[rates_by_mp, let(r_gain = max_gain * med, r_decrease = max_decrease * med), , on = 'mp']
  
  
  project_by_sample = unique(calls[, .(sample, project)])
  rates = merge.data.table(project_by_sample, change_counts[, .(project, mp, r_gain, r_decrease)], 
                            by = 'project', all = T, allow.cartesian = TRUE)[, .(sample, mp, r_gain, r_decrease)]
  rates$chr = chr
  window_size = rates$mp[2] - rates$mp[1]
  rates$mp2 = rates$mp + window_size - 1
  setkey(rates, chr, mp, mp2)
  
  setkey(gene_coord, chr, start, end)
  ol = foverlaps(rates, gene_coord, nomatch = NULL)
  
  r_gain = ol[, .(r = max(r_gain), rate_type = 'increase'), by = c('gene', 'sample')] # or mean, etc.
  r_decrease = ol[, .(r = max(r_decrease), rate_type = 'decrease'), by = c('gene', 'sample')]
  
  final_r = rbind(r_gain, r_decrease) 
  final_r[, r0 := 1 - r]
  return(final_r[])
}

# Divides segment set into nonoverlapping size bins and gets change counts by bin.
get_counts_by_size = function(calls, chr, sizes, fn = get_chr_change_counts, chr_ranges = chr_ranges) {
  min_sizes = c(0, sizes)
  max_sizes = c(sizes, Inf)
  by_size = mapply(function(min_size, max_size) {
    fn(calls = calls[end - start < max_size & end - start > min_size], 
       chr = chr, chr_ranges = chr_ranges, account = F)
  }, min_sizes, max_sizes, SIMPLIFY = FALSE)
  
  min_size_format = fcase(min_sizes == 0, '',
                          min_sizes < 1e6, paste0(min_sizes/1e3, ' kb < '),
                          min_sizes >= 1e6, paste0(min_sizes/1e6, ' Mb < '))
  max_size_format = fcase(max_sizes < 1e6, paste0(' < ', max_sizes/1e3, ' kb'),
                          max_sizes >= 1e6 & ! is.infinite(max_sizes), paste0(' < ', max_sizes/1e6, ' Mb'),
                          is.infinite(max_sizes), '')
  size_labels = unname(unlist(mapply(function(x, y) {
    paste0(x, 'L', y)
  }, min_size_format, max_size_format, SIMPLIFY = FALSE)))
  
  names(by_size) = size_labels
  all_sizes = rbindlist(by_size, idcol = 'group')
  all_sizes[, group := factor(group, levels = size_labels)]
  return(all_sizes)
}


cna_class_relative_rates = function(cna_calls, cna_recon, ploidy_calls) {
  # Get total burden for segments with decrease
  # Leave out the homdel groups: they are more biased by selection (and are small part of attribution)
  # CHOICE: could alternatively leave them in
  
  cna_calls = copy(cna_calls)
  cna_recon = copy(cna_recon)
  
  # Annotate with ploidy calls (2/4/8)
  ploidy_calls = copy(ploidy_calls)[sample %in% cna_recon$sample]
  cna_recon[ploidy_calls, ploidy := simple_ploidy, on = 'sample']
  stopifnot(! anyNA(cna_recon$ploidy))
  
  diploid_samples = ploidy_calls[simple_ploidy == 2, sample]
  
  # This will be moved into CES internal data later.
  sr_info = data.table(pretty_label = c('L < 100 kb', '100 kb < L < 1 Mb', '1 Mb < L < 10 Mb', '10 Mb < L < 40 Mb', '40 Mb < L'),
                       cosmic_label = c('0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb', '>40Mb'),
                       bin_id = c('r1', 'r2', 'r3', 'r4', 'r5'),
                       bin_order = 1:5)
  
  
  for_burden_calc = cna_recon[size_range %in% sr_info$cosmic_label & het != 'homdel']
  
  diploid_decrease_burden = for_burden_calc[ploidy == 2 & CN %in% c('0', '1'), 
                                            .(b = sum(value)), by = c('sample', 'size_range')]
  diploid_increase_burden = for_burden_calc[ploidy == 2 & CN %in% c('3-4', '5-8', '9+'), 
                                            .(b = sum(value)), by = c('sample', 'size_range')]
  
  # quick solution for CN 3-4 group
  wgd_samples = ploidy_calls[simple_ploidy == 4, sample]
  prop3 = cna_calls[sample %in% wgd_samples & total_copy %in% c('3', '4'), mean(total_copy == 3)]
  for_burden_calc[, val2 := value]
  for_burden_calc[ploidy == 4 & CN == '3-4', val2 := prop3 * value]
  
  tetraploid_decrease_burden = for_burden_calc[CN %in% c('0', '1', '2', '3-4') & sample %in% wgd_samples, 
                                               .(b = sum(val2)), by = c('sample', 'size_range')]
  tetraploid_increase_burden = for_burden_calc[CN %in% c('5-8', '9+') & sample %in% wgd_samples,
                                               .(b = sum(value)), by = c('sample', 'size_range')]
  
  wgt_samples = ploidy_calls[simple_ploidy == 8, sample]
  prop_not8 = cna_calls[sample %in% wgt_samples & total_copy > 4 & total_copy < 9, mean(total_copy != 8)]
  for_burden_calc[ploidy == 8 & CN == '5-8', val2 := prop_not8 * value]
  octaploid_decrease_burden = for_burden_calc[ploidy == 8 & CN %in% c('0', '1', '2', '3-4', '5-8') & sample %in% wgt_samples,
                                              .(b = sum(val2)), by = c('sample', 'size_range')]
  octaploid_increase_burden = for_burden_calc[CN == '9+' & sample %in% wgt_samples,
                                              .(b = sum(value)), by = c('sample', 'size_range')]
  
  decrease_burden = rbindlist(list(diploid_decrease_burden, tetraploid_decrease_burden, octaploid_decrease_burden))
  decrease_burden[, rel_burden := b/sum(b), by = 'size_range']
  
  increase_burden = rbindlist(list(diploid_increase_burden, tetraploid_increase_burden, octaploid_increase_burden))
  increase_burden[, rel_burden := b/sum(b), by = 'size_range']
  
  # all samples have entries across all size ranges
  stopifnot(decrease_burden[, uniqueN(sample), by = 'size_range'][, all(V1 == uniqueN(decrease_burden$sample))],
            ! anyNA(decrease_burden))
  
  stopifnot(increase_burden[, uniqueN(sample), by = 'size_range'][, all(V1 == uniqueN(increase_burden$sample))],
            ! anyNA(increase_burden))
  
  return(list(increase_burden = increase_burden, decrease_burden = decrease_burden))
}

prep_gene_probability_model = function(gene_name, event_type = 'increase', 
                                       prepped_data = NULL,
                                       rates,
                                       gene_coord = ces.refset.hg38$cancer_gene_coord,
                                       project_by_sample, r1_info = NULL) {
  cna_calls = prepped_data$prepped_calls$calls
  cna_burdens = prepped_data$cna_burdens
  segment_prop = copy(prepped_data$prepped_calls$segment_proportions[[event_type]])
  setnames(segment_prop, c('prop_of_inc', 'prop_of_dec'), c('prop_of_inc_or_dec', 'prop_of_inc_or_dec'),
           skip_absent = TRUE)
  
  setkey(cna_calls, chr, start, end)
  stopifnot(rlang::is_scalar_character(event_type))
  event_type = tolower(event_type)
  stopifnot(event_type %in% c('decrease', 'increase'))
  
  # This will be moved into CES internal data later.
  sr_info = data.table(pretty_label = c('L < 100 kb', '100 kb < L < 1 Mb', '1 Mb < L < 10 Mb', '10 Mb < L < 40 Mb', '40 Mb < L'),
                       cosmic_label = c('0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb', '>40Mb'),
                       bin_id = c('r1', 'r2', 'r3', 'r4', 'r5'),
                       bin_order = 1:5)
  
  groups = sr_info$cosmic_label # size ranges
  gene_prob = rates[gene == gene_name]
  
  gr = gene_coord[gene == gene_name, .(chr, start, end)]
  setkey(gr, chr, start, end)
  curr_ol = foverlaps(cna_calls, gr, nomatch = NULL)
  
  curr_samples = unique(curr_ol$sample)
  curr_n = length(curr_samples)
  
  if(event_type == 'increase') {
    curr_rates = cna_burdens$increase_burden[sample %in% curr_samples]
  } else {
    curr_rates = cna_burdens$decrease_burden[sample %in% curr_samples]
  }
  
  
  
  
  if(gene_prob[, .N] == 0) {
    warning('gene not in the chromosomes used; returning empty table')
    return(NULL)
  }
  

  gene_prob[segment_prop, prop_genome := sum(prop_genome), on = 'sample']
  gene_prob[is.na(prop_genome), prop_genome := 0 ] # hmm
  gene_prob[project_by_sample, project := project, on = 'sample']
  gene_prob[, r2 := .N * r * prop_genome/sum(prop_genome), by = 'project']
  # 
  # > gene_prob[, max(r2)]
  # [1] 0.1686747
  
  curr_rates[segment_prop, curr_prop := prop_of_inc_or_dec, on = c('sample', 'size_range')]
  curr_rates[is.na(curr_prop), curr_prop := 0]
  curr_rates[gene_prob, r2 := r2, on = 'sample']
  curr_rates[, r := r2 * curr_prop]
  
  
  curr_rates[sr_info, bin_id := bin_id, on = c(size_range = 'cosmic_label')]
  
  # If more than one event in a sample, use count the largest event
  if(event_type == 'decrease') {
    curr_ol = curr_ol[total_copy < exp_total_copy]
  } else if(event_type == 'loh') {
    curr_ol = curr_ol[nMinor == 0] # Yeah, total_copy might increase! Yet something was lost.
  } else if(event_type == 'increase') {
    curr_ol = curr_ol[total_copy > exp_total_copy]
  }
  
  if(! is.null(r1_info)) {
    
    curr_interval = gene_coord[gene == curr_gene_name, .(chr, int_start = start, int_end = end)]
    setkey(curr_interval, chr, int_start, int_end)
    other_interval = r1_info$interval
    si_other = r1_info$si
    curr_rel_neut = copy(calls)
    
    if(curr_interval$int_start > other_interval$end) {
      dist_to_other = curr_interval$int_start - other_interval$end
      rel_pos = 'right'
    } else if(curr_interval$int_end < other_interval$start) {
      dist_to_other = other_interval$start - curr_interval$int_end
      rel_pos = 'left'
    } else {
      stop('Gene intervals overlap.')
    }
    
    p_by_size = list()
    # order doesn't matter here
    for(curr_bin_id in sr_info$bin_id) {
      curr_size_range = sr_info[bin_id == curr_bin_id, cosmic_label]
      curr_rel_neut = calls[size_range == curr_size_range]
      curr_ol = foverlaps(curr_rel_neut, curr_interval, nomatch = NULL)
      ol_both = foverlaps(curr_ol, other_interval, nomatch = NULL)
      if(ol_both[, .N] == 0) {
        p_by_size[[curr_bin_id]] = 0
      } else {
        p_by_size[[curr_bin_id]] = ol_both[, .N]/curr_ol[, .N]
      }
    }
    
    p_by_size = setNames(as.numeric(p_by_size), names(p_by_size))
    
    curr_rates[, curr_p := p_by_size[bin_id]]
    curr_rates[, r := r * (1 - curr_p + si_other * curr_p)]
    
    min_by_type = curr_rates[r != 0, .(min(r)), by = 'bin_id']
    curr_rates[r == 0, r := min_by_type[bin_id, V1, on = 'bin_id']]
    rm(curr_rel_neut)
  }
  
  
  total_rate = curr_rates[, .(r_sum = sum(r)), by = 'sample']
  total_rate[, r0 := 1 - r_sum]
  with_none = setdiff(curr_samples, curr_ol$sample)
  
  samples_with = list()
  
  # We go from largest to smallest bin
  for(i in sort(sr_info$bin_order, decreasing = TRUE)) {
    curr_bin = sr_info[bin_order == i, bin_id]
    curr_size_range = sr_info[bin_order == i, cosmic_label]
    
    curr_samples_with = curr_ol[size_range == curr_size_range, unique(sample)]
    
    # Larger events in the same samples take precedence (already handled above, really)
    samples_with[[curr_bin]] = setdiff(curr_samples_with, unlist(samples_with))
  }
  rm(cna_calls)
  rm(segment_prop)
  rm(gene_prob)
  rm(prepped_data)
  rm(r1_info)
  lik = function(gamma) {
    sum_log_lik = 0
    if (length(with_none) > 0) {
      r0 = total_rate[with_none, r0, on = 'sample']
      prob = 1 - ((gamma * (1 - r0))/ (r0 + gamma * (1 - r0)))
      sum_log_lik = sum(log(prob))
    }
    
    for(bin_id in sr_info$bin_id) {
      curr_samples = samples_with[[bin_id]]
      if(length(curr_samples) > 0) {
        r0 = total_rate[curr_samples, r0, on = 'sample']
        curr_bin_id = bin_id
        r = curr_rates[bin_id == curr_bin_id][curr_samples, r, on = 'sample']
        prob = (gamma * r)/(r0 + gamma * (1 - r0))
        sum_log_lik = sum_log_lik + sum(log(prob))
      }
    }
    return(-1 * sum_log_lik)
  }
  
  formals(lik)[["gamma"]] = 1
  bbmle::parnames(lik) = "gamma"
  gene_interval = gr[, .(chr = chr[1], start = min(start), end = max(end))]
  output = list(lik = lik, 
                groups = samples_with,
                num_covering = length(curr_samples),
                interval = gene_interval,
                gene_name = gene_name)
  return(output)
}

# sig_def = system.file('extdata/COSMIC_v3.4_CN_GRCh37.txt', package ='ces.refset.hg38')
cn_signature_extraction = function(sig_def, cna_segments) {
  cn_sig = fread(sig_def)
  cn_sig = transpose(cn_sig, make.names = 'Type', keep.names = 'signature')
  cn_names = cn_sig$signature
  cn_sig = as.matrix(cn_sig[, -'signature'])
  rownames(cn_sig) = cn_names
  cn_sig = t(cn_sig)
  cn_mat = t(table(cna_segments$sample, 
                   factor(cna_segments$CNclass, levels = rownames(cn_sig))))
  
  # To-do: make multicore
  mp_out = MutationalPatterns::fit_to_signatures_strict(cn_mat, cn_sig)
  
  # Make a melted version
  recon = melt(as.data.table(mp_out$fit_res$reconstructed, keep.rownames = 'CNclass'),
               id.vars = 'CNclass', variable.name = 'sample', variable.factor = FALSE)
  recon[, c('CN', 'het', 'size_range') := tstrsplit(CNclass, ':')]
  recon[, prop_by_class := value/sum(value), by = c('sample', 'CN', 'het')] ## Needed?
  
  return(list(reconstruction = recon, mp_out = mp_out))
}

get_seg_increase_rates = function(cna_calls, increase_burden) {
  # For each size bin, get the median size event for each sample.
  typical_increase_sizes = cna_calls[nMajor + nMinor > exp_total_copy, .(typical_size = floor(median(end - start))), by = c('sample', 'size_range')]
  
  avg_typical_increase = typical_increase_sizes[, .(typical_size = floor(mean(typical_size))), by = 'size_range']
  
  # Rate over a 100kb window that events of a given size range hit that sample.
  increase_burden[typical_increase_sizes, typical_size := typical_size, on = c('sample', 'size_range')]
  
  # Some NAs in typical_size because some samples have no gains over ploidy for some size ranges. Plug in the averages.
  increase_burden[is.na(typical_size), typical_size := avg_typical_increase[size_range, typical_size, on = 'size_range']]
  
  # Note these are all nearly the same, ~2.9 Gb
  covered_size_by_sample = cna_calls[, .(cov_size = sum(end - start + 1)), by = 'sample']
  increase_burden[covered_size_by_sample, cov_size := cov_size, on = 'sample']
  
  # Ignoring window size because 1e5 << 2.9e9. Is that okay?
  # Note some unrealistic values for the >40Mb size range.
  increase_burden[, seg_rate := b * (typical_size/cov_size)]
  
  
  seg_rates = increase_burden[, .(sample, size_range, seg_rate)]
  total_rates = seg_rates[, .(seg_rate = sum(seg_rate)), by = 'sample']
  
  return(list(seg_rates = seg_rates, total_rates = total_rates))
}


get_disjoint_gene_coord = function(gene_coord) {
  gene_coord = copy(gene_coord)[, .(gene, cancer_anno, chr, 
                                    start = fcase(is.na(mane_start), start, 
                                                  default = mane_start),
                                    end = fcase(is.na(mane_end), end,
                                                default = mane_end))]
  gene_coord[, width := end - start + 1]
  setkey(gene_coord, chr, start, end)
  gene_ol = foverlaps(gene_coord, gene_coord, nomatch = NULL)[gene != i.gene]
  
  # Every overlap gets listed twice, once with gene/i.gene = gene A, gene B; once with gene/i.gene =
  # gene B, gene A.
  not_ol = gene_coord[! gene %in% gene_ol$gene]
  
  # Sort by cancer status then decreasing size (so that we favor cancer genes and longer genes) 
  gene_ol = gene_ol[(cancer_anno != 'noncancer' & i.cancer_anno == 'noncancer') | width > i.width]
  to_remove = gene_ol$i.gene
  deoverlapped = gene_ol[! gene %in% to_remove, .(gene, cancer_anno, chr, start, end)]

  gene_coord = rbind(deoverlapped, not_ol[, -"width"])
  gene_coord = unique(gene_coord[order(chr, start)])
  return(gene_coord[])
}

run_seg_multi = function(selection_loci, cna_calls = NULL, event_type = 'increase',
                         seg_rates = NULL, debug = FALSE) {
  if(is.null(seg_rates)) {
    stop('Need seg_rates (from get_seg_increase_rates(), for example).')
  }
  if(is.null(cna_calls)) {
    stop('Need cna_calls.')
  }
  if(! identical(event_type, 'increase')) {
    stop('Currently, must have event_type = "increase".')
  }
  if(! is.data.table(selection_loci)) {
    stop('selection_loci should be data.table.')
  }
  if(! all(c('gene', 'chr', 'start', 'end') %in% names(selection_loci))) {
    stop('selection_loci must have fields gene, chr, start, end.')
  }
  
  total_rates = copy(seg_rates$total_rates)
  bin_rates = copy(seg_rates$seg_rates)
  setkey(cna_calls, chr, start, end)
  
  # This will be moved into CES internal data later.
  sr_info = data.table(pretty_label = c('L < 100 kb', '100 kb < L < 1 Mb', '1 Mb < L < 10 Mb', '10 Mb < L < 40 Mb', '40 Mb < L'),
                       cosmic_label = c('0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb', '>40Mb'),
                       bin_id = c('r1', 'r2', 'r3', 'r4', 'r5'),
                       bin_order = 1:5)
  
  groups = sr_info$cosmic_label # size ranges
  
  stopifnot(all(c('sample', 'size_range', 'seg_rate') %in% names(bin_rates)))
  curr_rates = bin_rates[size_range %in% groups]
  
  if(curr_rates[, .N] == 0) {
    warning('No rates found; returning empty table')
    return(NULL)
  }
  
  # To-do: Insert check that loci are non-overlapping.
  gr = selection_loci
  stopifnot(uniqueN(gr$gene) == gr[, .N])
  setkey(gr, chr, start, end)
  curr_ol = foverlaps(cna_calls, gr, nomatch = NULL)
  
  num_genes = gr[, .N]
  incomplete_cov = curr_ol[, .(N = uniqueN(gene)), by = 'sample'][N < num_genes, sample]
  
  # To-do: implement better handling of incomplete coverage
  if(length(incomplete_cov) > 0) {
    curr_ol = curr_ol[! sample %in% incomplete_cov]
  }
  if(curr_ol[, .N] == 0) {
    stop('There are no samples with segments covering all loci!')
  }
  num_samples_used = uniqueN(curr_ol$sample)
  
  
  casted = dcast(curr_ol[, .(chr, i.start, i.end, size_range, gene, sample, is_inc = nMajor + nMinor > exp_total_copy)], 
                 ... ~ gene, value.var = 'gene', fun.aggregate = length)
  gene_names = gr$gene
  
  # length() in dcast made columns integer; need numeric.
  casted[, names(.SD) := lapply(.SD, as.numeric), .SDcols = gene_names]
  
  
  # What do to when a sample both has and doesn't have an increase/decrease? That is, there is at
  # least one increase/decrease segment overlapping the gr, and some other segment(s).
  # For now, we'll say that yes, we have an event. And we'll use just the largest segment of that class.
  
  # To-do: We could choose the smallest event, due to its rarity
  melted = melt(casted, measure.vars = gene_names, variable.name = 'gene')[value == 1]
  
  melted[sr_info, bin_order := bin_order, on = c(size_range = 'cosmic_label')]
  
  yes_seg_to_use = melted[is_inc == T, .SD[which.max(bin_order)], 
                          by = c('sample', 'gene')][, .(sample, chr, i.start, i.end, gene)]
  no_seg_to_use = melted[is_inc == F, .SD[which.max(bin_order)], 
                         by = c('sample', 'gene')][, .(sample, chr, i.start, i.end, gene)]
  
  # Might have a situation where yes/no segments both cover the gene
  no_seg_to_use = no_seg_to_use[! yes_seg_to_use, on = c('sample', 'gene')]
  
  ## Use this in unit testing
  #stopifnot((no_seg_to_use[, .N] + yes_seg_to_use[, .N]) / 4 == num_samples_used)
  
  yes_seg_to_use[, gene := NULL]
  yes_seg_to_use = unique(yes_seg_to_use)
  no_seg_to_use[, gene := NULL]
  no_seg_to_use = unique(no_seg_to_use)
  
  no_event = casted[no_seg_to_use, on = names(no_seg_to_use)]
  yes_event = casted[yes_seg_to_use, on = names(yes_seg_to_use)]

  
  yes_event[bin_rates, neutral_rate := seg_rate, on = c('sample', 'size_range')]
  no_event[total_rates, neutral_rate := seg_rate, on = 'sample']
  
  original_yes_event = copy(yes_event)
  original_no_event = copy(no_event)
  yes_event = yes_event[, .SD, .SDcols = c('sample', 'neutral_rate', gene_names)]
  no_event = no_event[, .SD, .SDcols = c('sample', 'neutral_rate', gene_names)]
  
  lik = function(si) {
    
    if(length(si) > 1) {
      yes_event[, r := 1]
      no_event[, r := 1]
      for(locus in names(si)) {
        setnames(yes_event, locus, 'curr_locus')
        setnames(no_event, locus, 'curr_locus')
        yes_event[, curr_r := fcase(curr_locus == 0, 1,
                                    default = si[locus])]
        no_event[, curr_r := fcase(curr_locus == 0, 1,
                                   default = si[locus])]
        setnames(yes_event, 'curr_locus', locus)
        setnames(no_event, 'curr_locus', locus)
        yes_event[, r := r * curr_r]
        no_event[, r := r * curr_r]
      }
    } else {
      yes_event[, r := si]
      no_event[, r := si]
    }

    
    # Note: log(-1 * expm1(-1 * r * neutral_rate))] is equivalent to 
    # log(1 - exp(-1 * r * neutral_rate))], and necessary to evaluate exp(x) when x 
    # is very close to 0.
    ll_yes = yes_event[, sum(log(-1 * expm1(-1 * r * neutral_rate)))]
    ll_no = no_event[, -1 * sum(r * neutral_rate)]
    # if(is.nan(ll_no + ll_yes) | is.infinite(ll_no + ll_yes)) {
    #   browser()
    # }
    yes_event[, r := NULL]
    no_event[, r := NULL]
    return(-1 * (ll_no + ll_yes))
  }
  formals(lik)[["si"]] = rep(1.0, length(gene_names))
  
  if(length(gene_names) == 1) {
    bbmle::parnames(lik) = 'si'
    par_init = formals(lik)
    names(par_init) = bbmle::parnames(lik)
  } else {
    bbmle::parnames(lik) = gene_names
    par_init = setNames(formals(lik)[["si"]], gene_names)
  }
  
  fit = bbmle::mle2(lik, start = par_init, vecpar = T, lower = 1e-6, upper = 1e9, method = 'L-BFGS-B')
  if(length(coef(fit)) == 1) {
    selection_loci[, si := coef(fit)['si']]
  } else {
    selection_loci[, si := coef(fit)[gene]]
  }
  selection_loci[, num_samples_used := num_samples_used]
  
  if(debug == TRUE) {
     return(list(si = selection_loci[],
                 yes = original_yes_event,
                 no = original_no_event))
  }
  return(selection_loci[])
}

## Deprecated
# prep_seg_lik = function(gene_name, cna_calls = NULL, event_type = 'increase', seg_rates,
#                         gene_coord = ces.refset.hg38$cancer_gene_coord) {
#   total_rates = copy(seg_rates$total_rates)
#   seg_rates = copy(seg_rates$seg_rates)
#   setkey(cna_calls, chr, start, end)
#   
#   # This will be moved into CES internal data later.
#   sr_info = data.table(pretty_label = c('L < 100 kb', '100 kb < L < 1 Mb', '1 Mb < L < 10 Mb', '10 Mb < L < 40 Mb', '40 Mb < L'),
#                        cosmic_label = c('0-100kb', '100kb-1Mb', '1Mb-10Mb', '10Mb-40Mb', '>40Mb'),
#                        bin_id = c('r1', 'r2', 'r3', 'r4', 'r5'),
#                        bin_order = 1:5)
#   
#   groups = sr_info$cosmic_label # size ranges
#   
#   stopifnot(all(c('sample', 'size_range', 'seg_rate') %in% names(seg_rates)))
#   curr_rates = seg_rates[size_range %in% groups]
#   
#   if(curr_rates[, .N] == 0) {
#     warning('No rates found; returning empty table')
#     return(NULL)
#   }
#   
#   gr = gene_coord[gene == gene_name, .(chr, start, end)]
#   setkey(gr, chr, start, end)
#   curr_ol = foverlaps(cna_calls, gr, nomatch = NULL)
#   
#   curr_samples = unique(curr_ol$sample)
#   curr_n = length(curr_samples)
#   
#   curr_ol[curr_rates, neutral_rate := seg_rate, on = c('sample', 'size_range')]
#   
#   #neutral_prob = curr_rates[, .(exp = (rate/n_nondiploid) * curr_n), by = 'size_range']
#   
#   
#   # If more than one event in a sample, we will count the largest event
#   if(event_type == 'decrease') {
#     curr_ol = curr_ol[total_copy < exp_total_copy]
#   } else if(event_type == 'loh') {
#     curr_ol = curr_ol[nMinor == 0] # Yeah, total_copy might increase! Yet something was lost.
#   } else if(event_type == 'increase') {
#     curr_ol = curr_ol[total_copy > exp_total_copy]
#   }
#   
#   neutral_prob = unique(curr_ol[, .(sample, size_range, neutral_rate)]) # sample might have two same-size events
#   neutral_prob[sr_info, bin_order := bin_order, on = c(size_range = 'cosmic_label')]
#   
#   # Get rid of smaller size entries where samples have multiple entries 
#   neutral_prob = unique(neutral_prob[order(-bin_order)], by = 'sample')
#   
#   samples_with = list()
#   
#   # We go from largest to smallest bin
#   for(i in sort(sr_info$bin_order, decreasing = TRUE)) {
#     curr_bin = sr_info[bin_order == i, bin_id]
#     curr_size_range = sr_info[bin_order == i, cosmic_label]
#     
#     curr_samples_with = curr_ol[size_range == curr_size_range, unique(sample)]
#     
#     # Larger events in the same samples take precedence (already handled above, really)
#     samples_with[[curr_bin]] = setdiff(curr_samples_with, unlist(samples_with))
#   }
#   samples_with$with_any = unlist(samples_with)
#   samples_with$with_none = setdiff(curr_samples, samples_with$with_any)
#   
#   
#   # for lik
#   combined_rate = setNames(total_rates$seg_rate, total_rates$sample)[curr_samples]
#   with_none = samples_with$with_none
#   neutral_prob[, total_rate := combined_rate[sample]]
#   cond_prob = setNames(neutral_prob[, neutral_rate/total_rate], neutral_prob$sample)
#   
#   lik = function(gamma) {
#     sum_log_lik = 0
#     if (length(with_none) > 0) {
#       sum_log_lik = -1 * sum(gamma * combined_rate[with_none])
#     }
#     
#     for(bin_id in sr_info$bin_id) {
#       curr_samples = samples_with[[bin_id]]
#       curr_combined_rate = combined_rate[curr_samples] * gamma
#       prob_any = 1 - exp(-curr_combined_rate)
#       prob_curr = prob_any * cond_prob[curr_samples]
#       sum_log_lik = sum_log_lik + sum(log(prob_curr))
#     }
#     return(-1 * sum_log_lik)
#   }
#   formals(lik)[["gamma"]] = 1
#   bbmle::parnames(lik) = "gamma"
#   gene_interval = gr[, .(chr = chr[1], start = min(start), end = max(end))]
#   output = list(lik = lik, 
#                 groups = samples_with,
#                 num_covering = length(curr_samples), n_neutral_dec = sum(curr_rates$rate),
#                 interval = gene_interval,
#                 gene_name = gene_name)
#   return(output)
# }





