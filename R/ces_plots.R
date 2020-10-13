#' Plot selection intensities
#' 
#' Lollipop plots can be used to compare highly selected variants
#' within or among groups of samples.
#' 
#' To keep the lollipops readable, no more than \code{max_sites} variants are shown for
#' any group. At the same time, all SIs that are within the limits of the output plot are
#' depicted (which means that specification of \code{max_sites} creates an implicit SI
#' floor). If the labels don't look good, you may have to reduce max_sites or break 
#' up the plot into pieces. (If you find a better way to get lots of variant labels to
#' display nicely, we would love to hear about it!)
#' 
#' 
#' @import ggplot2
#' @import ggrepel
#' @param si_list Named list of tables, usually created by subsetting from selection
#'   tables. Each table should have columns called "variant_name", "variant_type", and "gene",
#'   plus a selection intensity column
#' @param si_col Vector giving names of selection intensity column in each SI table,
#'               either 1-length or same length as si_list
#' @param max_sites Maximum number of variant sites to include per lollipop; if you try to
#'   include too many, you may have a challenge getting it to look good
#' @param title Plot title to pass to ggplot
#' @param ylab Y-axis label to pass to ggplot
#' @param label_size Text size for labels, either 1-length or same length as si_list.
#' @param merge_dist How close points must be to be eligible to have their labels combined
#'   (.04 = 4% of plot space). Either 1-length or same length as si_list. Try tweaking if
#'   labels are not looking good. Set to zero to prevent any labels being combined.
#' @return ggplot object with lollipops
#' @examples 
#' \dontrun{
#' # Compare cancer subtypes from different analyses
#' lollipops(si_list = list(luad = luad_cesa$selection, lusc = lusc_cesa$selection))
#' 
#' # Compare results for two genes
#' tp53 = cesa$selection[gene == 'TP53']
#' egfr = cesa$selection[gene == 'EGFR']
#' lollipops(si_list = list(EGFR = egfr, TP53 = tp53))
#' }
#' @export
lollipops = function(si_list, si_col = "selection_intensity", title = "My SIs", ylab = "selection intensity", 
                     max_sites = 50, label_size = 3.0, merge_dist = .04) {
  if (! is(si_list, "list") || ! all(sapply(si_list, is.data.table))) {
    stop("si_list should be a named list of data.tables.")
  }
  num_groups = length(si_list)
  group_names = names(si_list)
  
  if(is.null(group_names)) {
    stop("si_list should be a named list (for labeling lollipops)")
  }
  
  if(! is(title, "character") | length(title) > 1) {
    stop("title should be 1-length character")
  }
  
  if(! is(max_sites, "numeric")) {
    stop("max_sites should be numeric")
  }
  
  if(! is(si_col, "character")) {
    stop("si_col should be character")
  }
  if (length(si_col) == 1) {
    si_col = rep(si_col, length(si_list))
  } else if (length(si_col) != length(si_list)) {
    stop("si_col should be 1-length or same length as si_list")
  }
  
  
  if (! is(label_size, "numeric") || ! length(label_size) %in% c(1, num_groups) || any(label_size < 0)) {
    stop("Illegal label_size value.")
  }
  if (length(label_size) < num_groups) {
    label_size = rep(label_size, num_groups)
  }
  
  if (! is(merge_dist, "numeric") || ! length(merge_dist) %in% c(1, num_groups) || any(merge_dist < 0) || any(merge_dist > 1)) {
    stop("Illegal merge_dist value.")
  }
  if (length(merge_dist) < num_groups) {
    merge_dist = rep(merge_dist, num_groups)
  }
  
  
  # for each group, get the highest N SIs (set by max_sites), then pick the
  # greatest of all of these as the SI floor so that no more than max_sites SIs print per lollipop
  min_si = 0
  for (i in 1:num_groups) {
    curr = copy(si_list[[i]])
    curr_cols = colnames(curr)
    curr_si_col = si_col[i]
    if (! curr_si_col %in% curr_cols) {
      stop("Selection intensity column not found in input table ", i, " (set with si_col).")
    }
    if (! all(c("variant_name", "gene", "variant_type") %in% curr_cols)) {
      stop("Missing required columns in input table ", i, " (see help).")
    }
    setnames(curr, curr_si_col, "selection_intensity")
    curr = curr[, .(variant_name, variant_type, gene, selection_intensity)]
    
    if (curr[, .N] > max_sites) {
      lowest_passing = curr[order(selection_intensity, decreasing = T), selection_intensity[max_sites]]
      if (lowest_passing > min_si) {
        min_si = lowest_passing
      }
    }
    si_list[[i]] = curr
  }
  
  si_list = lapply(si_list, function(x) x[selection_intensity >= min_si])
  
  # variants from genes with multiple variants will get colors in output
  all_genes = na.omit(unlist(lapply(si_list, function(x) x$gene)))
  multi_variant_genes = names(which(table(all_genes) > 1))


  placeholder = data.table(group = numeric(), selection_intensity = numeric(), variant_name_print = character(),
                           display_color = factor())
  xlim_low = .8
  xlim_high = num_groups + 1
  lolli = ggplot(placeholder, mapping = aes(group, selection_intensity, label = variant_name_print, color = display_color)) + 
          scale_y_log10() + ylab(ylab) +
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            legend.position = "none",
            plot.title = element_text(size = 20, hjust = .5)
          ) + ggtitle(label = title) +  scale_colour_discrete() + 
          scale_x_continuous(breaks = 1:num_groups, labels = group_names, limits = c(xlim_low, xlim_high))
  for (i in 1:num_groups) {
    top_hits = si_list[[i]][order(selection_intensity, decreasing = T)]
    if (top_hits[, .N] == 0) {
      warning("Group ", i, " has nothing to show in the plot space.", call. = F)
      next
    }
    top_hits[, log_si := log10(selection_intensity)]
    top_hits[, variant_name_print := gsub('_', ' ', variant_name)]
    curr_label_size = label_size[i]
    curr_merge_dist = merge_dist[i]
    
    # leave out gene names if all AACs have same gene
    # otherwise, keep gene names and add them to SNVs as well
    if (top_hits[variant_type == "aac", .N] > 0 & top_hits[variant_type == "aac", length(unique(gene)) == 1]) {
      regex = paste0('^', top_hits[variant_type == "aac", gene], '')
      old_names = top_hits[variant_type == "aac", variant_name_print]
      new_names = mapply(function(r, n) gsub(r, '', n), regex, old_names)
      top_hits[variant_type == "aac", variant_name_print := new_names]
    }
    top_hits[variant_type == "snv" & ! is.na(gene), variant_name_print := paste0('(', gene, ') ', variant_name_print)] 
      # usually there will be an intergenic annotation column
  
    top_hits[, group := i]
    
    # since intergenic variants have no gene, specify display color
    top_hits[, display_color := NA_character_]
    top_hits[gene %in% multi_variant_genes, display_color := factor(gene)]

    
    # Shows genes in different colors if more than one
    multiple_genes = length(unique(top_hits$gene)) > 1
    
    # Greedily merge labels of tightly clustered points for readability
    log_si_max = max(top_hits$log_si)
    log_si_range = log_si_max - min(top_hits$log_si)
    label_data = top_hits
    num_labels = label_data[, .N]
    wrap_width = (3 / curr_label_size) * 42
    min_wrap_width = wrap_width * .6 
    if (num_labels > 9) {
      label_data[, merge_eligible := T]
      label_data[, screen_y := (log_si_max - log_si) / log_si_range]
      while(TRUE) {
        label_data[, num_char := nchar(variant_name_print)]
        label_data[, num_near := sum(abs(label_data$screen_y - screen_y) < curr_merge_dist & merge_eligible == T), by = "variant_name"]
        l_center = label_data[merge_eligible == T & num_near == max(num_near), which = T][1]
        curr_screen_y = label_data[l_center, screen_y]
        l_to_add = setdiff(label_data[abs(screen_y - curr_screen_y) < curr_merge_dist & merge_eligible == T, which = T], l_center)
        if (length(l_to_add) == 0) {
          break
        }
        
        labels_to_merge = sort(c(l_center, l_to_add))
        labels_text = label_data$variant_name_print[labels_to_merge]
        num_to_add = length(labels_text)
        
        combined_label = paste(labels_text, collapse = ", ")
        combined_label = strwrap(combined_label, width = wrap_width)
        num_lines = length(combined_label)

        if (num_lines > 1 & nchar(combined_label[num_lines]) < min_wrap_width) {
          while(length(combined_label) == num_lines) {
            num_to_add = num_to_add - 1
            labels_to_merge = labels_to_merge[1:num_to_add]
            labels_text = label_data$variant_name_print[labels_to_merge]
            combined_label = paste(labels_text, collapse = ", ")
            combined_label = strwrap(combined_label, width = wrap_width)
          }
        }
        combined_label = paste0(combined_label, collapse = "\n")
        label_data[labels_to_merge[1], variant_name_print := combined_label]
        label_data[labels_to_merge[1], merge_eligible := F]
        
        # When merging labels fails due to unworkable nchar, you can end up with num_to_add == 1
        if(num_to_add > 1) {
          to_drop = labels_to_merge[2:num_to_add]
          label_data = label_data[! to_drop]
        }
      }
    }
    
    # each group i gets SIs plotted on the x = i line, with labels allowed roughly in the (i, i + 1) space
    lolli = lolli + 
            geom_vline(xintercept = i, color = "grey") +
            geom_label_repel(data = label_data, nudge_x = .1, direction = 'both', 
                             segment.size = .3, size = curr_label_size, color = "black",  xlim = c(i, i + .98))
    # silly: scale_color_discrete fails on all NAs despite NA color option
    if (all(is.na(top_hits$display_color))) {
      lolli = lolli + geom_point(data = top_hits, alpha = .7, color = "gray50")
    } else {
      lolli = lolli + geom_point(data = top_hits, alpha = .7)
    }
  }
  return(lolli)
}


