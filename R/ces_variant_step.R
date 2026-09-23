#' Calculate cancer effects of variants across ordered progression stages
#'
#' Like \code{ces_variant()}, but under a step-specific model of selection: instead of assuming
#' a variant's scaled selection coefficient is constant across all included samples,
#' \code{ces_variant_step()} allows it to differ across an ordered sequence of tumor progression
#' stages (for example, normal tissue and primary tumor, or an arbitrary number of stages).
#'
#' This function is a standalone counterpart to \code{ces_variant()} (it does not call it, and
#' \code{ces_variant()} is unaffected by anything here): the underlying model of selection is
#' fundamentally different (multiple stage-specific selection coefficients, rather than one), so
#' the mechanics of confidence interval calculation and output differ in ways not supported by
#' \code{ces_variant()}'s \code{model} argument.
#'
#' Setting up a run requires two additional pieces of information beyond a normal
#' \code{ces_variant()} call:
#' \itemize{
#'   \item Which stage each sample belongs to, in the form of a \code{sample_index} table (see
#'     \code{assign_stage_index()}, or supply \code{stage_col}/\code{stage_order} here and it will
#'     be built for you).
#'   \item For each variant's gene, what proportion of the gene's mutation rate is estimated to
#'     accumulate during each stage, in the form of a \code{stage_mut_prop} table (see
#'     \code{stage_mutation_proportions()}).
#' }
#' Because stage mutation-rate proportions are gene-level, each variant (or compound variant) run
#' through \code{ces_variant_step()} must map unambiguously to exactly one gene.
#'
#' The output selection table has one selection intensity column per stage, named
#' \code{si_<stage name>} (and, when \code{conf} is set, matching \code{ci_low_<conf*100>_si_<stage
#' name>}/\code{ci_high_<conf*100>_si_<stage name>} columns), rather than the single
#' \code{selection_intensity} column of \code{ces_variant()}. It does not include
#' \code{included_with_variant}, \code{included_total}, or \code{uncovered} columns, since
#' sample accounting is inherently stage-specific here. Use \code{step_selection_LRT()}
#' to test whether a step-specific model fits significantly better than a plain \code{ces_variant()}
#' run, and \code{plot_effects_step()} to visualize output.
#'
#' @param cesa CESAnalysis object
#' @param variants Which variants to estimate effects for, specified with a variant table such as
#'   from \code{[CESAnalysis]$variants} or \code{select_variants()}, or a \code{CompoundVariantSet}
#'   from \code{define_compound_variants()}. Defaults to all recurrent mutations; that is,
#'   \code{[CESAnalysis]$variants[maf_prevalence > 1]}. To include all variants, set to
#'   \code{[CESAnalysis]$variants}. Every variant (or compound variant) must map to exactly one
#'   gene present in \code{stage_mut_prop}.
#' @param stage_mut_prop A data.table of per-gene, per-stage mutation rate proportions, as
#'   produced by \code{stage_mutation_proportions()}. Its \code{p_<stage name>} columns must cover
#'   every stage named in \code{sample_index} (or implied by \code{stage_order}).
#' @param stage_col Name of a sample-level data column recording each sample's stage. Supply this
#'   and \code{stage_order} to have \code{sample_index} built automatically; alternatively, supply
#'   \code{sample_index} directly.
#' @param stage_order Values of \code{stage_col} in earliest-to-latest stage order (see
#'   \code{assign_stage_index()}); used only when \code{sample_index} is not supplied directly.
#' @param sample_index A pre-built data.table associating samples with stages, as produced by
#'   \code{assign_stage_index()}. Supply this, or \code{stage_col}/\code{stage_order}, but not both.
#' @param samples Which samples to include in inference. Defaults to all samples. Can be a vector
#'   of Unique_Patient_Identifiers, or a data.table containing rows from the CESAnalysis sample
#'   table. Every included sample must be covered by \code{sample_index}.
#' @param run_name Optionally, a name to identify the current run.
#' @param optimizer_args Named list of arguments to pass to the optimizer, bbmle::mle2. Defaults to
#'   bounding stage selection intensities to [1e-3, 1e9] via L-BFGS-B.
#' @param return_fit TRUE/FALSE (default FALSE): Embed model fit for each variant in a "fit"
#'   attribute of the selection results table. Use \code{attr(selection_table, 'fit')} to access
#'   the list of fitted models. Defaults to FALSE to save memory.
#' @param hold_out_same_gene_samples When finding likelihood of each variant, hold out samples that
#'   lack the variant but have any other mutations in the same gene. By default, TRUE when running
#'   with single variants, FALSE with a CompoundVariantSet.
#' @param cores Number of cores to use for processing variants in parallel (not useful for Windows
#'   systems).
#' @param conf Confidence interval width for stage selection intensities (NULL skips calculation,
#'   speeds runtime).
#' @return CESAnalysis object with selection results appended to the selection output list
#' @export
ces_variant_step = function(cesa = NULL,
                            variants = select_variants(cesa, min_freq = 2),
                            stage_mut_prop = NULL,
                            stage_col = NULL,
                            stage_order = NULL,
                            sample_index = NULL,
                            samples = character(),
                            run_name = "auto",
                            optimizer_args = list(method = 'L-BFGS-B', lower = 1e-3, upper = 1e9),
                            return_fit = FALSE,
                            hold_out_same_gene_samples = "auto",
                            cores = 1,
                            conf = .95)
{
  if (! is.numeric(cores) || length(cores) != 1 || cores - as.integer(cores) != 0 || cores < 1) {
    stop('cores should be 1-length positive integer')
  }
  if (! rlang::is_bool(return_fit)) {
    stop('return_fit should be TRUE/FALSE.')
  }
  validate_optimizer_args(optimizer_args)
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (length(hold_out_same_gene_samples) == 1) {
    if (! is.logical(hold_out_same_gene_samples)) {
      if (identical(hold_out_same_gene_samples, "auto")) {
        hold_out_same_gene_samples = is(variants, "data.table")
      } else {
        stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
      }
    }
  } else {
    stop("hold_out_same_gene_samples should be TRUE/FALSE or left \"auto\".")
  }

  if (! is(run_name, "character") || length(run_name) != 1) {
    stop("run_name should be 1-length character")
  }
  if (run_name %in% names(cesa@selection_results)) {
    stop("The run_name you chose has already been used. Please pick a new one.")
  }
  if (! grepl('^[a-z]', tolower(run_name), perl = T) || grepl('\\s\\s', run_name)) {
    stop("Invalid run name. The name must start with a latter and contain no consecutive spaces.")
  }
  if (run_name == "auto") {
    run_number = length(cesa@selection_results) + 1
    run_name = paste0('selection.', run_number)
    while (run_name %in% names(cesa@selection_results)) {
      run_number = run_number + 1
      run_name = paste0('variant_effects_', run_number)
    }
  }

  if (! is.null(conf)) {
    if (! is(conf, "numeric") || length(conf) > 1 || conf <= 0 || conf >= 1) {
      stop("conf should be 1-length numeric (e.g., .95 for 95% confidence intervals), or NULL to skip CI calculation.", call. = F)
    }
  }

  samples = select_samples(cesa, samples)
  if (samples[, .N] < cesa@samples[, .N]) {
    num_excluded = cesa@samples[, .N] - samples[, .N]
    pretty_message(paste0("Note that ", format(num_excluded, big.mark = ','), " samples are being excluded from selection inference."))
  }

  # Build or validate sample_index, which assigns every included sample to an ordered stage
  if (is.null(sample_index)) {
    if (is.null(stage_col) || is.null(stage_order)) {
      stop("Supply either sample_index directly, or both stage_col and stage_order so it can be ",
           "built automatically (see assign_stage_index()).")
    }
    sample_index = assign_stage_index(cesa = cesa, stage_col = stage_col, stage_order = stage_order, samples = samples)
  } else {
    if (! is.null(stage_col) || ! is.null(stage_order)) {
      stop("Supply sample_index, or stage_col/stage_order, but not both.")
    }
    if (! is(sample_index, "data.table") ||
        ! all(c("Unique_Patient_Identifier", "group_index", "group_name") %in% names(sample_index))) {
      stop("sample_index should be a data.table with columns Unique_Patient_Identifier, group_index, ",
           "and group_name (see assign_stage_index()).")
    }
  }
  num_stages = sample_index[, uniqueN(group_index)]
  if (num_stages < 2) {
    stop("sample_index describes fewer than 2 stages; step-specific selection requires at least 2.")
  }
  missing_from_index = setdiff(samples$Unique_Patient_Identifier, sample_index$Unique_Patient_Identifier)
  if (length(missing_from_index) > 0) {
    stop(length(missing_from_index), " sample(s) being used for inference are missing from sample_index.")
  }

  # Validate stage_mut_prop, the per-gene table of stage-specific mutation rate proportions
  if (! is(stage_mut_prop, "data.table")) {
    stop("stage_mut_prop should be a data.table, as produced by stage_mutation_proportions().")
  }
  if ("gene" %in% names(stage_mut_prop)) {
    id_col = "gene"
  } else if ("pid" %in% names(stage_mut_prop)) {
    id_col = "pid"
  } else {
    stop("stage_mut_prop should have a gene or pid identifier column (see stage_mutation_proportions()).")
  }
  # Match stage proportion columns to sample_index's stages by name, not just position, to avoid
  # silent misalignment if the two tables were built with different stage orderings.
  stage_names_in_order = unique(sample_index[order(group_index)], by = "group_index")$group_name
  prop_cols = paste0("p_", stage_names_in_order)
  missing_prop_cols = setdiff(prop_cols, names(stage_mut_prop))
  if (length(missing_prop_cols) > 0) {
    stop("stage_mut_prop is missing column(s) matching sample_index's stage names: ",
         paste(missing_prop_cols, collapse = ", "), '. Column names should be "p_<stage name>" ',
         "(see stage_mutation_proportions()).")
  }
  stage_mut_prop = data.table::copy(stage_mut_prop)
  data.table::setkeyv(stage_mut_prop, id_col)

  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  mutations = cesa@mutations

  running_compound = FALSE
  if (is(variants, "data.table")) {
    if (! "variant_id" %in% names(variants)) {
      stop("variants table is missing a variant_id column. Typically, variants is generated using select_variants().")
    }
    nonoverlapping = attr(variants, "nonoverlapping")
    if (is.null(nonoverlapping)) {
      if ('variant_id' %in% names(variants)) {
        pretty_message('Taking variants from variant_id column of input table....')
      }
    } else if (! identical(nonoverlapping, TRUE)) {
      stop("Input variants table may contain overlapping variants; re-run select_variants() to get a non-overlapping table.")
    }
    variants = select_variants(cesa, variant_ids = variants[, variant_id])
  } else if (is(variants, "CompoundVariantSet")) {
    running_compound = TRUE
    if (cesa@advanced$uid != variants@cesa_uid) {
      stop("Input CompoundVariantSet does not appear to derive from the input CESAnalysis.")
    }
    if (cesa@samples[, .N] != variants@cesa_num_samples) {
      stop("The number of samples in the CESAnalysis has changed since the CompoundVariantSet was created. ",
           "Please re-generate it.")
    }
    compound_variants = variants
    variants = select_variants(cesa, variant_ids = compound_variants@snvs$snv_id, include_subvariants = TRUE)
    variants = variants[compound_variants@snvs, compound_name := compound_name, on = c(variant_id = "snv_id")]
    if (variants[, .N] != compound_variants@snvs[, .N]) {
      stop("Internal error: select_variants() didn't return variant info 1-to-1 with compound input.")
    }
    variants[compound_variants@compounds, covered_in := shared_cov, on = "compound_name"]
  } else {
    stop("variants expected to be a variant table (from select_variants(), usually) or a CompoundVariantSet")
  }
  if (variants[, .N] == 0) {
    stop("There are no variants in the input!")
  }
  aac_ids = variants[variant_type == "aac", variant_id]
  noncoding_snv_ids = variants[variant_type == "snv", variant_id]
  if (length(aac_ids) + length(noncoding_snv_ids) == 0) {
    stop("No variants pass filters, so there are no SIs to calculate.", call. = F)
  }

  maf = cesa@maf[samples$Unique_Patient_Identifier, on = "Unique_Patient_Identifier", nomatch = NULL]
  tmp = unique(maf[, .(gene = unlist(genes)), by = "Unique_Patient_Identifier"])[, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene = tmp$samples
  names(tumors_with_variants_by_gene) = tmp$gene
  tumors_with_variants_by_gene = list2env(tumors_with_variants_by_gene)

  snv_aac_of_interest = cesa@mutations$aac_snv_key[aac_ids, on = 'aac_id']
  tmp = maf[snv_aac_of_interest, .(variant_id), on = c(variant_id = 'snv_id'),
            by = "Unique_Patient_Identifier", nomatch = NULL]
  tmp[snv_aac_of_interest, aac_id := aac_id, on = c(variant_id = 'snv_id')]
  tmp = tmp[, .(samples = list(unique(Unique_Patient_Identifier))), by = "aac_id"]
  samples_by_aac = setNames(tmp$samples, tmp$aac_id)

  setkey(maf, "variant_id")
  tmp = maf[noncoding_snv_ids, variant_id, by = "Unique_Patient_Identifier", nomatch = NULL][, .(samples = list(Unique_Patient_Identifier)), by = "variant_id"]
  samples_by_snv = tmp$samples
  names(samples_by_snv) = tmp$variant_id
  samples_by_variant = list2env(c(samples_by_aac, samples_by_snv))

  setkey(samples, "covered_regions")
  genome_wide_cov_samples = samples["genome", Unique_Patient_Identifier, nomatch = NULL]

  selection_results = NULL
  all_coverage = rbind(cesa@mutations$snv[, .(variant_id = snv_id, covered_in)],
                       cesa@mutations$amino_acid_change[, .(variant_id = aac_id, covered_in)])
  variants[all_coverage, covered_in := covered_in, on = 'variant_id']
  coverage_groups = unique(variants$covered_in)
  num_coverage_groups = length(coverage_groups)

  setkey(maf, "variant_id")
  setkey(variants, "variant_id")

  selection_results = lapply(1:length(coverage_groups), function(i) {
    coverage_group = coverage_groups[[i]]
    if (length(coverage_group) == 1 && is.na(coverage_group)) {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, NA_character_)))]
    } else {
      curr_variants = variants[which(sapply(variants$covered_in, function(x) identical(x, coverage_group)))]
    }
    message(sprintf("Preparing to calculate cancer effects (batch %i of %i)...", i, num_coverage_groups))
    if (is.null(coverage_group)) {
      covered_samples = genome_wide_cov_samples
    } else {
      covered_samples = c(samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL], genome_wide_cov_samples)
    }
    variants[curr_variants$variant_id, num_covered_and_in_samples := length(covered_samples), on = 'variant_id']

    if (length(covered_samples) == 0) {
      return(list(data.table(), NULL))
    }
    work_size = length(covered_samples) * curr_variants[, .N] * 8
    num_proc_groups = ceiling(work_size / 1e9)
    curr_variants[, subgroup := ceiling(num_proc_groups * 1:.N / .N)]

    if (running_compound) {
      curr_variants[, subgroup := rep.int(subgroup[1], .N), by = "compound_name"]
      num_proc_groups = max(curr_variants$subgroup)
    }

    curr_results = lapply(1:num_proc_groups, function(j) {
      if (num_proc_groups > 1) {
        message(sprintf("Working on sub-batch %i of %i...", j, num_proc_groups))
      }
      curr_subgroup = curr_variants[subgroup == j]
      aac_ids = curr_subgroup[variant_type == "aac", variant_id]
      snv_ids = curr_subgroup[variant_type == "snv", variant_id]

      baseline_rates = baseline_mutation_rates(cesa, aac_ids = aac_ids, snv_ids = snv_ids, samples = covered_samples)
      gene_lookup = curr_subgroup[, all_genes]
      names(gene_lookup) = curr_subgroup[, variant_id]
      gene_lookup = list2env(gene_lookup)

      process_variant = function(variant_id) {
        if (running_compound) {
          compound_id = variant_id
          tumors_with_variant = intersect(compound_variants@sample_calls[[compound_id]], covered_samples)
          current_snvs = compound_variants@snvs[compound_name == compound_id]
          all_genes = current_snvs[, unique(unlist(genes))]
          variant_id = current_snvs$snv_id
          rates = baseline_rates[, ..variant_id]
          rates = rowSums(rates)
        } else {
          tumors_with_variant = samples_by_variant[[variant_id]]
          all_genes = gene_lookup[[variant_id]]
          rates = baseline_rates[, ..variant_id][[1]]
        }
        names(rates) = baseline_rates[, Unique_Patient_Identifier]

        if (hold_out_same_gene_samples) {
          if (length(all_genes) == 1) {
            if (is.na(all_genes)) {
              tumors_with_gene_mutated = tumors_with_variant
            } else {
              tumors_with_gene_mutated = tumors_with_variants_by_gene[[all_genes]]
            }
          } else {
            tumors_with_gene_mutated = unique(unlist(sapply(all_genes, function(x) tumors_with_variants_by_gene[[x]])))
          }
          tumors_without = setdiff(covered_samples, tumors_with_gene_mutated)
        } else {
          tumors_without = setdiff(covered_samples, tumors_with_variant)
        }
        rates_tumors_with = rates[tumors_with_variant]
        rates_tumors_without = rates[tumors_without]

        # Step-specific selection requires exactly one gene per variant, so its stage_mut_prop
        # row can be looked up unambiguously.
        display_id = if (running_compound) compound_id else variant_id
        if (length(all_genes) != 1 || anyNA(all_genes)) {
          stop("ces_variant_step() requires each variant (or compound variant) to map to exactly ",
               "one gene, so its stage_mut_prop row can be identified; variant/compound \"",
               display_id[1], "\" maps to ", length(all_genes), " gene(s).", call. = FALSE)
        }
        this_prop_row = stage_mut_prop[all_genes, on = id_col, nomatch = NULL]
        if (this_prop_row[, .N] == 0) {
          stop("Gene \"", all_genes, "\" (for variant/compound \"", display_id[1],
               "\") is not present in stage_mut_prop.", call. = FALSE)
        }
        this_stage_mut_prop = as.numeric(unlist(this_prop_row[1, ..prop_cols], use.names = FALSE))
        if (anyNA(this_stage_mut_prop)) {
          stop("Gene \"", all_genes, "\" (for variant/compound \"", display_id[1],
               "\") has NA stage_mut_prop values.", call. = FALSE)
        }

        lik_args = list(sample_index = sample_index, stage_mut_prop = this_stage_mut_prop,
                        rates_tumors_with = rates_tumors_with, rates_tumors_without = rates_tumors_without)
        fn = do.call(step_selection_lik, lik_args)
        par_init = formals(fn)[[1]]
        names(par_init) = bbmle::parnames(fn)
        final_optimizer_args = c(list(minuslogl = fn, start = par_init, vecpar = T), optimizer_args)

        withCallingHandlers(
          {
            fit = do.call(bbmle::mle2, final_optimizer_args)
          },
          warning = function(w) {
            if (startsWith(conditionMessage(w), "some parameters are on the boundary")) {
              invokeRestart("muffleWarning")
            }
            if (grepl(x = conditionMessage(w), pattern = "convergence failure")) {
              invokeRestart("muffleWarning")
            }
          }
        )

        selection_intensity = bbmle::coef(fit)
        loglikelihood = as.numeric(bbmle::logLik(fit))

        if (running_compound) {
          variant_id = compound_id
        }
        variant_output = c(list(variant_id = variant_id),
                           as.list(selection_intensity),
                           list(loglikelihood = loglikelihood))

        if (! is.null(conf)) {
          min_value = ifelse(is.null(final_optimizer_args$lower), -Inf, final_optimizer_args$lower)
          max_value = ifelse(is.null(final_optimizer_args$upper), Inf, final_optimizer_args$upper)
          variant_output = c(variant_output,
                             univariate_si_conf_ints(fit, fn, min_value, max_value, conf))
        }
        return(list(summary = variant_output, fit = if (return_fit) fit else NULL))
      }

      if (running_compound) {
        variants_to_run = curr_subgroup[, unique(compound_name)]
      } else {
        variants_to_run = curr_subgroup$variant_id
      }
      message("Calculating cancer effects...")
      subgroup_results = pbapply::pblapply(variants_to_run, process_variant, cl = cores)
      subgroup_selection = rbindlist(lapply(subgroup_results, '[[', 1))
      subgroup_fit = lapply(subgroup_results, '[[', 2)
      return(list(subgroup_selection, subgroup_fit))
    })
    group_selection = rbindlist(lapply(curr_results, '[[', 1))
    group_fit = lapply(curr_results, '[[', 2)
    return(list(group_selection, group_fit))
  })

  fits = unlist(lapply(selection_results, '[[', 2))
  selection_results = rbindlist(lapply(selection_results, '[[', 1))

  if (selection_results[, .N] == 0) {
    msg = paste0("No selection inference was performed, so returning CESAnalysis unaltered without selection output. ",
                 "Perhaps none of the variants had coverage in the specified samples?")
    pretty_message(msg)
    return(cesa)
  }

  if (running_compound) {
    selection_results[, variant_type := "compound"]
    setnames(selection_results, "variant_id", "variant_name")
    selection_results = selection_results[compound_variants@compounds, on = c(variant_name = "compound_name")]
    setattr(selection_results, "is_compound", TRUE)
    setcolorder(selection_results, c("variant_name", "variant_type"))
  } else {
    selection_results[variants, c("variant_type", "variant_name", "gene", "intergenic") :=
                        list(variant_type, variant_name, gene, intergenic), on = "variant_id"]
    selection_results[intergenic == T, gene := NA]
    selection_results$intergenic = NULL
    setattr(selection_results, "is_compound", FALSE)
    setcolorder(selection_results, c("variant_name", "variant_type", "gene"))
    setcolorder(selection_results, c(setdiff(names(selection_results), 'variant_id'), 'variant_id'))
  }

  if (hold_out_same_gene_samples == FALSE) {
    selection_results[, held_out := 0]
  }

  if (return_fit) {
    fits = lapply(fits, function(x) {
      x@call.orig = call('[not shown]')
      parent.env(environment(x)) = emptyenv()
      return(x)
    })
    setattr(selection_results, 'fit', fits)
  }
  curr_results = list(selection_results)
  names(curr_results) = run_name

  cesa@selection_results = c(cesa@selection_results, curr_results)
  return(cesa)
}
