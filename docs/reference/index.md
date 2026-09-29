# Package index

## Create and manage projects

Create, load, and save projects

- [`CESAnalysis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/CESAnalysis.md)
  : Create a cancereffectsizeR analysis
- [`load_cesa()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/load_cesa.md)
  : Load a previously saved CESAnalysis
- [`save_cesa()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/save_cesa.md)
  : Save a CESAnalysis in progress
- [`set_refset_dir()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/set_refset_dir.md)
  : Set reference data directory

## Obtain and prep MAF data

- [`preload_maf()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/preload_maf.md)
  : Read and verify MAF somatic mutation data
- [`check_sample_overlap()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/check_sample_overlap.md)
  : Catch duplicate samples
- [`get_TCGA_project_MAF()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/get_TCGA_project_MAF.md)
  : Get MAF data from TCGA cohort
- [`vcfs_to_maf_table()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/vcfs_to_maf_table.md)
  : Read a VCF into an MAF-like table
- [`make_PathScore_input()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/make_PathScore_input.md)
  : Make a PathScore input file from MAF data
- [`lift_bed()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/lift_bed.md)
  : Convert BED intervals between genome builds

## Load and manage variants

- [`load_maf()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/load_maf.md)
  : Load MAF somatic mutation data
- [`select_variants()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/select_variants.md)
  : Select and filter variants
- [`variant_counts()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/variant_counts.md)
  : Assess variant prevalence and coverage
- [`samples_with()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/samples_with.md)
  : Find samples with specified variants
- [`add_variants()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/add_variants.md)
  : Add variant annotations
- [`add_covered_regions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/add_covered_regions.md)
  : add_covered_regions
- [`baseline_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/baseline_mutation_rates.md)
  : Baseline mutation rate calculation

## Load sample-level data

- [`load_sample_data()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/load_sample_data.md)
  : Add sample data
- [`clear_sample_data()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/clear_sample_data.md)
  : Clear sample data

## Compound variants

Combine variants into arbitrary batches and test for batch-level
selection

- [`define_compound_variants()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/define_compound_variants.md)
  : Divide batches of variants into a CompoundVariantSet
- [`CompoundVariantSet()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/CompoundVariantSet.md)
  : Create CompoundVariantSet from variant IDs

## Trinucleotide signatures and rates

Mutational signature extraction and inference of context-specific
mutation rates

- [`trinuc_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/trinuc_mutation_rates.md)
  : Calculate relative rates of trinucleotide-context-specific mutations
  by extracting underlying mutational processes
- [`suggest_cosmic_signature_exclusions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/suggest_cosmic_signature_exclusions.md)
  : Tissue-specific mutational signature exclusions
- [`trinuc_snv_counts()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/trinuc_snv_counts.md)
  : Tabulate SNVs by trinucleotide context
- [`convert_signature_weights_for_mp()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/convert_signature_weights_for_mp.md)
  : Get MutationalPatterns contributions matrix
- [`clear_trinuc_rates_and_signatures()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/clear_trinuc_rates_and_signatures.md)
  : Clear mutational signature attributions and related mutation rate
  information
- [`set_signature_weights()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/set_signature_weights.md)
  : Set SNV signature weights
- [`set_trinuc_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/set_trinuc_rates.md)
  : Assign pre-calculated relative trinucleotide mutation rates
- [`assign_group_average_trinuc_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/assign_group_average_trinuc_rates.md)
  : Skip mutational signature analysis and assign group average relative
  trinucleotide-context-specific mutation rates to all samples

## Gene mutation rates

Calculate neutral gene mutation rates

- [`gene_mutation_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/gene_mutation_rates.md)
  : Use dNdScv with tissue-specific covariates to calculate gene-level
  mutation rates
- [`set_gene_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/set_gene_rates.md)
  : Assign pre-calculated regional mutation rates
- [`clear_gene_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/clear_gene_rates.md)
  : Clear regional mutation rates

## Cancer effect sizes

Quantify selection for somatic variants

- [`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md)
  : Calculate cancer effects of variants
- [`ces_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_epistasis.md)
  : Variant-level pairwise epistasis
- [`ces_gene_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_gene_epistasis.md)
  : Gene-level epistasis
- [`mutational_signature_effects()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/mutational_signature_effects.md)
  : Attribute cancer effects to mutational signatures
- [`clear_effect_output()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/clear_effect_output.md)
  : Clear variant effect output
- [`clear_epistasis_output()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/clear_epistasis_output.md)
  : Clear epistasis output

## Continuous covariate selection

Functions for estimating and visualizing selection as a function of a
continuous sample-level covariate.

- [`ces_variant_linear()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_linear.md)
  : Estimate variant effects with a linear continuous-covariate model
- [`ces_variant_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_logistic.md)
  : Estimate variant effects with a logistic continuous-covariate model
- [`ces_variant_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_sigmoid.md)
  : Estimate variant effects with a generalized sigmoid
  continuous-covariate model
- [`plot_effects_continuous()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_effects_continuous.md)
  : Plot continuous-covariate cancer effects
- [`sswm_age_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik.md)
  : Linear continuous-covariate SSWM likelihood
- [`sswm_age_lik_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_logistic.md)
  : Logistic continuous-covariate SSWM likelihood
- [`sswm_age_lik_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_sigmoid.md)
  : Generalized sigmoid continuous-covariate SSWM likelihood

## Step-specific selection

Quantify selection separately across ordered tumor progression stages

- [`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md)
  : Calculate cancer effects of variants across ordered progression
  stages
- [`assign_stage_index()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/assign_stage_index.md)
  : Assign samples to ordered progression stages
- [`stage_mutation_proportions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/stage_mutation_proportions.md)
  : Calculate per-gene, per-stage mutation rate proportions
- [`step_selection_LRT()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/step_selection_LRT.md)
  : Likelihood ratio test for step-specific selection
- [`plot_effects_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_effects_step.md)
  : Plot step-specific cancer effects

## Visualization

Display and compare variant effect sizes

- [`plot_effects()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_effects.md)
  : Plot cancer effects
- [`plot_signature_effects()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_signature_effects.md)
  : Plot mutational source and effect attributions
- [`plot_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_epistasis.md)
  : Plot pairwise epistasis
- [`epistasis_plot_schematic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/epistasis_plot_schematic.md)
  : Get epistatic effect schematic

## Selection models

Likelihood function generators for various models of selection

- [`sswm_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_lik.md)
  : sswm_lik
- [`step_selection_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/step_selection_lik.md)
  : step_selection_lik
- [`pairwise_epistasis_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/pairwise_epistasis_lik.md)
  : pairwise_epistasis_lik

## Explore reference data

- [`cosmic_signature_info()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/cosmic_signature_info.md)
  : Get COSMIC signature descriptions
- [`get_PathScore_coding_regions()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/get_PathScore_coding_regions.md)
  : Get PathScore coding regions
- [`list_ces_refsets()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/list_ces_refsets.md)
  : list_ces_refsets
- [`list_ces_covariates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/list_ces_covariates.md)
  : list_ces_covariates
- [`list_ces_signature_sets()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/list_ces_signature_sets.md)
  : list_ces_signature_sets
- [`get_ces_signature_set()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/get_ces_signature_set.md)
  : get_ces_signature_set

## Create custom reference data

Build your own reference data set for almost any genome or tissue type

- [`create_refset()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/create_refset.md)
  : Create a custom refset
- [`build_RefCDS()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/build_RefCDS.md)
  : cancereffectsizeR's RefCDS builder
- [`validate_signature_set()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/validate_signature_set.md)
  : validate_signature_set

## Accessors

Data accessors that you probably won’t need (use cesa\$maf,
cesa\$samples, etc. instead)

- [`maf_records()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/maf_records.md)
  : View data loaded into CESAnalysis
- [`excluded_maf_records()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/excluded_maf_records.md)
  : View excluded MAF data
- [`get_sample_info()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/get_sample_info.md)
  : View sample metadata
- [`get_trinuc_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/get_trinuc_rates.md)
  : Get estimated relative rates of trinucleotide-specific SNV mutation
- [`get_signature_weights()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/get_signature_weights.md)
  : Get table of signature attributions
- [`get_gene_rates()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/get_gene_rates.md)
  : Get table of neutral gene mutation rates
- [`snv_results()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/snv_results.md)
  : View results from ces_variant
- [`epistasis_results()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/epistasis_results.md)
  : View output from epistasis functions
