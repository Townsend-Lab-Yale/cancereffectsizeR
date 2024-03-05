
# <font style = "opacity:0">cancereffectsizeR 4.0.0</font>
Patch releases (as in, x.y.1 &#8594; x.y.2) have minor bug fixes or small improvements that do not significantly affect the numerical output of cancer effect analyses. Minor/major updates may change some outputs due to bug fixes or methodological tweaks, as described in these version notes.<br><br>

# cancereffectsizeR 3.0.0
* Doublet base substitutions (DBS) are included in variant annotation tables.
* select_variants() and variant_counts() are simpler and more powerful.
* Breaking changes: Pre-3.0 analyses can only be loaded in a read-only summary format. Some table column names and function arguments have changed throughout the package, so scripts may need to be updated.

# cancereffectsizeR 2.9.0
* plot_signature_effects() visualizes the relative contributions of mutational signatures to mutation and selection.
* Change to how mutational_signature_effects() calculates cohort-averaged signature effect shares. See the function's updated documentation for clarification of how outputs are calculated.
* When available, canonical transcripts (MANE transcripts) are now favored in determining the consequences of substitutions. Alternate transcripts will still be favored if they show more dramatic alterations. For example, if a substitution is a missense on the canonical transcript and splice-site-disrupting on an alternate transcript, the alternate transcript will be reported in the [CESAnalysis]$variants table. Additionally, a new annotation column (is_MANE) identifies if the transcript being reported is canonical. These features require the new reference data package release, ces.refset.hg38 v1.3.0.

# cancereffectsizeR 2.8.1
* get_PathScore_coding_regions() provides access to the CDS definitions used by PathScore (and by make_PathScore_input()).
* Various small fixes.

# cancereffectsizeR 2.8.0
* plot_effects() visualizes cancer effect inferences with custom labels, colors, variant groupings, and more; see website for examples.
* make_PathScore_input() converts MAF data into an input file for PathScore, a web tool that identifies significantly mutated pathways in cancer. See https://pathscore.publichealth.yale.edu/.
* ces_variant() gains the argument optimizer_args, giving you precise control over how
selection models are fit. This is intended for use with custom models of selection.

# cancereffectsizeR 2.7.0
* mutational_signature_effects() identifies the relative contributions of mutational signatures to mutation burden and cancer effect, within samples and across cohorts.
* Epistatic effect inferences now report p-values for overall significance of epistatic effects (over a null model of no epistasis) and significance of epistatic effects for each variant/gene in tested pairs. Additionally, output has been beautified with clearer column names and additional documentation.
* Removed deprecated "sample_groups" argument from CESAnalysis(). Package functions allow flexible sample grouping via other means.

# cancereffectsizeR 2.6.5
* The "nonsilent" method in ces_gene_epistasis() now uses nonsynonymous coding variants and essential splice-site variants (coding or not). Previously, noncoding splice-site variants weren't included.
* In epistasis output, the fitted models for each inference are now available for exploration. See the `return_fit`
option in ces_epistasis() and ces_gene_epistasis().
* Updated tutorial.

# cancereffectsizeR 2.6.4
* load_maf(): Optional maf_name argument makes it easy to see which samples in a CESAnalysis come from which MAF.
* More flexible support for user-supplied gene mutation rates.
* Changed par_init for selection inference on ces_variant()'s default model. This change will cause trivially small differences in effect estimates from prior versions.
* Various other small fixes.

# cancereffectsizeR 2.6.3
* vcfs_to_maf_table(): Read in VCF files from somatic variant calling pipelines and convert into a combined MAF-like data.table.
* Removed sequential model of selection pending improvements to methodology.
* Improved user guidance for get_TCGA_project_MAF().

# cancereffectsizeR 2.6.2
* Miscellaneous minor improvements, including better handling of chromosome naming styles.

# cancereffectsizeR 2.6.1
* Compatibility fix for Windows systems.
* Bug fixes for features introduced in 2.6.0.

# cancereffectsizeR 2.6.0
* get_TCGA_project_MAF(): Conveniently download MAF data from the Genomic Data Commons and create MAF files covering full TCGA projects.
* Change to relative trinucleotide mutation rate calculation for low-mutation samples for consistency with other samples. Impact of the change will typically be small.
* ces_gene_epistasis() now accepts a list of gene pairs to analyze.
* Conversion of adjacent SNVs to DBS/MNV has moved from load_maf() to preload_maf(), and is now optional (but typically recommended).

# cancereffectsizeR 2.5.0
* Fixed bug where combining WGS data with generic WXS data (without specific coverage definitions) resulted in all variants being annotated as exonic.

# cancereffectsizeR 2.4.0
* ces_gene_epistasis() gains a variants argument to customize which variants in each gene are used in epistasis inference.
* Gene mutation rate estimates now come with 95% confidence intervals.
* add_variants() can now be used to load annotations for amino-acid-changing substitutions into a CESAnalysis.
* Fixed a bug in create_refset() that caused custom reference data sets to have incorrect values in gene trinucleotide proportions (gene_trinuc_comp).
* Added utility functions clear_effect_output() and clear_epistasis_output().

# cancereffectsizeR 2.3.4
* New website and expanded tutorial.
* Small fixes, including a tweak to variant prioritization.

# cancereffectsizeR 2.3.3
* New quickstart instructions on the website, and a new tutorial in progress. See the site for the latest version, as it will continue to be updated.
* For simplicity, load_maf() no longer supports non-standard MAF column names or genome assembly conversion. Instead,call preload_maf() first to prep data.
* Renamed function: suggest_cosmic_signature_exclusions() replaces suggest_cosmic_signatures_to_remove(), and the related argument to trinuc_mutation_rates() is now signature_exclusions.
* Small fixes and tweaks.

# cancereffectsizeR 2.3.1
* ces_variant() output includes, for each input variant, the number of samples included in the inference and the number of those samples with the variant. (For various reasons, including the default exclusion of samples with same-gene variants, the total number of samples in a given selection inference can be less than the number of samples in the data set.)
* Miscellaneous small improvements and fixes.

# cancereffectsizeR 2.3.0
* Load sample-level data into a CESAnalysis sample table using add_sample_data(), or during load_maf() with the new sample_data_cols argument.
* trinuc_mutation_rates(), gene_mutation_rates(), and selection inference functions can be run on arbitrary subsets of samples. The less-flexible CESAnalysis "sample_groups" functionality has been deprecated.
* samples_with() makes it easy to see which samples have various mutations.
* variant_counts() provides variant prevalence and coverage information, with the option to break counts down into groups based on sample table columns.  
* check_sample_overlap() now accepts a list of MAFs.
* In loaded MAF data, columns top_consequence and top_gene give the most significant annotated coding changes for each mutation record. Annotation precedence is determined by MAF prevalence (usually equal), essential splice status, premature stop codon, nonsilent status, MAF mutation prevalence across the transcript (often favors longer transcripts), and finally alphabetical order. The columns are recalculated when more data is loaded, so changes in MAF prevalence can change which variants appear.
* Because ces.refset.hg38 is now available, CESAnalysis() and preload_maf() now require reference data sets to be specified by the user (formerly, they defaulted to ces.refset.hg19).

# cancereffectsizeR 2.2.2
* View dNdScv results more easily, and support for passing custom parameters to dNdScv.
* Simplified signature output to two tables: "raw attributions" and "biological weights" (see docs).
* Miscellaneous small fixes and improvements.

# cancereffectsizeR 2.2.1
* Reorganization of signature weights data to raw_weights, adjusted_weights, biological_weights. Raw weights come straight from signature extraction.
* Variant annotation optimized for better performance on large data sets.
* Fixes for preload_maf/load_maf usability issues.

# cancereffectsizeR 2.2.0
* Support for using MutationalPatterns for signature extraction. The MutationalPatterns fit_to_signatures_strict method is now the default used by trinuc_mutation_rates().
* Helper function trinuc_snv_counts eases exploratory signature analysis, and the counts used internally by trinuc_mutation_rates() are now viewable.
* convert_signature_weights_for_mp() makes it easy to input mutational signature weights into MutationalPatterns functions, such as visualizations.
* Custom reference data sets can now contain multiple transcripts per gene.
* Fixed bug that that inflated selection intensities of noncoding SNVs when using default model/workflow.

# cancereffectsizeR 2.1.4
* Minor fixes and improvements.

# cancereffectsizeR 2.1.3
* Improvements to liftOver support.
* preload_maf() gives more information about problematic records.

# cancereffectsizeR 2.1.2
* load_maf() detects adjacent SNV records and re-annotates as doublet base substitutions; improved annotation of insertions, deletions, and more complex multi-nucleotide variants.
* Performance improvements and minor bug fixes.

# cancereffectsizeR 2.1.1
* Minor improvements and documentation updates.


# cancereffectsizeR 2.1.0
* Simple workflow for building custom tissue covariates data, to inform calculation of gene mutation rates.
* preload_maf() detects problems in MAF data and informs quality filtering of records.
* check_sample_overlap() offers a simple way of spotting possible sample duplicates.
* Bug fixes.

# cancereffectsizeR 2.0.0
* Supports targeted, whole-exome, and whole-genome sequencing data with proper handling of covered regions.
* Define sample groups and run analyses with different parameters for each group.
* Provides powerful and flexible mutation rate calculation methods; also support for loading in arbitrary rates.
* Supports arbitrary genome/reference information, with tools for building custom reference data sets (refsets).
* Reference data has been separated into its own package (see ces.refset.hg19).
* Supports signature extraction with built-in or custom SNV signature definitions.
* Adds a model of stage/grade-specific variant selection and models of variant- and gene-level epistatic selection.
* Variant annotations are revamped and new functions assist in filtering variants and viewing annotations.
* Collections of variants can be batched into "compound variants" and treated like single variants.
* Variant effect sizes can be visualized with plotting functionality.
* Save, share, and reload entire analyses.
* Numerous performance and usability improvements.
