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
