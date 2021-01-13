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
