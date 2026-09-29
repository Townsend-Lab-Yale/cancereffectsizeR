# Attribute cancer effects to mutational signatures

Within patients and across the cohort, calculate mutational source
probabilities and the share of cancer effects attributable to each
source, where sources are biologically-associated mutational signatures.
See [Attribution of Cancer Origins to Endogenous, Exogenous, and
Preventable Mutational
Processes](https://doi.org/10.1093/molbev/msac084) for background and
applications.

## Usage

``` r
mutational_signature_effects(cesa = cesa, effects = NULL, samples = NULL)
```

## Arguments

- cesa:

  CESAnalysis with cancer effects calculated for the variants of
  interest.

- effects:

  A table of cancer effect estimates for a set of variants, as produced
  with
  [`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md).
  Different sets of variants (or parameter choices in the `ces_variant`
  run) will affect output.

- samples:

  Samples for which to calculate mutational sources and effect shares;
  defaults to all samples. Reported averages apply to the samples
  included.

## Value

A nested list containing...

- **mutational_sources** (list):

  - **source_probabilities** (data.table): For each variant in each
    sample, the probability that each signature was the source of the
    variant.

  - **average_by_variant** (data.table): For each distinct variant, the
    source probabilities averaged over all samples with the variant.

  - **average_source_shares** (numeric): For each signature, the average
    proportion of each sample's mutations that are attributable to the
    signature. Calculated by averaging
    `[CESAnalysis]$mutational_signatures$biological_weights` of included
    samples. Compare with average_effect_shares (described below) to
    identify signatures with disproportionate contributions to
    oncogenesis.

- **effect_shares** (list):

  - **by_sample** (data.table): The share of each sample's cancer effect
    (summed across the sample's variants) attributable to each
    signature.

  - **average_effect_shares** (numeric): For each signature, the average
    proportion of each sample's cancer effects that are attributable to
    signature. Compare with average_source_shares, above, to identify
    signatures with disproportionate contributions to oncogenesis. Only
    variants present in `effects` are included in this calculation.
