# sswm_sequential_lik

As in sswm_lik, selection intensities are calculated at variant sites
under a "strong selection, weak mutation" assumption. In this version,
each sample is assigned to one of an ordered set of disease progression
states, and selection is assumed to vary across states. For example, in
a two-state local/metastatic model, each variant has two independent
selection intensities. Metastatic samples could have acquired the
variant while in their current state or at some earlier time, while the
local state selection intensity applied.

## Usage

``` r
sswm_sequential_lik(rates_tumors_with, rates_tumors_without, sample_index)
```

## Arguments

- rates_tumors_with:

  named vector of site-specific mutation rates for all tumors with
  variant

- rates_tumors_without:

  named vector of site-specific mutation rates for all eligible tumors
  without variant

- sample_index:

  data.table with columns Unique_Patient_Identifier, group_name,
  group_index

## Details

All arguments to this likelihood function factory are automatically
supplied by
[`ces_variant()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant.md).
