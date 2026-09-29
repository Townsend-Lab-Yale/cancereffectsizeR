# Plot mutational source and effect attributions

Compare the extent to which mutational signatures contribute mutations
(mutational source share) to the degree to which they contribute
high-effect mutations (cancer effect share).

## Usage

``` r
plot_signature_effects(
  mutational_effects = NULL,
  signature_groupings = "auto",
  viridis_option = NULL,
  num_sig_groups = 7,
  other_color = "white"
)
```

## Arguments

- mutational_effects:

  Output from
  [mutational_signature_effects](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/mutational_signature_effects.md)().
  To compare groups of samples, supply a named list with each element
  corresponding to output from a separate run of
  mutational_signature_effects().

- signature_groupings:

  A data.table of signature names and descriptions; signatures with
  identical descriptions are grouped together. Only signatures present
  in the data get displayed. Setting to "auto" (the default) uses the
  table returned by
  [`cosmic_signature_info()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/cosmic_signature_info.md),
  which only makes sense when using COSMIC signatures. A custom table
  should have columns "name", "short_name", and "description".
  Additional options:

  - To force a signature group to appear in the plot even if it has a
    low effect share, add a column called "prioritize" and set to TRUE
    where desired.

  - To make a signature appear in its own group, make its description
    unique.

  - Add a "color" column to manually specify colors for each group.

  Alternatively, setting `signature_groupings = "cannataro"` applies the
  same signature groupings and color palette as [Cannataro et al.
  2022](https://academic.oup.com/mbe/article/39/5/msac084/6570859). You
  can use Cannataro signature groupings with a different color palette
  by specifying `viridis_option`.

- viridis_option:

  A viridis color mapping, specified with a single letter ('A' to 'H').
  By default, map 'G' (mako) unless using Cannataro signature groupings.

- num_sig_groups:

  How many groups of signatures to display. Groups are ordered by their
  highest effect shares, and the rest get lumped into an "other
  signatures" group.

- other_color:

  Color to use for signatures in the "other" group, supplied as a scalar
  character.
