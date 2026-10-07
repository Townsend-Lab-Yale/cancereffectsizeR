# Plot pairwise epistasis

Visualize pairwise epistatic scaled selection coefficients, as
calculated by
[`ces_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_epistasis.md)
or
[`ces_gene_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_gene_epistasis.md).
For each variant, the isolated site effect–the overall scaled selection
at the site without regard for the mutation status of the other site–is
also depicted.

## Usage

``` r
plot_epistasis(
  epistatic_effects,
  pairs_per_row = 8,
  x_title = "Site pairs",
  variant_label_size = 6.5,
  dodge_labels = "auto",
  alternating_colors = c("#7cb4de", "#7DD3AF"),
  include_schematic = TRUE,
  schematic_label_size = 2.5,
  significance_levels = c(0.05, 0.01, 0.001),
  significance_cols = list(A_change = "p_A_change", B_change = "p_B_change", model =
    "p_epistasis"),
  inference_floor = 0.001
)
```

## Arguments

- epistatic_effects:

  Epistatic effects inference table, as produced by
  [`ces_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_epistasis.md)
  or
  [`ces_gene_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_gene_epistasis.md).

- pairs_per_row:

  How many epistatic pairs to show in each plot row. The provided value
  is incremented if need to prevent the legend schematic from being
  isolated on its own row.

- x_title:

  X-axis label. Set NULL for no label.

- variant_label_size:

  Text size for the variant labels.

- dodge_labels:

  TRUE/FALSE on using n.dodge (height staggering) on variant labels.
  Defaults to "auto"; you can try setting manually if labels are not
  looking good.

- alternating_colors:

  Colors, provided as character vector, to use on epistatic effect
  arrows. It's recommended to supply one or two colors, but more will
  work.

- include_schematic:

  TRUE/FALSE on whether to include the schematic that shows how to
  interpret the plot. If you need to put the schematic somewhere
  special, set to FALSE and then get your own copy of it with
  [`epistasis_plot_schematic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/epistasis_plot_schematic.md).

- schematic_label_size:

  Text size of labels in the schematic.

- significance_levels:

  A vector of 1-3 distinct numeric values on (0, 1) in descending order
  to use for significance annotations.

- significance_cols:

  A named list of column names that give significance values for nonzero
  change in selection for each pair of sites A and B, and for the
  performance of the epistatic selection model over a model that ignores
  epistatic interactions. List elements must be named A_change (default
  "p_A_change") B_change ("p_B_change"), and model ("p_epistasis"). The
  purpose of this option is to support the use of different columns when
  multiple testing correction is performed.

- inference_floor:

  Numeric value of the optimization floor used in epistatic effect
  inference. Typically, should be left at the default value, which
  matches cancereffectsizeR epistatic inference defaults. For plot
  legibility, there will be a dashed horizontal line in the output plot,
  higher than `inference_floor` and lower than any non-minimized
  parameter inference, which indicates that all arrows pointing below
  the line have estimates at the optimization floor.

## Details

Variant pairs for which the epistatic selection model is not
significantly better than a no-epistasis model (in which
`ces_A0 = ces_A_on_B` and `ces_B0 = ces_B_on_A`) are depicted with faded
arrows; the threshold is `p_epistasis < .05` (or, instead of
p_epistasis, whatever corrected significance column is specified with
`significance_cols`).
