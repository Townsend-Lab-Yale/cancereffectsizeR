# Get epistatic effect schematic

Get a copy of the explanatory schematic that appears in epistatic effect
plots (see
[`plot_epistasis()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_epistasis.md)).
May be useful for putting the schematic in custom locations when
assembling complex figures.

## Usage

``` r
epistasis_plot_schematic(
  title = "Types of effects",
  schematic_label_size = 3,
  with_border = TRUE
)
```

## Arguments

- title:

  Schematic title text.

- schematic_label_size:

  Text size of labels in the schematic (title gets size + 1).

- with_border:

  TRUE/FALSE on the appearance of a thin visible border around the
  schematic.

## Value

The schematic as a ggplot.
