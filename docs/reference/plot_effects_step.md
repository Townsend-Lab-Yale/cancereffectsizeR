# Plot step-specific cancer effects

Visualize cancer effects estimated separately across ordered progression
stages, as produced by
[`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md).
Works for any number of stages, not just two: one panel is drawn per
variant (or gene, or other grouping column), with stage on the x-axis.

## Usage

``` r
plot_effects_step(
  effects,
  group_by = "variant_name",
  stage_order = NULL,
  stage_labels = NULL,
  show_ci = TRUE,
  lrt = NULL,
  ncol = NULL,
  color_by = "stage",
  viridis_option = "cividis",
  title = "",
  x_title = NULL,
  y_title = NULL,
  legend.position = "right"
)
```

## Arguments

- effects:

  Cancer effects table produced by
  [`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md).
  You can combine multiple such tables via rbind() to plot multiple runs
  together.

- group_by:

  Facet the plot by this column: "variant_name" (default), "gene", or
  another column present in `effects`.

- stage_order:

  Optional character vector giving stage display order (matching the
  `si_<stage>` column suffixes in `effects`). Defaults to the order
  those columns appear in the table, which matches the stage order used
  in the original
  [`ces_variant_step()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_step.md)
  call.

- stage_labels:

  Optional named character vector to relabel stages for display, e.g.
  `c(Pre = "Normal tissue", Pri = "Primary tumor")`. Names should match
  the stage names (i.e., `si_<stage>` suffixes, or `stage_order` if
  supplied).

- show_ci:

  TRUE/FALSE to depict confidence intervals (error bars) in the plot
  (default TRUE; ignored with a message if `effects` has no
  ci_low/ci_high columns).

- lrt:

  Optional output of
  [`step_selection_LRT()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/step_selection_LRT.md).
  When supplied, each facet is annotated with a significance code
  (`*`/`**`/`***`; `ns` if not significant) based on p_value. Requires
  `group_by` to identify a single row of `lrt` per facet (true when
  `group_by` is "variant_name", or "gene" when each gene has exactly one
  variant/compound in `effects`).

- ncol:

  Number of facet columns (default: ggplot chooses).

- color_by:

  Set to "stage" (default) to color points by stage (viridis discrete
  scale), a single R color to use throughout, or NULL to disable
  coloring.

- viridis_option:

  Viridis color map option, used when `color_by = "stage"`.

- title:

  Main plot title (default none).

- x_title, y_title:

  Axis titles.

- legend.position:

  Passed to ggplot's legend.position (none, left, right, top, bottom).

## Value

A ggplot
