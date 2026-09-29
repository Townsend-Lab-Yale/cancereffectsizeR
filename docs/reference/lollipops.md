# Plot selection intensities

Lollipop plots can be used to compare highly selected variants within or
among groups of samples.

## Usage

``` r
lollipops(
  si_list,
  title = "My SIs",
  ylab = "selection intensity",
  max_sites = 50,
  label_size = 3,
  merge_dist = 0.04,
  si_col = "auto",
  group_names = "auto"
)
```

## Arguments

- si_list:

  Selection output table or named list of selection output tables. Each
  table should have columns called "variant_name", "variant_type", and
  "gene", plus a selection intensity column (or columns, and then each
  column will have its values plotted on a separate lollipop.)

- title:

  Plot title to pass to ggplot.

- ylab:

  Y-axis label to pass to ggplot.

- max_sites:

  Maximum number of variant sites to include per lollipop; if you try to
  include too many, you may have a challenge getting it to look good.

- label_size:

  Text size for labels, either 1-length or same length as si_list.

- merge_dist:

  How close points must be to be eligible to have their labels combined
  (.04 = 4 percent of plot space). Either 1-length or same length as
  si_list. Try tweaking if labels are not looking good. Set to zero to
  prevent any labels being combined.

- si_col:

  Vector giving names of selection intensity column in each SI table,
  either 1-length or same length as si_list.

- group_names:

  Names to use for labeling each group of selection intensities.

## Value

ggplot object with lollipops

## Details

To keep the lollipops readable, no more than `max_sites` variants are
shown for any group. At the same time, all SIs that are within the
limits of the output plot are depicted (which means that specification
of `max_sites` creates an implicit SI floor). If the labels don't look
good, you may have to reduce max_sites or break up the plot into pieces.
(If you find a better way to get lots of variant labels to display
nicely, we would love to hear about it!)
