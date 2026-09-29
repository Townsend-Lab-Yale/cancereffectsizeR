# Plot continuous-covariate cancer effects

Visualize cancer-effect estimates from the linear, logistic, and
generalized sigmoid continuous-selection models. The function can
display all three continuous models together, select the best-fitting
model by AIC, or display any user-specified subset of the three models.
It can return the continuous plot alone or place a discrete-cohort plot
beside the continuous plot.

## Usage

``` r
plot_effects_continuous(
  cesa,
  linear_run_name,
  logistic_run_name,
  sigmoid_run_name,
  cohort_run_names = NULL,
  variants = NULL,
  covariate_col,
  covariate_range = NULL,
  x_range = c(0, 100),
  n_covariate_points = 100,
  output = c("continuous", "cohort_continuous"),
  continuous_models = "best",
  title = NULL,
  x_title = NULL,
  y_title = "Strength of selection",
  cohort_x_title = "Cohort",
  legend.position = c(0.82, 0.58),
  cohort_colors = c("#E71F19", "#F4A016", "#2B6A99"),
  model_colors = c(Linear = "#C2A5CF", Logistic = "#2B6A99", Sigmoid = "#F4A016")
)
```

## Arguments

- cesa:

  A `CESAnalysis` object containing the requested selection runs.

- linear_run_name:

  Run name for point estimates from
  [`ces_variant_linear()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_linear.md).

- logistic_run_name:

  Run name for point estimates from
  [`ces_variant_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_logistic.md).

- sigmoid_run_name:

  Run name for point estimates from
  [`ces_variant_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_sigmoid.md).

- cohort_run_names:

  Optional character vector of run names for discrete cohort-specific
  effect estimates. Required when `output = "cohort_continuous"`. If the
  vector is named, its names are used as the cohort labels and their
  order is preserved; otherwise the run names themselves are used as
  labels. Cohort runs must contain `ci_low_95` and `ci_high_95` for the
  plotted 95% confidence intervals.

- variants:

  Optional character vector of variant IDs or variant names to plot. If
  `NULL`, all variants present in all three continuous-model runs are
  plotted.

- covariate_col:

  Name of the numeric continuous-covariate column in `cesa@samples`.

- covariate_range:

  Optional numeric vector of length two giving the minimum and maximum
  covariate values over which to draw the fitted curves. If `NULL`, the
  observed range of `covariate_col` is used.

- x_range:

  Numeric vector of length two giving the displayed x-axis limits for
  the continuous panel. Defaults to `c(0, 100)`.

- n_covariate_points:

  Number of points used to draw each continuous curve.

- output:

  Either `"continuous"` for the continuous plot alone or
  `"cohort_continuous"` for discrete-cohort and continuous panels side
  by side.

- continuous_models:

  Which continuous models to display. Use `"best"` for the lowest-AIC
  model, `"all"` for all three models, or a character vector containing
  one or more of `"linear"`, `"logistic"`, and `"sigmoid"`.

- title:

  Optional plot title. With multiple variants, a supplied title is used
  as a prefix; otherwise each variant name is used.

- x_title:

  X-axis title for the continuous plot. If `NULL`, `covariate_col` is
  used.

- y_title:

  Y-axis title used for both continuous and cohort plots.

- cohort_x_title:

  X-axis title for the discrete-cohort plot.

- legend.position:

  Position of the continuous-model line legend. The default places the
  legend inside the continuous panel.

- cohort_colors:

  Colors for the discrete cohorts. The default provides three colors;
  supply a vector with at least one color per cohort when more cohorts
  are plotted. A named vector may be used to assign colors by cohort
  label.

- model_colors:

  Named character vector giving colors for `"Linear"`, `"Logistic"`, and
  `"Sigmoid"`.

## Value

For one variant, a ggplot object when `output = "continuous"` or a
patchwork object when `output = "cohort_continuous"`. For multiple
variants, a named list of such objects. The best-fitting model and AIC
values are stored as attributes on each returned plot.

## Details

Point-estimate runs are read directly from `cesa@selection_results`
using the supplied run names.

The three continuous models are

\$\$\gamma(x) = \beta_0 + \beta_1 x\$\$

for the linear model,

\$\$\gamma(x) = \frac{\exp(L)}{1 + \exp\[-k(x-m)\]}\$\$

for the logistic model, and

\$\$ \gamma(x) = \exp(C) + \[\exp(L)-\exp(C)\] \frac{x^s}{x^s+m^s} \$\$

for the generalized sigmoid model.

When `continuous_models = "best"`, AIC is calculated separately for each
variant as \\2k - 2\ell\\, using 2, 3, and 4 fitted parameters for the
linear, logistic, and sigmoid models, respectively.

When `output = "cohort_continuous"`, the cohort and continuous panels
use the same y-axis limits. The shared range is determined from the
cohort point estimates and confidence intervals and from the continuous
curve(s) that are actually displayed.

## See also

[`plot_effects`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/plot_effects.md),
[`ces_variant_linear`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_linear.md),
[`ces_variant_logistic`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_logistic.md),
[`ces_variant_sigmoid`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_sigmoid.md)

Other continuous selection models:
[`ces_variant_linear()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_linear.md),
[`ces_variant_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_logistic.md),
[`ces_variant_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_sigmoid.md)
