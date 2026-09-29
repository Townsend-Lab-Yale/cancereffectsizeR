# Calculate uniroot CIs on selection intensities

Given a model fit, calculate univariate confidence intervals for each
parameter. Returns a list of low/high bounds.

## Usage

``` r
univariate_si_conf_ints(fit, lik_fn, min_si, max_si, conf)
```

## Arguments

- fit:

  From bbmle

- lik_fn:

  likelihood function

- min_si:

  lower limit on SI/CI

- max_si:

  upper limit on SI/CI

- conf:

  e.g., .95 -\> 95% CIs
