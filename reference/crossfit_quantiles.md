# Cross-fitted quantile grids for the debiased mean scan

Splits the rows into folds and predicts, for each fold from models
fitted on the other folds, a grid of conditional quantiles of \\Y\\
given \\R = 1\\ and \\X\\, and \\\hat r\\, as
[`sharp_scan_mean_ci()`](https://mchlbckr.github.io/nonignorable/reference/sharp_scan_mean_ci.md)
requires. The learner is yours. `fit_q` receives only the observed rows
of the training folds and the levels `tau` = \\(k - 1/2)/K\\.

## Usage

``` r
crossfit_quantiles(
  data,
  response,
  fit_q,
  fit_r,
  K = 200L,
  folds = 5L,
  seed = NULL,
  clip = 1e-06
)
```

## Arguments

- data:

  A data frame with all rows, observed and unobserved.

- response:

  Name of the response indicator column (0/1 or logical).

- fit_q:

  `function(train, newdata, tau)` returning a matrix with one row per
  row of `newdata` and one column per level in `tau`; `train` holds
  observed rows.

- fit_r:

  `function(train, newdata)` returning \\P(R = 1 \mid X)\\.

- K:

  Number of quantile levels. Default 200.

- folds:

  Number of folds. Default 5.

- seed:

  Optional seed for the fold assignment.

- clip:

  \\\hat r\\ is clipped from below at `clip`. Default 1e-6.

## Value

A list with `q` (an \\n \times K\\ matrix), `r` and `fold`, one per row
of `data`.

## Details

A smooth location and scale model is a good choice: fit \\\mu(x)\\ and
\\\sigma(x)\\ (for instance with
`mgcv::gam(..., family = mgcv::gaulss())`), and return \\\hat\mu(x) +
\hat\sigma(x) \hat F^{-1}(\tau)\\, with \\\hat F^{-1}\\ the quantiles of
the standardised training residuals.
[`quantile_grid()`](https://mchlbckr.github.io/nonignorable/reference/quantile_grid.md)
can be used inside `fit_q`.

## See also

[`sharp_scan_mean_ci()`](https://mchlbckr.github.io/nonignorable/reference/sharp_scan_mean_ci.md),
[`crossfit_inputs()`](https://mchlbckr.github.io/nonignorable/reference/crossfit_inputs.md)
