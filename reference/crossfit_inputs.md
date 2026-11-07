# Cross-fitted scan inputs for any learner

Splits the rows into folds and predicts \\\hat\eta\_{\mathrm{obs}}\\ and
\\\hat r\\ for each fold from models fitted on the other folds, as
[`sharp_scan_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_ci.md)
requires. The learner is yours: `fit_eta` and `fit_r` are functions that
fit on a training frame and return probabilities for new rows. `fit_eta`
receives only the observed rows of the training folds.

## Usage

``` r
crossfit_inputs(
  data,
  response,
  fit_eta,
  fit_r,
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

- fit_eta:

  `function(train, newdata)` returning \\P(Y = 1 \mid X, R = 1)\\ for
  `newdata`; `train` holds observed rows.

- fit_r:

  `function(train, newdata)` returning \\P(R = 1 \mid X)\\.

- folds:

  Number of folds. Default 5.

- seed:

  Optional seed for the fold assignment.

- clip:

  Predictions are clipped to `[clip, 1 - clip]`. Default 1e-6.

## Value

A list with `eta_obs`, `r` (one per row of `data`) and `fold`.

## See also

[`sharp_scan_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_ci.md)
