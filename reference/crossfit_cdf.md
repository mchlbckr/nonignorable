# Cross-fitted inputs for the distribution-function scan

Splits the rows into folds and predicts, for each fold from models
fitted on the other folds, \\\hat\eta_t(x) = P(Y \le t \mid R = 1, X =
x)\\ at every value of `t`, by a classification of \\Z = 1(Y \le t)\\
among the observed rows, and \\\hat r\\, as
[`sharp_scan_cdf()`](https://buecker.ms/nonignorable/reference/sharp_scan_cdf.md)
requires. The learner is yours: `fit_eta(train, newdata)` receives the
observed rows of the training folds with a 0/1 column `z` added and
returns probabilities for `newdata`; it is called once per fold and
value of `t`. Where all training rows share one value of `z` (a `t`
below or above the observed outcomes), the prediction is that value and
`fit_eta` is not called.

## Usage

``` r
crossfit_cdf(
  data,
  response,
  outcome,
  t,
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

  Name of the 0/1 response column.

- outcome:

  Name of the outcome column.

- t:

  Increasing grid of values.

- fit_eta, fit_r:

  Functions `(train, newdata)` returning probabilities.

- folds:

  Number of folds. Default 5.

- seed:

  Optional seed for the fold assignment.

- clip:

  Bound for \\\hat r\\ and \\\hat\eta_t\\, as in
  [`crossfit_inputs()`](https://buecker.ms/nonignorable/reference/crossfit_inputs.md).

## Value

A list with `eta_t` (an \\n \times T\\ matrix), `r`, `fold` and `t`.

## See also

[`sharp_scan_cdf()`](https://buecker.ms/nonignorable/reference/sharp_scan_cdf.md),
[`crossfit_inputs()`](https://buecker.ms/nonignorable/reference/crossfit_inputs.md)
