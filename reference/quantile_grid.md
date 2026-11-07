# Quantile grid for the mean scan

`quantile_grid()` produces the \\n \times K\\ matrix of conditional
quantiles of \\Y\\ given \\R = 1, X = x_i\\ at the levels \\(k -
1/2)/K\\ that
[`sharp_scan_mean()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean.md)
takes as `q`. It is a convenience like
[`scan_inputs()`](https://buecker.ms/nonignorable/reference/scan_inputs.md),
not a commitment to a learner, and accepts three kinds of `fit`:

## Usage

``` r
quantile_grid(fit, newdata = NULL, K = 200L)
```

## Arguments

- fit:

  An `rq`/`rqs` object, a `function(tau, newdata)`, or a numeric matrix
  of draws.

- newdata:

  Rows at which the quantiles are wanted (for `rq` and function fits).

- K:

  Number of levels. Default 200.

## Value

A numeric matrix, `nrow(newdata)` (or `nrow(fit)` for draws) by `K`,
rows sorted, with attributes `"levels"`, `"crossing"` and `"shift"`.

## Details

- a quantile regression from quantreg (`rq` or `rqs`), fitted on the
  observed rows. It is refitted at the `K` levels with
  [`update()`](https://rdrr.io/r/stats/update.html), so the data it was
  fitted on must still be reachable, and predicted at `newdata`;

- a function `function(tau, newdata)` returning, for a vector of levels
  `tau`, a matrix with one row per row of `newdata` and one column per
  level – for any other conditional quantile model;

- a numeric matrix of draws from the conditional distribution, one row
  per row of interest and any number of draws per row – for a
  distributional or generative model. `newdata` is then not used.

Quantile curves estimated level by level can cross. Each row is sorted
(rearranged), which yields a proper quantile function and never moves it
further from the truth (Chernozhukov, Fernández-Val and Galichon, 2010).
With many levels nearly every row crosses somewhere by sampling noise
alone, so the share of rows that needed it (attribute `"crossing"`) says
little. Attribute `"shift"` says how much was moved: the median over
rows of the largest change by sorting, relative to the row's range. A
few percent is noise; a large shift says the quantile model is poor.

## References

Chernozhukov, V., Fernández-Val, I. and Galichon, A. (2010). Quantile
and probability curves without crossing. *Econometrica*, 78(3),
1093–1125. [doi:10.3982/ECTA7880](https://doi.org/10.3982/ECTA7880)

## See also

[`sharp_scan_mean()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean.md)

## Examples

``` r
dc <- simulate_shadow(3000, tilt = -1, outcome = "continuous", seed = 1)
obs <- dc[dc$observed == 1, ]
rej <- dc[dc$observed == 0, ]
# any conditional quantile model, here a location model fitted by lm()
m <- lm(y ~ x1 + x2 + z_shadow + z_policy + z_score, data = obs)
res <- residuals(m)
q <- quantile_grid(function(tau, newdata)
  outer(predict(m, newdata), quantile(res, tau, names = FALSE), "+"),
  rej, K = 100)
sm <- sharp_scan_mean(q, rej$obs_rate, threshold = 0.5)
sm
#> <sharp_scan_mean> sharp identified set for a mean under a propensity ratio <= exp(lambda)
#>   rows: 1803   aggregation: union 
#>   lambda = 0 (model on observed rows): -0.0117 
#>   lambda = 2 (propensity ratio 7.39): sharp [-0.718, 0.762], outer [-0.741, 0.775]
#>   threshold 0.5 reached by
#>     sharp set:   lambda = 1.31 (propensity ratio 3.71) 
#>     outer band:  lambda = 1.31 (propensity ratio 3.71) 
if (requireNamespace("quantreg", quietly = TRUE)) {
  fq <- quantreg::rq(y ~ x1 + x2 + z_shadow + z_policy + z_score,
                     tau = 0.5, data = obs)
  q2 <- quantile_grid(fq, rej, K = 50)
  attr(q2, "crossing")
}
#> [1] 0.452579
```
