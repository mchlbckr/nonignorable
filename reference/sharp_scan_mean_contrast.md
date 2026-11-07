# Debiased scan for the difference of two portfolio means

`sharp_scan_mean_contrast()` is
[`sharp_scan_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_contrast.md)
for an outcome that is not binary: it compares the means of two disjoint
portfolios chosen by \\X\\ under the assumption of
[`sharp_scan_mean()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean.md),
with each endpoint estimated by
[`sharp_scan_mean_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean_ci.md).
The sharp set for \\A - B\\ is \\\[\mathrm{lower}\_A -
\mathrm{upper}\_B,\\ \mathrm{upper}\_A - \mathrm{lower}\_B\]\\; because
the targets are disjoint, the standard error of each endpoint is the
root of the sum of the squared standard errors of the arms. By default
each target is a whole segment of rows, observed and unobserved alike;
`among = "unobserved"` compares the means of the missing outcomes of the
two segments, \\E(Y \mid R = 0)\\ in each.

## Usage

``` r
sharp_scan_mean_contrast(
  q,
  r,
  response,
  outcome,
  target_a,
  target_b,
  lambda = seq(0, 2, length.out = 81),
  weights = NULL,
  threshold = 0,
  level = 0.9,
  trim = 0,
  r_bounds = c(0.01, 0.99),
  labels = c("A", "B"),
  among = c("all", "unobserved")
)
```

## Arguments

- q:

  An \\n \times K\\ matrix of conditional quantiles of \\Y\\ given \\R =
  1, X = x_i\\ for **every** row, fitted out of fold, at the levels
  \\(k - 1/2)/K\\.

- r:

  Cross-fitted \\\hat r(x_i)\\ for every row, in (0, 1\].

- response:

  Response indicator for every row, 0/1 or logical.

- outcome:

  Outcome; ignored (and may be `NA`) where `response` is 0.

- target_a, target_b:

  Logical vectors, one per row, selecting the two portfolios by \\X\\.
  They must be disjoint and non-empty.

- lambda:

  Grid of \\\lambda\\ values; must include 0.

- weights:

  Optional non-negative row weights, a function of \\X\\.

- threshold:

  The difference at which the decision changes. Default 0.

- level:

  Confidence level. Default 0.9.

- trim:

  Rows with \\\hat r \<\\ `trim` leave the target. Default 0.

- r_bounds:

  \\\hat r\\ is clipped into this range before use. Default
  `c(0.01, 0.99)`. The number of clipped target rows is reported, below
  and above apart: clipping at the top only stabilises the cap \\\pi \le
  1\\ for rows that are almost always observed, clipping at the bottom
  limits inverse weights.

- labels:

  Names of the two portfolios, for printing and plotting.

- among:

  `"all"` (default) or `"unobserved"`, as in
  [`sharp_scan_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_contrast.md).

## Value

An object of class
`c("sharp_scan_contrast", "sharp_scan_ci", "sharp_scan")` as from
[`sharp_scan_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_contrast.md),
with `of_means = TRUE`.

## See also

[`sharp_scan_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_contrast.md),
[`sharp_scan_mean_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean_ci.md),
[`crossfit_quantiles()`](https://buecker.ms/nonignorable/reference/crossfit_quantiles.md)

## Examples

``` r
d <- simulate_shadow(3000, tilt = -1, outcome = "continuous", seed = 3)
q <- obs_quantiles(d, K = 50)   # the true quantile grids of the simulation
k <- sharp_scan_mean_contrast(q, d$obs_rate, d$observed, d$y,
                               target_a = d$x1 > 0, target_b = d$x1 <= 0,
                               lambda = seq(0, 2, length.out = 21),
                               labels = c("x1 > 0", "x1 <= 0"))
k
#> <sharp_scan_contrast> difference of means x1 > 0 - x1 <= 0, debiased, with 90% confidence set
#>   rows: 3000   x1 > 0: 1507 rows   x1 <= 0: 1493 rows
#>   clipped: 5 rows entering the sets with r-hat below 0.01, 0 above 0.99
#>   r-hat < 0.01: 5 target rows, 0 observed rows; smallest r-hat observed 0.025
#>   lambda = 0: 1.21 [1.05, 1.37]
#>   lambda = 2 (propensity ratio 7.39): set [-0.41, 2.74], confidence set [-0.599, 2.92]
#>   difference 0 reached by
#>     estimated set:        lambda = 1.45 (propensity ratio 4.26) 
#>     confidence set:       lambda = 1.24 (propensity ratio 3.44) 
#>     outer band (plug-in): lambda = 1.64 (propensity ratio 5.15) 
plot(k)
```
