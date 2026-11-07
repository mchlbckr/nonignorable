# Extract the two scan inputs from fitted models

A convenience, not an architectural commitment: calls
[`predict()`](https://rdrr.io/r/stats/predict.html) with
`type = "response"` on each model at `newdata`. For a classification
forest that returns class probabilities (e.g. `ranger` with
`probability = TRUE`) the column for class `"1"` (or the last column) is
used. Anything else, compute the two vectors yourself.

## Usage

``` r
scan_inputs(eta_model, r_model, newdata)
```

## Arguments

- eta_model:

  Model of \\Y\\ fitted on the observed rows.

- r_model:

  Model of \\R\\ fitted on all rows.

- newdata:

  Rows at which the set is wanted.

## Value

A list with numeric vectors `eta_obs` and `r`.

## Examples

``` r
d <- simulate_shadow(3000, seed = 2)
f <- ~ x1 + x2 + z_shadow + z_policy + z_score
m_eta <- glm(update(f, y ~ .), binomial, data = d[d$observed == 1, ])
m_r <- glm(update(f, observed ~ .), binomial, data = d)
inp <- scan_inputs(m_eta, m_r, d[d$observed == 0, ])
sharp_scan(inp$eta_obs, inp$r)
#> <sharp_scan> sharp identified set for a rate under a propensity ratio <= exp(lambda)
#>   rows: 1347   aggregation: union 
#>   lambda = 0 (model on observed rows): 0.217 
#>   lambda = 2 (propensity ratio 7.39): sharp [0.0701, 0.533], outer [0.0511, 0.544]
#>   cap pi <= 1 at lambda = 2 binds in 71% of the rows entering the set
#>     for the lower bound, 14% for the upper
```
