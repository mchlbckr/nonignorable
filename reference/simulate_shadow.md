# Simulate a nonignorably selected outcome with known truth

Draws covariates, an outcome \\Y\\ and an observation indicator \\R\\
whose propensity depends on \\Y\\, with the ratio of response
propensities between outcome values bounded by
\\e^{\|\mathrm{tilt}\|}\\. Because the mechanism is known, the target,
the observed-data law and the conditional observation rate \\r(x) = P(R
= 1 \mid X = x)\\ are returned exactly. The simulator exists so that
examples and tests depend on no outside data.

## Usage

``` r
simulate_shadow(
  n = 5000,
  tilt = -1,
  outcome = c("binary", "continuous"),
  seed = 1L
)

obs_quantiles(sim, K = 200L)
```

## Arguments

- n:

  Number of rows.

- tilt:

  The log of the largest ratio \\b\\ of response propensities. Negative
  values mean that a high outcome (a default, \\Y = 1\\, or a large
  loss) lowers the chance of being observed, as when such cases are more
  likely to have been rejected. `0` is ignorable selection.

- outcome:

  `"binary"` (default) or `"continuous"`.

- seed:

  Seed for the random number generator. The caller's RNG state is
  restored on exit.

- sim:

  A data frame from `simulate_shadow(outcome = "continuous")`.

- K:

  Number of quantile levels, \\(k - 1/2)/K\\.

## Value

A data frame with columns `x1`, `x2`, `z_shadow`, `z_policy`, `z_score`,
`y` (the outcome, `NA` where unobserved), `observed` (0/1) and the exact
truth columns described above. Attributes `"tilt"` and `"outcome"`.

`obs_quantiles()`: the exact \\n \times K\\ matrix of conditional
quantiles of \\Y\\ given \\R = 1, X = x_i\\.

## Details

The covariates play fixed roles, so that
[`shadow_screen()`](https://buecker.ms/nonignorable/reference/shadow_screen.md)
has something to find:

- `x1`, `x2`: ordinary covariates entering both outcome and selection.

- `z_shadow`: enters the outcome only, a valid shadow variable.

- `z_policy`: enters the selection only, high exclusion violation and no
  relevance.

- `z_score`: enters both strongly, like a bureau score that is also used
  in the acceptance decision.

Linear index of the outcome: \\\mu = 0.8 x_1 - 0.5 x_2 + 0.6
z\_{\mathrm{shadow}} + 1.0 z\_{\mathrm{score}}\\; of selection: \\g =
0.8 - 0.8 x_1 + 0.4 x_2 + 1.0 z\_{\mathrm{policy}} - 1.0
z\_{\mathrm{score}}\\. All truth columns condition on all five
covariates.

**Binary** (`outcome = "binary"`): \\\eta = \mathrm{expit}(\mu - 1.5)\\,
\\\pi(0,x) = \mathrm{expit}(g)\min(1, e^{-b})\\ and \\\pi(1,x) =
\pi(0,x)e^{b}\\ with \\b\\ = `tilt`. Truth: `eta`, `eta_obs`,
`obs_rate`.

**Continuous** (`outcome = "continuous"`): \\Y = \mu + \varepsilon\\,
\\\varepsilon \sim N(0,1)\\, and \\\pi(y,x) = c(x)\exp\\b\\
\Phi(\varepsilon)\\\\ with \\c(x) = \mathrm{expit}(g)\min(1, e^{-b})\\.
Since \\\Phi(\varepsilon)\\ is uniform, the ratio of propensities is at
most \\e^{\|b\|}\\, \\r(x) = c(x)(e^{b} - 1)/b\\, and the quantile
function of \\Y\\ given \\R = 1, X = x\\ is \\\mu(x) +
\Phi^{-1}\[\log\\1 + p(e^{b} - 1)\\/b\]\\, all in closed form. Truth:
`mu` \\= E(Y \mid x)\\, `mean_obs` \\= E(Y \mid R = 1, x)\\, `obs_rate`;
`obs_quantiles()` returns the quantile grid that
[`sharp_scan_mean()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean.md)
takes.

## Examples

``` r
d <- simulate_shadow(2000, tilt = -1, seed = 1)
head(d)
#>           x1          x2   z_shadow    z_policy    z_score  y observed
#> 1 -0.6264538 -0.88614959 -1.1346302 -0.61882708  0.2637034  0        1
#> 2  0.1836433 -1.92225490  0.7645571 -1.10942196 -0.8294518 NA        0
#> 3 -0.8356286  1.61970074  0.5707101 -2.17033523 -1.4616348 NA        0
#> 4  1.5952808  0.51926990 -1.3516939 -0.03130307  1.6839902 NA        0
#> 5  0.3295078 -0.05584993 -2.0298855 -0.26039848 -1.5443243 NA        0
#> 6 -0.8204684  0.69641761  0.5904787  0.53443047 -0.1908871  0        1
#>          eta     eta_obs   obs_rate
#> 1 0.12183541 0.048560624 0.47629521
#> 2 0.31806567 0.146455571 0.32143524
#> 3 0.01634164 0.006074504 0.79510116
#> 4 0.59617762 0.351959635 0.07534637
#> 5 0.01851066 0.006890315 0.84786452
#> 6 0.08776912 0.034185063 0.87019366
dc <- simulate_shadow(2000, tilt = -1, outcome = "continuous", seed = 1)
head(dc)
#>           x1          x2   z_shadow    z_policy    z_score         y observed
#> 1 -0.6264538 -0.88614959 -1.1346302 -0.61882708  0.2637034        NA        0
#> 2  0.1836433 -1.92225490  0.7645571 -1.10942196 -0.8294518        NA        0
#> 3 -0.8356286  1.61970074  0.5707101 -2.17033523 -1.4616348        NA        0
#> 4  1.5952808  0.51926990 -1.3516939 -0.03130307  1.6839902        NA        0
#> 5  0.3295078 -0.05584993 -2.0298855 -0.26039848 -1.5443243        NA        0
#> 6 -0.8204684  0.69641761  0.5904787  0.53443047 -0.1908871 -1.366172        1
#>           mu   mean_obs   obs_rate
#> 1 -0.4751630 -0.7532270 0.32619803
#> 2  0.7373245  0.4592605 0.25431793
#> 3 -2.5975620 -2.8756260 0.50784578
#> 4  1.8895635  1.6114995 0.07643177
#> 5 -2.4707244 -2.7487884 0.54229801
#> 6 -0.8411834 -1.1192475 0.58237802
dim(obs_quantiles(dc, K = 50))
#> [1] 2000   50
```
