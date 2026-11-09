# nonignorable 0.2.0

First release on CRAN.

* `autoplot()` methods for ggplot2 draw the screen, the scans (with their
  confidence sets and contrasts) and the distribution-function scan as ggplot
  objects, with the same content as the `plot()` methods. ggplot2 is only
  suggested: the methods are registered when it is loaded, and the base
  `plot()` methods are unchanged.

# nonignorable 0.1.0

First release, on GitHub and Zenodo.

## Screen

* `shadow_screen()` reports, for each candidate shadow variable, its
  relevance for the outcome among the observed rows and its exclusion
  violation in the selection equation, both as incremental McFadden R². It
  has `print()` and `plot()` methods.
* `numeric = "classes"` codes numeric characteristics as atom-free classes
  (`screen_classes()`) instead of linearly. Linear coding understates the
  exclusion violation when selection cuts hard on a characteristic.
* The screen reports the cases of the rarer class per parameter of the
  outcome and of the response model (`epp_relevance`, `epp_exclusion`) and
  warns below 10, where extra parameters can manufacture relevance or
  violation.
* A ratio whose exclusion chi2/df is below 1 divides by noise. The screen
  flags it in `ratio_readable`; `print()` names such candidates and `plot()`
  draws them hollow without a ratio.
* A numeric outcome with more than two values is screened too: relevance
  comes from a proportional-odds model (`MASS::polr`) on `outcome_classes`
  classes of equal mass, five by default; a factor outcome with more than two
  levels by an ordinal (ordered factor) or a multinomial model (`nnet`). The
  screen measures relevance; it does not certify identification.
* `print()` shows the verdict and its reasons in one column and fits 80
  columns; `plot()` draws rejected candidates grey and labels them.
* `citation("nonignorable")` lists the package and the two papers behind it.

## Scan

* `sharp_scan()` gives the sharp identified set for the rate of a binary
  outcome over a set of rows when the propensity ratio is at most
  `exp(lambda)`: the average of eta(x), or, with `unobserved = TRUE`, the
  rows' own missing outcomes, E(Y | R = 0). Rows are aggregated by `"union"`
  by default; `"common"` is available as an option.
* `sharp_scan_mean()` gives the same for the mean of a non-binary outcome,
  also with `unobserved`,
  taking a grid of conditional quantiles. `quantile_grid()` builds that grid
  from a `quantreg::rq` fit, from any quantile function, or from draws.
* `sharp_scan_ci()` adds a debiased estimate and a confidence set for the
  binary scan, with cross-fitted inputs from `crossfit_inputs()`. The
  standard error defaults to the larger of an empirical and a conditional
  one. It reports the rows with small r-hat and warns when the
  inverse weights blow up; `crossfit_inputs()` warns when the response model
  separates observed rows.
* `sharp_scan_mean_ci()` adds a debiased estimate and a confidence set for
  the mean scan, for the mean over all rows, over a segment chosen by X, or
  of the missing outcomes.
  Its score is Neyman-orthogonal in both regimes of the cap, so no margin
  condition is needed. `crossfit_quantiles()` cross-fits the quantile grids
  and r-hat for any learner.
* `sharp_scan_contrast()` compares two disjoint portfolios chosen by X,
  such as treatment arms. It gives the sharp set for the difference of their
  rates, a confidence set, and the selection strength from which "no
  difference" is no longer ruled out.
* `sharp_scan_mean_contrast()` does the same for the difference of two
  portfolio means.
* `sharp_scan_ci()` and `sharp_scan_mean_ci()` both target all rows by
  default. `target = "unobserved"` targets the missing outcomes themselves,
  E(Y | R = 0), such as the default rate of rejected applications, with the
  sharp set {psi - E(RY)} / pr(R = 0) and its influence function; a
  `segment` chosen by X restricts it. Both contrasts take
  `among = "unobserved"` to compare the missing outcomes of two segments.
* `sharp_scan_cdf()` gives the sharp set for the distribution function
  P(Y <= t), or P(Y <= t | R = 0), at each value of a grid, with a debiased
  estimate and a confidence set pointwise in t; `quantile()` inverts the
  bounds into bounds for quantiles, widened outwards on the grid. The sets
  are sharp for each t, not as a band, and no uniform band is offered. The
  output reports the share of rows near the kink of each bound, where the
  normal approximation needs the margin condition. `crossfit_cdf()` fits
  P(Y <= t | R = 1, X) as one classification per t for any learner.
* Debiased bounds and confidence sets are kept monotone in lambda.
* The bounds without the cap pi <= 1 are called the outer band throughout
  (`naive_*` columns). Scans report the share of the rows entering the set where the cap
  binds and, for rates, the lambda from which it holds a bound flat; the
  plot notes such a bound in its legend, which it places clear of the
  curves, and leaves the flip labels a free strip below them.
  `sharp_scan_mean_ci()` counts rows with r-hat clipped below and above
  apart.
* `simulate_shadow()`, `obs_quantiles()` and `scan_inputs()` provide data and
  inputs for examples and tests.
