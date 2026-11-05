# shadowscan 0.1.0

First release.

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
* `citation("shadowscan")` lists the package and the two papers behind it.

## Scan

* `shadow_scan()` gives the sharp identified set for the rate of a binary
  outcome among the unobserved rows when the propensity ratio is at most
  `exp(lambda)`. Rows are aggregated by `"union"` by default; `"common"` is
  available as an option.
* `shadow_scan_mean()` gives the same for the mean of a non-binary outcome,
  taking a grid of conditional quantiles. `quantile_grid()` builds that grid
  from a `quantreg::rq` fit, from any quantile function, or from draws.
* `shadow_scan_ci()` adds a debiased estimate and a confidence set for the
  binary scan, with cross-fitted inputs from `crossfit_inputs()`. The
  standard error defaults to the larger of an empirical and a conditional
  one. It reports the rows with small r-hat and warns when the
  inverse weights blow up; `crossfit_inputs()` warns when the response model
  separates observed rows.
* `shadow_scan_mean_ci()` adds a debiased estimate and a confidence set for
  the mean scan, for the mean over all rows or over a segment chosen by X.
  Its score is Neyman-orthogonal in both regimes of the cap, so no margin
  condition is needed. `crossfit_quantiles()` cross-fits the quantile grids
  and r-hat for any learner.
* `shadow_scan_contrast()` compares two disjoint portfolios chosen by X,
  such as treatment arms. It gives the sharp set for the difference of their
  rates, a confidence set, and the selection strength from which "no
  difference" is no longer ruled out.
* `simulate_shadow()`, `obs_quantiles()` and `scan_inputs()` provide data and
  inputs for examples and tests.
