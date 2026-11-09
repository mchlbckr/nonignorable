# Package index

## Screen

Is a candidate fit to serve as a shadow variable?

- [`shadow_screen()`](https://buecker.ms/nonignorable/reference/shadow_screen.md)
  : Screen candidate shadow variables
- [`plot(`*`<shadow_screen>`*`)`](https://buecker.ms/nonignorable/reference/plot.shadow_screen.md)
  : Relevance against exclusion violation
- [`autoplot(`*`<shadow_screen>`*`)`](https://buecker.ms/nonignorable/reference/autoplot.shadow_screen.md)
  : Relevance against exclusion violation, with ggplot2
- [`screen_classes()`](https://buecker.ms/nonignorable/reference/screen_classes.md)
  : Classes for a numeric characteristic, cut only between distinct
  values

## Scan of a rate

Sharp identified set for the rate of a binary outcome under a bound on
the ratio of response propensities, with known inputs.

- [`sharp_scan()`](https://buecker.ms/nonignorable/reference/sharp_scan.md)
  [`as.data.frame(`*`<sharp_scan>`*`)`](https://buecker.ms/nonignorable/reference/sharp_scan.md)
  : Sharp identified set for a binary outcome under a bounded ratio of
  response propensities
- [`plot(`*`<sharp_scan>`*`)`](https://buecker.ms/nonignorable/reference/plot.sharp_scan.md)
  : Sensitivity curve of a scan
- [`autoplot(`*`<sharp_scan>`*`)`](https://buecker.ms/nonignorable/reference/autoplot.sharp_scan.md)
  : Sensitivity curve of a scan, with ggplot2

## Scan of a mean

The same for the mean of any outcome, from a grid of conditional
quantiles.

- [`sharp_scan_mean()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean.md)
  : Sharp identified set for a mean under a bounded ratio of response
  propensities
- [`quantile_grid()`](https://buecker.ms/nonignorable/reference/quantile_grid.md)
  : Quantile grid for the mean scan

## Distribution function and quantiles

The scan at each value of a grid, pointwise in t, and bounds for
quantiles.

- [`sharp_scan_cdf()`](https://buecker.ms/nonignorable/reference/sharp_scan_cdf.md)
  [`quantile(`*`<sharp_scan_cdf>`*`)`](https://buecker.ms/nonignorable/reference/sharp_scan_cdf.md)
  [`plot(`*`<sharp_scan_cdf>`*`)`](https://buecker.ms/nonignorable/reference/sharp_scan_cdf.md)
  : Debiased scan for the distribution function and its quantiles
- [`autoplot(`*`<sharp_scan_cdf>`*`)`](https://buecker.ms/nonignorable/reference/autoplot.sharp_scan_cdf.md)
  : Distribution-function scan, with ggplot2
- [`crossfit_cdf()`](https://buecker.ms/nonignorable/reference/crossfit_cdf.md)
  : Cross-fitted inputs for the distribution-function scan

## Sampling error

Debiased estimates and confidence sets from cross-fitted inputs.

- [`sharp_scan_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_ci.md)
  : Confidence set for the binary scan, by debiased estimation
- [`sharp_scan_mean_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean_ci.md)
  : Debiased scan for the mean of a non-binary outcome, with a
  confidence set
- [`crossfit_inputs()`](https://buecker.ms/nonignorable/reference/crossfit_inputs.md)
  : Cross-fitted scan inputs for any learner
- [`crossfit_quantiles()`](https://buecker.ms/nonignorable/reference/crossfit_quantiles.md)
  : Cross-fitted quantile grids for the debiased mean scan

## Comparing two groups

- [`sharp_scan_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_contrast.md)
  : Debiased scan for the difference between two portfolios
- [`sharp_scan_mean_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean_contrast.md)
  : Debiased scan for the difference of two portfolio means

## Simulated data

- [`simulate_shadow()`](https://buecker.ms/nonignorable/reference/simulate_shadow.md)
  [`obs_quantiles()`](https://buecker.ms/nonignorable/reference/simulate_shadow.md)
  : Simulate a nonignorably selected outcome with known truth
- [`scan_inputs()`](https://buecker.ms/nonignorable/reference/scan_inputs.md)
  : Extract the two scan inputs from fitted models

## Package

- [`nonignorable`](https://buecker.ms/nonignorable/reference/nonignorable-package.md)
  [`nonignorable-package`](https://buecker.ms/nonignorable/reference/nonignorable-package.md)
  : nonignorable: shadow-variable screens and sharp sensitivity scans
