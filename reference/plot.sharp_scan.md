# Sensitivity curve of a scan

Plots the portfolio rate (or mean, for
[`sharp_scan_mean()`](https://mchlbckr.github.io/nonignorable/reference/sharp_scan_mean.md))
against selection strength: the sharp identified set as a filled band
with solid edges, the outer band, the bounds without the cap \\\pi \le
1\\ (for a rate \\\mathrm{expit}(\ell \pm \lambda)\\), dashed and
lighter, the point that a model on the observed rows asserts (\\\lambda
= 0\\) marked, and, if the scan has a `threshold`, the decision
threshold with a drop line where the sharp set first reaches it. That
crossing is the reason to run the scan: it shows at which strength of
selection the decision flips.

## Usage

``` r
# S3 method for class 'sharp_scan'
plot(
  x,
  scale = c("ratio", "log"),
  colors = scan_colors,
  text_col = "grey20",
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A `"sharp_scan"` or `"sharp_scan_mean"` object.

- scale:

  `"ratio"` (default) puts \\e^\lambda\\, the largest admissible ratio
  of response propensities, on the horizontal axis, starting at 1;
  `"log"` puts \\\lambda\\.

- colors:

  Named colours `sharp` and `naive`. The defaults are validated for
  colour-vision deficiency; on a dark background use
  `c(sharp = "#3987e5", naive = "#d95926")`.

- text_col:

  Colour for all text and reference lines. Text is never drawn in a
  series colour.

- xlab, ylab, main:

  Axis labels and title.

- ...:

  Passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## Details

Where \\\hat r\\ is small, as for rejected applications, the upper edge
of the sharp set coincides with the outer one; see
[`sharp_scan()`](https://mchlbckr.github.io/nonignorable/reference/sharp_scan.md).

The legend goes into the first corner, from top left, that no curve,
band or annotation touches; if none is free, the vertical axis is
extended upwards to make room. Where the cap \\\pi \le 1\\ binds in
every target row from some \\\lambda\\ on, that edge of the set stays
flat from there, and the legend says from where (see `held` in
[`sharp_scan()`](https://mchlbckr.github.io/nonignorable/reference/sharp_scan.md)).

## Examples

``` r
d <- simulate_shadow(5000, tilt = -1, seed = 1)
rej <- d[d$observed == 0, ]
plot(sharp_scan(rej$eta_obs, rej$obs_rate, threshold = 0.35))
```
