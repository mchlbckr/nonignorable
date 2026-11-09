# Sensitivity curve of a scan, with ggplot2

The ggplot2 version of
[`plot.sharp_scan()`](https://buecker.ms/nonignorable/reference/plot.sharp_scan.md),
with the same information: the sharp identified set as a filled band
with solid edges, the outer band (without the cap \\\pi \le 1\\) dashed,
the confidence set dotted where the scan has one, the point that a model
on the observed rows asserts (\\\lambda = 0\\), and, if the scan has a
`threshold`, the decision threshold with a drop line and a label where
the set first reaches it. Where the cap \\\pi \le 1\\ holds an edge of
the set flat from some \\\lambda\\ on (see `held` in
[`sharp_scan()`](https://buecker.ms/nonignorable/reference/sharp_scan.md)),
a diamond on that edge marks where, and the legend gives the value.

## Usage

``` r
# S3 method for class 'sharp_scan'
autoplot(
  object,
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

- object:

  A `"sharp_scan"` object: from
  [`sharp_scan()`](https://buecker.ms/nonignorable/reference/sharp_scan.md),
  [`sharp_scan_mean()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean.md),
  [`sharp_scan_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_ci.md),
  [`sharp_scan_mean_ci()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean_ci.md),
  [`sharp_scan_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_contrast.md)
  or
  [`sharp_scan_mean_contrast()`](https://buecker.ms/nonignorable/reference/sharp_scan_mean_contrast.md).

- scale:

  `"ratio"` (default) puts \\e^\lambda\\ on the horizontal axis, `"log"`
  puts \\\lambda\\, as in
  [`plot.sharp_scan()`](https://buecker.ms/nonignorable/reference/plot.sharp_scan.md).

- colors:

  Named colours `sharp` and `naive`, as in
  [`plot.sharp_scan()`](https://buecker.ms/nonignorable/reference/plot.sharp_scan.md).

- text_col:

  Colour for text and reference lines.

- xlab, ylab, main:

  Axis labels and title; `NULL` gives the same defaults as
  [`plot.sharp_scan()`](https://buecker.ms/nonignorable/reference/plot.sharp_scan.md).
  `+ ggplot2::labs()` works too.

- ...:

  Not used.

## Value

A ggplot object.

## Details

The legend sits below the panel. The result is an ordinary ggplot
object: add layers, scales,
[`labs()`](https://ggplot2.tidyverse.org/reference/labs.html) or a
[`theme()`](https://ggplot2.tidyverse.org/reference/theme.html) to it as
usual. The methods are registered for
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
when ggplot2 is loaded; ggplot2 is not required to install or use the
rest of the package.

The curve layers carry their data in long form:
[`layer_data()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)
or `p$layers[[i]]$data` gives `lambda`, the horizontal position `x`, the
value `y`, the `series` (`"sharp"`, `"naive"`, `"ci"`) and the `edge`.

## See also

[`plot.sharp_scan()`](https://buecker.ms/nonignorable/reference/plot.sharp_scan.md)
for base graphics.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  d <- simulate_shadow(5000, tilt = -1, seed = 1)
  rej <- d[d$observed == 0, ]
  s <- sharp_scan(rej$eta_obs, rej$obs_rate, threshold = 0.35)
  autoplot(s)
  autoplot(s, scale = "log") + theme_minimal()
}
```
