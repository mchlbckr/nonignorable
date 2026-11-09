# Distribution-function scan, with ggplot2

The ggplot2 version of the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method of
[`sharp_scan_cdf()`](https://buecker.ms/nonignorable/reference/sharp_scan_cdf.md):
the sharp set for \\P(Y \le t)\\ at the scan's \\\lambda\\ as a filled
band with solid edges, the outer band dashed, the confidence set dotted
and the estimate under MAR (\\\lambda = 0\\) as a line, all pointwise in
\\t\\. The curve layers carry `t`, `y`, `series` and `edge` in long
form.

## Usage

``` r
# S3 method for class 'sharp_scan_cdf'
autoplot(
  object,
  colors = scan_colors,
  text_col = "grey20",
  xlab = "t",
  ylab = NULL,
  main = NULL,
  ...
)
```

## Arguments

- object:

  A `"sharp_scan_cdf"` object.

- colors:

  Named colours `sharp` and `naive`, as in
  [`plot.sharp_scan()`](https://buecker.ms/nonignorable/reference/plot.sharp_scan.md).

- text_col:

  Colour for text and reference lines.

- xlab, ylab, main:

  Axis labels and title; `ylab = NULL` gives the same default as
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- ...:

  Not used.

## Value

A ggplot object.

## See also

[`sharp_scan_cdf()`](https://buecker.ms/nonignorable/reference/sharp_scan_cdf.md)
for [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
