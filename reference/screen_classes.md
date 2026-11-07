# Classes for a numeric characteristic, cut only between distinct values

`screen_classes()` turns a numeric characteristic into a factor of about
`k` classes of roughly equal mass for
[`shadow_screen()`](https://mchlbckr.github.io/nonignorable/reference/shadow_screen.md).
Unlike quantile classes from [`cut()`](https://rdrr.io/r/base/cut.html)
and [`quantile()`](https://rdrr.io/r/stats/quantile.html), a class
boundary never falls on a value that occurs in the data, so no
observation depends on whether intervals are closed on the left or the
right. That matters for characteristics with heaped values: loan amounts
cluster on round numbers, and quantiles then land exactly on those
heaps.

## Usage

``` r
screen_classes(x, k = 10, support = range(x, na.rm = TRUE))
```

## Arguments

- x:

  A numeric vector. `NA` stays `NA`.

- k:

  Target number of classes inside the support. Default 10.

- support:

  Length-two range; values outside it get the classes `"below"` and
  `"above"`. Default: the range of `x`, i.e. none.

## Value

A factor, levels ordered from below to above, with attribute
`"heaping"`: the share of rows inside the support that lie on a value at
which [`quantile()`](https://rdrr.io/r/stats/quantile.html)-based
`k`-classes would place a boundary. Above a few percent, quantile
classes would assign a material part of the sample by convention.

## Details

Values outside `support` get classes of their own, `"below"` and
`"above"`. Pass as `support` the range of the characteristic among the
rows whose outcome is observed: a lender that never accepts below a
score cut-off then shows up as one class with no accepted row, which is
the step that a linear term in the response model cannot represent.
[`shadow_screen()`](https://mchlbckr.github.io/nonignorable/reference/shadow_screen.md)
does this itself with `numeric = "classes"`.

Classes are built greedily from the smallest value up: a class is closed
once it holds at least its share \\1/k\\ of the rows inside the support.
A single value holding more than that share forms a class of its own, so
class masses vary; a last class smaller than half its share is merged
into the one before. With at most `k` distinct values every value is its
own class.

## See also

[`shadow_screen()`](https://mchlbckr.github.io/nonignorable/reference/shadow_screen.md)

## Examples

``` r
amount <- c(rep(10000, 30), rep(5000, 20), round(runif(50, 1000, 35000), -2))
f <- screen_classes(amount, k = 5)
table(f)
#> f
#>   [1,300, 5,000]  [6,300, 10,000] [10,900, 24,100] [24,400, 34,700] 
#>               27               37               20               16 
attr(f, "heaping")
#> [1] 0.5
# a cut-off: nothing below 662 was ever accepted
score <- round(runif(200, 600, 800))
table(screen_classes(score, k = 4, support = c(662, 800)))
#> 
#>      below [663, 702] [704, 735] [736, 771] [772, 800] 
#>         56         37         36         37         34 
```
