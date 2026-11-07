# nonignorable: shadow-variable screens and sharp sensitivity scans

If an outcome is missing not at random, every model fitted on the
observed rows silently asserts *one* value for how observability depends
on the outcome: none. Nothing in the data selects that value. The
package offers two tools for two different questions.

## Details

- [`shadow_screen()`](https://mchlbckr.github.io/nonignorable/reference/shadow_screen.md):
  *is this variable fit to serve as a shadow variable?* Relevance for
  the outcome against violation of the exclusion restriction, on one
  scale. This belongs to a point-identified shadow-variable model.

- [`sharp_scan()`](https://mchlbckr.github.io/nonignorable/reference/sharp_scan.md)
  and
  [`sharp_scan_mean()`](https://mchlbckr.github.io/nonignorable/reference/sharp_scan_mean.md):
  *how far can the answer move?* The sharp identified set under a bound
  on the ratio of response propensities between outcome values, with no
  shadow variable at all; for a binary outcome and for the mean of any
  outcome.

The two rest on different assumptions and are not two halves of one
procedure. The help pages state the formulas, so that every number can
be checked by hand.
[`simulate_shadow()`](https://mchlbckr.github.io/nonignorable/reference/simulate_shadow.md)
provides data with known truth.

## See also

Useful links:

- <https://github.com/mchlbckr/nonignorable>

- Report bugs at <https://github.com/mchlbckr/nonignorable/issues>

## Author

**Maintainer**: Michael B\<c3\>\<bc\>cker
<michael.buecker@fh-muenster.de> \[copyright holder\]

Authors:

- Michael B\<c3\>\<bc\>cker <michael.buecker@fh-muenster.de> \[copyright
  holder\]
