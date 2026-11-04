#' Classes for a numeric characteristic, cut only between distinct values
#'
#' `screen_classes()` turns a numeric characteristic into a factor of about
#' `k` classes of roughly equal mass for [shadow_screen()]. Unlike quantile
#' classes from [cut()] and [quantile()], a class boundary never falls on a
#' value that occurs in the data, so no observation depends on whether
#' intervals are closed on the left or the right. That matters for
#' characteristics with heaped values: loan amounts cluster on round numbers,
#' and quantiles then land exactly on those heaps.
#'
#' Values outside `support` get classes of their own, `"below"` and
#' `"above"`. Pass as `support` the range of the characteristic among the rows
#' whose outcome is observed: a lender that never accepts below a score
#' cut-off then shows up as one class with no accepted row, which is the step
#' that a linear term in the response model cannot represent.
#' [shadow_screen()] does this itself with `numeric = "classes"`.
#'
#' Classes are built greedily from the smallest value up: a class is closed
#' once it holds at least its share \eqn{1/k} of the rows inside the support.
#' A single value holding more than that share forms a class of its own, so
#' class masses vary; a last class smaller than half its share is merged into
#' the one before. With at most `k` distinct values every value is its own
#' class.
#'
#' @param x A numeric vector. `NA` stays `NA`.
#' @param k Target number of classes inside the support. Default 10.
#' @param support Length-two range; values outside it get the classes
#'   `"below"` and `"above"`. Default: the range of `x`, i.e. none.
#' @return A factor, levels ordered from below to above, with attribute
#'   `"heaping"`: the share of rows inside the support that lie on a value at
#'   which [quantile()]-based `k`-classes would place a boundary. Above a few
#'   percent, quantile classes would assign a material part of the sample by
#'   convention.
#' @seealso [shadow_screen()]
#' @export
#' @examples
#' amount <- c(rep(10000, 30), rep(5000, 20), round(runif(50, 1000, 35000), -2))
#' f <- screen_classes(amount, k = 5)
#' table(f)
#' attr(f, "heaping")
#' # a cut-off: nothing below 662 was ever accepted
#' score <- round(runif(200, 600, 800))
#' table(screen_classes(score, k = 4, support = c(662, 800)))
screen_classes <- function(x, k = 10, support = range(x, na.rm = TRUE)) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (!is.numeric(k) || length(k) != 1L || k < 1)
    stop("`k` must be a single number >= 1.", call. = FALSE)
  if (!is.numeric(support) || length(support) != 2L || anyNA(support) ||
      support[1] > support[2])
    stop("`support` must be an increasing range of two numbers.", call. = FALSE)
  k <- as.integer(k)
  ok <- !is.na(x)
  inside <- ok & x >= support[1] & x <= support[2]
  xi <- x[inside]
  out <- rep(NA_character_, length(x))
  heaping <- 0
  levs <- character()
  if (length(xi)) {
    u <- sort(unique(xi))                      # exact values, not table() names
    m <- tabulate(match(xi, u), nbins = length(u))
    if (length(u) <= k) {
      cls <- seq_along(u)
    } else {
      target <- length(xi) / k
      cls <- integer(length(u)); j <- 1L; acc <- 0
      for (i in seq_along(u)) {
        cls[i] <- j; acc <- acc + m[i]
        if (acc >= target && i < length(u)) { j <- j + 1L; acc <- 0 }
      }
      last <- cls == max(cls)
      if (max(cls) > 1L && sum(m[last]) < target / 2) cls[last] <- max(cls) - 1L
    }
    lo <- tapply(u, cls, min); hi <- tapply(u, cls, max)
    fmt <- function(v) vapply(v, function(a) format(signif(a, 6), trim = TRUE,
                              big.mark = ",", scientific = FALSE,
                              drop0trailing = TRUE), "")
    lab <- ifelse(lo == hi, fmt(lo), paste0("[", fmt(lo), ", ", fmt(hi), "]"))
    lab <- make.unique(unname(lab))
    out[inside] <- lab[cls[match(xi, u)]]
    levs <- lab
    qb <- stats::quantile(xi, seq_len(k - 1L) / k, names = FALSE)
    heaping <- mean(xi %in% qb)
  }
  below <- ok & x < support[1]; above <- ok & x > support[2]
  out[below] <- "below"; out[above] <- "above"
  levs <- c(if (any(below)) "below", levs, if (any(above)) "above")
  f <- factor(out, levels = levs)
  attr(f, "heaping") <- heaping
  f
}
