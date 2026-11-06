#' Debiased scan for the distribution function and its quantiles
#'
#' `sharp_scan_cdf()` estimates, at each value `t` of a grid, the sharp
#' identified set for the distribution function \eqn{F(t) = P(Y \le t)} of an
#' outcome that may be missing not at random, under the assumption of
#' [sharp_scan_mean()]: the ratio of the response propensities of any two
#' outcome values is at most \eqn{e^\lambda}. It adds a debiased estimate and
#' a confidence set at each `t`. [quantile()] inverts the estimated bounds
#' into bounds for quantiles.
#'
#' @section The set:
#' Fix \eqn{t} and write \eqn{Z = 1(Y \le t)} and \eqn{\eta_t(x) = P(Y \le t
#' \mid R = 1, X = x)}. The identified set for \eqn{F(t \mid x)} is that of the
#' binary scan for \eqn{Z} (Bücker 2026, Corollary S2):
#' \deqn{[\max\{\mathrm{expit}(\ell_t - \lambda),\, r\eta_t\},\;
#'   \min\{\mathrm{expit}(\ell_t + \lambda),\, r\eta_t + 1 - r\}],
#'   \qquad \ell_t = \mathrm{logit}\,\eta_t,}
#' and the set for \eqn{F(t)} is its mean over rows. Each endpoint is
#' estimated as in [sharp_scan_ci()], whose score for the binary outcome
#' \eqn{Z} it uses; `target`, `segment`, `weights`, `trim` and `se_type` have
#' the same meaning, so `target = "unobserved"` gives the distribution
#' function of the missing outcomes, \eqn{P(Y \le t \mid R = 0)}.
#'
#' **The sets are sharp for each \eqn{t} separately, not as a band.** In
#' general no single law compatible with the data attains the upper bound at
#' every \eqn{t}; the curves are the pointwise envelope. Likewise the
#' confidence sets hold pointwise in \eqn{t}; bands uniform in \eqn{t} are not
#' offered, because the rows where the two pieces of a bound cross move with
#' \eqn{t}. To cover the set at a few values of \eqn{t} at once, pass only
#' those values and `level = 1 - (1 - a) / k` for \eqn{k} values and joint
#' level `a` (Bonferroni); this is valid and conservative.
#'
#' @section Inference and the margin condition:
#' Each endpoint is the mean of a minimum (or maximum) of two pieces that
#' cross at an angle, so, unlike for the mean, the normal approximation needs
#' a margin condition: few rows may lie close to where the pieces cross, for
#' the upper end where \eqn{r(1 + \kappa\eta_t) = 1}, \eqn{\kappa = e^\lambda
#' - 1}. The output reports, at each `t`, the weighted share of rows whose
#' two pieces lie within 0.25 on the logit scale (`near_lower`,
#' `near_upper`); where that share is large the confidence set should be read
#' with caution. Nuisance estimates must be cross-fitted and converge faster
#' than \eqn{n^{-1/4}}.
#'
#' Estimate \eqn{\eta_t(x)} directly, as a classification of \eqn{Z} on
#' \eqn{X} among the observed rows, one fit per `t`; [crossfit_cdf()] does
#' that for any learner. Reading \eqn{\eta_t} off a grid of conditional
#' quantiles is possible but adds the error of the grid.
#'
#' @section Quantiles:
#' For a level \eqn{\alpha}, the bounds \eqn{Q^- = \inf\{t : F^+(t) \ge
#' \alpha\}} and \eqn{Q^+ = \inf\{t : F^-(t) \ge \alpha\}} contain the
#' \eqn{\alpha}-quantile of every compatible law. They are its sharp set when
#' \eqn{Y} given \eqn{R = 1} has a density that is positive on an open
#' interval and zero outside it and \eqn{\eta_t(x)} is continuous in \eqn{t};
#' otherwise they are an outer set. `quantile(x, probs)` computes them from
#' the estimated curves on the grid `t`, after sorting each curve into
#' increasing order (a rearrangement, which leaves a monotone curve unchanged).
#' On the grid they are widened outwards: \eqn{Q^-} is reported as the last
#' grid value at which \eqn{F^+} is still below \eqn{\alpha}, \eqn{Q^+} as
#' the first at which \eqn{F^-} reaches it, so a coarse grid gives a wider
#' interval, never a narrower one. They carry no confidence statement.
#'
#' @param eta_t An \eqn{n \times T} matrix: column \eqn{j} holds the
#'   cross-fitted \eqn{\hat\eta_{t_j}(x_i)} for every row, as from
#'   [crossfit_cdf()]. Values are clipped into `[clip, 1 - clip]`.
#' @param r Cross-fitted \eqn{\hat r(x_i)} for every row, in (0, 1].
#' @param response Response indicator for every row, 0/1 or logical.
#' @param outcome Numeric outcome; ignored (and may be `NA`) where `response`
#'   is 0.
#' @param t Increasing grid of values, one per column of `eta_t`.
#' @param lambda A single non-negative \eqn{\lambda}.
#' @param target,segment,weights,level,trim,se_type As in [sharp_scan_ci()].
#' @param clip \eqn{\hat\eta_t} is kept inside `[clip, 1 - clip]`.
#' @return An object of class `"sharp_scan_cdf"` with `curve`, a data frame
#'   with one row per `t`: `t`, `point` (the estimate at \eqn{\lambda = 0}),
#'   `lower`, `upper`, `se_lower`, `se_upper`, `ci_lower`, `ci_upper`,
#'   `naive_lower`, `naive_upper` (the plug-in outer band, without the cap
#'   \eqn{\pi \le 1}), `near_lower`, `near_upper`; and `lambda`, `level`,
#'   `target`, `n`, `n_target`.
#' @references
#' Bücker, M. (2026). Sharp intervals for a mean under bounded selection
#' ratios when the response rate is known. Manuscript.
#'
#' Chernozhukov, V., Fernández-Val, I. and Galichon, A. (2009). Improving
#' point and interval estimators of monotone functions by rearrangement.
#' *Biometrika*, 96(3), 559--575. \doi{10.1093/biomet/asp030}
#' @seealso [crossfit_cdf()], [sharp_scan_ci()], [sharp_scan_mean_ci()]
#' @export
#' @examples
#' d <- simulate_shadow(2000, tilt = -1, outcome = "continuous", seed = 4)
#' xs <- c("x1", "x2", "z_shadow", "z_policy", "z_score")
#' fit <- function(lhs) function(train, newdata)
#'   predict(glm(reformulate(xs, lhs), binomial, train), newdata,
#'           type = "response")
#' tg <- seq(-2, 2, by = 0.5)
#' cf <- crossfit_cdf(d, "observed", "y", tg, fit_eta = fit("z"),
#'                    fit_r = fit("observed"), folds = 3, seed = 1)
#' k <- sharp_scan_cdf(cf$eta_t, cf$r, d$observed, d$y, t = tg, lambda = 0.5)
#' k
#' quantile(k, c(0.25, 0.5, 0.75))
#' plot(k)
sharp_scan_cdf <- function(eta_t, r, response, outcome, t, lambda = 1,
                           target = "all", weights = NULL, level = 0.9,
                           trim = 0, se_type = c("max", "conditional", "empirical"),
                           segment = NULL, clip = 1e-6) {
  se_type <- match.arg(se_type)
  if (is.null(dim(eta_t))) eta_t <- matrix(eta_t, ncol = 1L)
  eta_t <- as.matrix(eta_t)
  n <- nrow(eta_t)
  if (!is.numeric(eta_t) || anyNA(eta_t) || any(eta_t < 0 | eta_t > 1))
    stop("`eta_t` must be a numeric matrix with values in [0, 1], without NA.",
         call. = FALSE)
  if (!is.numeric(t) || length(t) != ncol(eta_t) || anyNA(t) || any(!is.finite(t)))
    stop("`t` must be finite numbers, one per column of `eta_t`.", call. = FALSE)
  if (is.unsorted(t, strictly = TRUE))
    stop("`t` must be strictly increasing.", call. = FALSE)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
    stop("`lambda` must be a single non-negative number.", call. = FALSE)
  if (!is.numeric(clip) || length(clip) != 1L || !(clip > 0 && clip < 0.5))
    stop("`clip` must be a single number in (0, 0.5).", call. = FALSE)
  if (is.logical(response)) response <- as.integer(response)
  if (!is.numeric(response) || length(response) != n || anyNA(response) ||
      !all(response %in% c(0, 1)))
    stop("`response` must be 0/1 or logical, one per row of `eta_t`, without NA.",
         call. = FALSE)
  if (length(outcome) != n || !is.numeric(outcome))
    stop("`outcome` must be numeric, one value per row.", call. = FALSE)
  R <- response
  if (anyNA(outcome[R == 1]))
    stop("`outcome` must not be NA where `response` is 1.", call. = FALSE)
  if (!is.numeric(r) || length(r) != n)
    stop("`r` must have one value per row of `eta_t`.", call. = FALSE)

  # rows whose eta_t(x) enters the set: all rows of the target or segment
  # (also the observed ones for "unobserved"), less the trimmed ones
  entering <- (if (is.logical(target)) target
               else if (!is.null(segment)) segment
               else rep(TRUE, n)) & r >= trim
  lams <- unique(c(0, lambda))
  i_lam <- length(lams)
  warns <- character()
  rows <- vector("list", length(t))
  first <- NULL
  for (j in seq_along(t)) {
    z <- ifelse(R == 1, as.integer(outcome <= t[j]), NA_integer_)
    e <- pmin(pmax(eta_t[, j], clip), 1 - clip)
    k <- withCallingHandlers(
      sharp_scan_ci(e, r, R, z, target = target, lambda = lams, weights = weights,
                    level = level, trim = trim, se_type = se_type, segment = segment),
      warning = function(w) {
        warns <<- c(warns, sprintf("at t = %s: %s", format(t[j], digits = 4),
                                   conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    if (is.null(first)) first <- k
    cv <- k$curve[i_lam, ]
    rows[[j]] <- data.frame(
      t = t[j], point = k$point, lower = cv$lower, upper = cv$upper,
      se_lower = cv$se_lower, se_upper = cv$se_upper,
      ci_lower = cv$ci_lower, ci_upper = cv$ci_upper,
      naive_lower = cv$naive_lower, naive_upper = cv$naive_upper,
      near_lower = NA_real_, near_upper = NA_real_)
    near <- near_kink(e, r, lambda, entering, weights)
    rows[[j]]$near_lower <- near[["lower"]]
    rows[[j]]$near_upper <- near[["upper"]]
  }
  if (length(warns))
    for (w in utils::head(warns, 3L)) warning(w, call. = FALSE)
  structure(list(
    curve = do.call(rbind, rows), lambda = lambda, level = level,
    target = first$target, se_type = se_type, n = n, n_target = first$n_target
  ), class = "sharp_scan_cdf")
}

# Weighted share of entering rows whose two pieces of each bound lie within
# `width` on the logit scale: the mass near the kink.
near_kink <- function(e, r, lambda, entering, weights, width = 0.25) {
  w <- if (is.null(weights)) rep(1, length(e)) else weights
  w <- w * entering
  ell <- stats::qlogis(e)
  lo_piece <- stats::qlogis(pmin(pmax(r * e, 1e-12), 1 - 1e-12))
  up_piece <- stats::qlogis(pmin(pmax(r * e + 1 - r, 1e-12), 1 - 1e-12))
  c(lower = sum(w * (abs(ell - lambda - lo_piece) <= width)) / sum(w),
    upper = sum(w * (abs(ell + lambda - up_piece) <= width)) / sum(w))
}

#' @rdname sharp_scan_cdf
#' @param x A `"sharp_scan_cdf"` object.
#' @param probs Levels \eqn{\alpha} in (0, 1).
#' @param ... Ignored.
#' @return `quantile()`: a data frame with `prob`, `point` (the quantile of
#'   the estimate at \eqn{\lambda = 0}, interpolated linearly between grid
#'   values), `lower` (\eqn{Q^-}) and `upper` (\eqn{Q^+}). `NA` means the
#'   quantile is not bracketed by the grid: at or below its first value, or
#'   above its last.
#' @importFrom stats quantile
#' @export
quantile.sharp_scan_cdf <- function(x, probs = c(0.25, 0.5, 0.75), ...) {
  if (!is.numeric(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1))
    stop("`probs` must lie in (0, 1).", call. = FALSE)
  cv <- x$curve
  tg <- cv$t
  # first grid index where the rearranged curve reaches a, NA if none
  hit <- function(f, a) { j <- which(sort(f) >= a); if (length(j)) j[1] else NA_integer_ }
  lo <- function(a) {               # Q- lies in (t[j-1], t[j]]: report t[j-1]
    j <- hit(cv$upper, a)
    if (is.na(j) || j == 1L) NA_real_ else tg[j - 1L]
  }
  up <- function(a) {               # Q+ lies in (t[j-1], t[j]]: report t[j]
    j <- hit(cv$lower, a)
    if (is.na(j) || j == 1L) NA_real_ else tg[j]
  }
  pt <- function(a) {               # linear interpolation of the MAR curve
    f <- sort(cv$point); j <- hit(f, a)
    if (is.na(j) || j == 1L) return(NA_real_)
    tg[j - 1L] + (a - f[j - 1L]) / (f[j] - f[j - 1L]) * (tg[j] - tg[j - 1L])
  }
  data.frame(prob = probs, point = vapply(probs, pt, 0),
             lower = vapply(probs, lo, 0), upper = vapply(probs, up, 0))
}

#' @export
print.sharp_scan_cdf <- function(x, digits = 3, ...) {
  cv <- x$curve
  cat(sprintf(paste0("<sharp_scan_cdf> sharp set for the distribution function, ",
                     "debiased, with %g%% confidence set\n"), 100 * x$level))
  cat(sprintf("  rows: %d   target: %s (%d rows)\n", x$n, x$target, x$n_target))
  cat(sprintf("  lambda = %s (propensity ratio %s); %d values of t, sets pointwise in t\n",
              format(x$lambda, digits = digits), format(exp(x$lambda), digits = digits),
              nrow(cv)))
  show <- unique(round(seq(1, nrow(cv), length.out = min(nrow(cv), 7))))
  tab <- data.frame(t = cv$t, `MAR` = cv$point, lower = cv$lower, upper = cv$upper,
                    ci_lower = cv$ci_lower, ci_upper = cv$ci_upper,
                    check.names = FALSE)[show, ]
  print(format(tab, digits = digits), row.names = FALSE)
  nk <- pmax(cv$near_lower, cv$near_upper)
  if (max(nk) > 0.1) {
    i <- which.max(nk)
    cat(strwrap(sprintf(paste0(
      "near the kink: at t = %s, %.0f%% of the rows have both pieces of a bound ",
      "within 0.25 on the logit scale; the normal approximation needs few such ",
      "rows (margin condition), so read the confidence set there with caution"), format(cv$t[i], digits = digits), 100 * nk[i]),
      indent = 2, exdent = 4), sep = "\n")
  }
  invisible(x)
}

#' @rdname sharp_scan_cdf
#' @param colors,text_col Colours, as for [plot.sharp_scan()].
#' @param xlab,ylab,main Axis labels and title.
#' @export
plot.sharp_scan_cdf <- function(x, colors = scan_colors, text_col = "grey20",
                                xlab = "t", ylab = NULL, main = NULL, ...) {
  cv <- x$curve
  if (is.null(ylab))
    ylab <- if (startsWith(x$target, "unobserved")) "P(Y <= t | R = 0)" else "P(Y <= t)"
  ylim <- range(0, 1, cv$ci_lower, cv$ci_upper, cv$naive_lower, cv$naive_upper)
  graphics::plot(cv$t, cv$upper, type = "n", ylim = ylim, xlab = "", ylab = "",
                 axes = FALSE, ...)
  graphics::axis(1, col = text_col, col.axis = text_col)
  graphics::axis(2, las = 1, col = text_col, col.axis = text_col)
  graphics::box(bty = "l", col = text_col)
  graphics::title(xlab = xlab, ylab = ylab, main = main, col.lab = text_col,
                  col.main = text_col)
  naive_col <- grDevices::adjustcolor(colors[["naive"]], 0.75)
  graphics::polygon(c(cv$t, rev(cv$t)), c(cv$lower, rev(cv$upper)),
                    col = grDevices::adjustcolor(colors[["sharp"]], 0.18), border = NA)
  graphics::lines(cv$t, cv$naive_lower, lty = 2, lwd = 1.5, col = naive_col)
  graphics::lines(cv$t, cv$naive_upper, lty = 2, lwd = 1.5, col = naive_col)
  graphics::lines(cv$t, cv$lower, lwd = 2, col = colors[["sharp"]])
  graphics::lines(cv$t, cv$upper, lwd = 2, col = colors[["sharp"]])
  graphics::lines(cv$t, cv$ci_lower, lwd = 1.3, lty = 3, col = colors[["sharp"]])
  graphics::lines(cv$t, cv$ci_upper, lwd = 1.3, lty = 3, col = colors[["sharp"]])
  graphics::lines(cv$t, cv$point, lwd = 1.5, col = text_col)
  graphics::legend("topleft", bty = "n", text.col = text_col, cex = 0.85,
                   legend = c(sprintf("Sharp set at e^lambda = %s, pointwise in t",
                                      format(exp(x$lambda), digits = 3)),
                              "Outer band (without pi <= 1)",
                              sprintf("%g%% confidence set, pointwise", 100 * x$level),
                              "Model on observed rows (lambda = 0)"),
                   fill = c(grDevices::adjustcolor(colors[["sharp"]], 0.18), NA, NA, NA),
                   border = c(colors[["sharp"]], NA, NA, NA),
                   lty = c(NA, 2, 3, 1), lwd = c(NA, 1.5, 1.3, 1.5),
                   col = c(NA, naive_col, colors[["sharp"]], text_col), merge = TRUE)
  invisible(x)
}

#' Cross-fitted inputs for the distribution-function scan
#'
#' Splits the rows into folds and predicts, for each fold from models fitted
#' on the other folds, \eqn{\hat\eta_t(x) = P(Y \le t \mid R = 1, X = x)} at
#' every value of `t`, by a classification of \eqn{Z = 1(Y \le t)} among the
#' observed rows, and \eqn{\hat r}, as [sharp_scan_cdf()] requires. The
#' learner is yours: `fit_eta(train, newdata)` receives the observed rows of
#' the training folds with a 0/1 column `z` added and returns probabilities
#' for `newdata`; it is called once per fold and value of `t`. Where all
#' training rows share one value of `z` (a `t` below or above the observed
#' outcomes), the prediction is that value and `fit_eta` is not called.
#'
#' @param data A data frame with all rows, observed and unobserved.
#' @param response Name of the 0/1 response column.
#' @param outcome Name of the outcome column.
#' @param t Increasing grid of values.
#' @param fit_eta,fit_r Functions `(train, newdata)` returning probabilities.
#' @param folds Number of folds. Default 5.
#' @param seed Optional seed for the fold assignment.
#' @param clip Bound for \eqn{\hat r} and \eqn{\hat\eta_t}, as in
#'   [crossfit_inputs()].
#' @return A list with `eta_t` (an \eqn{n \times T} matrix), `r`, `fold` and
#'   `t`.
#' @seealso [sharp_scan_cdf()], [crossfit_inputs()]
#' @export
crossfit_cdf <- function(data, response, outcome, t, fit_eta, fit_r, folds = 5L,
                         seed = NULL, clip = 1e-6) {
  if (!is.function(fit_eta) || !is.function(fit_r))
    stop("`fit_eta` and `fit_r` must be functions(train, newdata).", call. = FALSE)
  if (!outcome %in% names(data))
    stop("Column `", outcome, "` not found.", call. = FALSE)
  if (!is.numeric(t) || !length(t) || anyNA(t) || is.unsorted(t, strictly = TRUE))
    stop("`t` must be a strictly increasing numeric vector.", call. = FALSE)
  if ("z" %in% names(data))
    stop("`data` must not have a column `z`; it is used for 1(outcome <= t).",
         call. = FALSE)
  cs <- crossfit_setup(data, response, folds, seed)
  n <- cs$n; R <- cs$R; fold <- cs$fold
  y <- data[[outcome]]
  if (!is.numeric(y) || anyNA(y[R == 1]))
    stop("Column `", outcome, "` must be numeric and not NA where observed.",
         call. = FALSE)
  eta <- matrix(NA_real_, n, length(t))
  r <- rep(NA_real_, n)
  for (k in seq_len(cs$folds)) {
    tr <- data[fold != k, , drop = FALSE]
    te <- data[fold == k, , drop = FALSE]
    pr <- as.numeric(fit_r(tr, te))
    if (length(pr) != nrow(te))
      stop("`fit_r` must return one prediction per row of `newdata`.", call. = FALSE)
    r[fold == k] <- pr
    obs <- tr[R[fold != k] == 1, , drop = FALSE]
    for (j in seq_along(t)) {
      obs$z <- as.integer(obs[[outcome]] <= t[j])
      pe <- if (length(unique(obs$z)) == 1L) rep(obs$z[1], nrow(te))
            else as.numeric(fit_eta(obs, te))
      if (length(pe) != nrow(te))
        stop("`fit_eta` must return one prediction per row of `newdata`.", call. = FALSE)
      eta[fold == k, j] <- pe
    }
  }
  if (anyNA(eta) || anyNA(r))
    stop("A fit function returned missing predictions.", call. = FALSE)
  warn_separation(R, r, clip)
  list(eta_t = pmin(pmax(eta, clip), 1 - clip), r = pmin(pmax(r, clip), 1),
       fold = fold, t = t)
}
