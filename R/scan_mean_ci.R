# Orthogonal score for the upper bound of the sharp set of a mean, per row
# (mnar-sharp notes/09, step 5). `side` is a sorted grid from mean_prep():
# `ys` in decreasing order and `S`, the tail sums U(j/K). With
# s° = min(s*, sbar), q = Q(1 - s°), theta = g(s°), D = 1 + kappa s°,
#   phi = theta + R {Y + kappa (Y - q)_+ + kappa q s° - theta D} / (r D)
#         - 1{sbar < s*} (q - theta)(R - r) / (r^2 D).
# The lower bound is the same for -Y, negated. Y is only used where R = 1.
mean_score <- function(side, m, r, Y, R, lam) {
  kap <- expm1(lam)
  n <- length(m)
  if (kap == 0)
    return(list(phi = m + R * (Y - m) / r, theta = m, cap = rep(FALSE, n)))
  K <- ncol(side$ys)
  P <- (0:K) / K
  g <- (m + kap * side$S) / (1 + kap * rep(P, each = n))
  # s* = jstar / K, exactly: the first kink within rounding of the maximum, so
  # that a flat g (a point mass, where s* is not defined) gives s* = 0 instead
  # of whichever kink floating-point noise favours
  gmax <- g[cbind(seq_len(n), max.col(g, ties.method = "first"))]
  jstar <- max.col(g >= gmax - 1e-10 * pmax(1, abs(gmax)), ties.method = "first") - 1L
  sbar <- pmin(1, (1 - r) / (r * kap))
  cap <- sbar < jstar / K
  so <- ifelse(cap, sbar, jstar / K)
  j <- ifelse(cap, pmin(floor(sbar * K), K - 1L), pmin(jstar, K - 1L))
  idx <- cbind(seq_len(n), j + 1L)
  q <- side$ys[idx]
  U <- side$S[idx] + (so - j / K) * q
  D <- 1 + kap * so
  theta <- (m + kap * U) / D
  phi <- theta + R * (Y + kap * pmax(Y - q, 0) + kap * q * so - theta * D) / (r * D) -
    cap * (q - theta) * (R - r) / (r^2 * D)
  list(phi = phi, theta = theta, cap = cap)
}

#' Debiased scan for the mean of a non-binary outcome, with a confidence set
#'
#' `sharp_scan_mean_ci()` estimates the sharp identified set of
#' [sharp_scan_mean()] -- the mean of an outcome that may be missing not at
#' random, when the ratio of the response propensities of any two outcome
#' values is at most \eqn{e^\lambda} -- with a debiased estimator and a
#' confidence set. It is the counterpart of [sharp_scan_ci()] for an outcome
#' that is not binary.
#'
#' @section The estimator:
#' Each endpoint of the set is a mean over rows of \eqn{\mu^\pm(x)}, which
#' depends on \eqn{r(x)} and on the quantile function of \eqn{Y} given
#' \eqn{R = 1, X = x}. The plug-in, [sharp_scan_mean()], inherits the
#' first-order error of both estimates; in simulations it was typically biased
#' towards the inside of the set. The estimator here adds to each row a correction term that
#' removes that first-order error: an inverse-probability weighted residual
#' of the observed outcome, and, where the constraint \eqn{\pi \le 1} binds,
#' a term in \eqn{R - \hat r}. The resulting score is Neyman-orthogonal in
#' both regimes, so the estimate is insensitive to small errors in the
#' nuisances, and the correction in \eqn{r} vanishes continuously where the
#' regime changes (on the grid of quantile levels below: up to an error of
#' order \eqn{1/K}), so no margin condition is needed. The confidence set is
#' the Wald interval from the empirical standard deviation of the score.
#'
#' The quantiles enter as a grid of \eqn{K} levels \eqn{(k - 1/2)/K}, as in
#' [sharp_scan_mean()]. An estimator that uses the full empirical
#' distribution of residuals instead differs from this one by
#' \eqn{O(1/K)}.
#'
#' @section What is covered:
#' At each \eqn{\lambda} separately, the interval
#' \eqn{[\hat\psi^- - z\,\widehat{se}^-,\ \hat\psi^+ + z\,\widehat{se}^+]}
#' covers the identified set \eqn{[\psi^-, \psi^+]} with probability at
#' least `level` asymptotically. It is a confidence set for the set, not for
#' the mean inside it, and it holds pointwise in \eqn{\lambda}. The theory
#' requires overlap, \eqn{r(x)} bounded away from 0 and 1, which is why
#' \eqn{\hat r} is kept inside `r_bounds`; a density of \eqn{Y} given
#' \eqn{R = 1, X = x} that is bounded below on an inner range of quantiles;
#' no point mass at the optimum, so that \eqn{s^*} lies strictly inside
#' \eqn{(0, 1)} (for a flat grid the code takes \eqn{s^* = 0}); a bounded
#' conditional second moment of \eqn{Y}; a non-degenerate score variance;
#' and nuisance estimates whose errors
#' vanish faster than \eqn{n^{-1/4}}, cross-fitted as by
#' [crossfit_quantiles()].
#'
#' In the package's own simulation (`simulate_shadow(outcome = "continuous",
#' tilt = -1)`, learners as in the example below with smooth scale, `K = 200`,
#' five folds, `level = 0.9`), the interval covered the whole identified set
#' in 88.7%, 91.1% and 92.7% of 1000 samples of 5000 rows at
#' \eqn{\lambda} = 0.5, 1 and 2, and in 87.3%, 89.0% and 91.0% of 300 samples
#' of 2000 rows.
#'
#' The learner matters for width, not for validity alone. A smooth location
#' and scale model for \eqn{Y} given \eqn{X}, such as `mgcv::gaulss()`, with
#' the quantiles of the standardised residuals, is a good default. In
#' simulations, a linear model for the location only, which cannot follow
#' heteroskedasticity, made the set markedly too wide, and quantile forests
#' gave valid but conservative sets.
#'
#' The target is the mean over all rows (`target = "all"`), over a segment
#' chosen by \eqn{X} (a logical vector), or the mean of the missing outcomes,
#' \eqn{E(Y \mid R = 0)} (`target = "unobserved"`, within a `segment` if
#' given). For the last, \eqn{E(Y \mid R = 0) = \{E(Y) - E(RY)\}/\mathrm{pr}(R
#' = 0)}, so the sharp set is \eqn{\theta^\pm = \{\psi^\pm -
#' E(RY)\}/\mathrm{pr}(R = 0)} with \eqn{\psi^\pm} the endpoints for all rows,
#' estimated with sample means in place of the expectations; the influence
#' function is \eqn{\{\varphi^\pm - RY - \theta^\pm(1 - R)\}/\mathrm{pr}(R =
#' 0)}, and the standard error the sample standard deviation of its
#' estimate, under the same conditions as above (Bücker 2026, corollary on
#' the mean of the nonrespondents). It is not the mean of \eqn{E(Y \mid X =
#' x)} over the unobserved rows. In the package's simulation (as above, five
#' folds, 200 samples each of 2,000 and 5,000 rows, \eqn{\lambda} from 0 to 2),
#' the interval covered the set for \eqn{E(Y \mid R = 0)} in 89.5% to 92.5%
#' of samples, and the standard error matched the spread of the estimates
#' to within 6%.
#'
#' @param q An \eqn{n \times K} matrix of conditional quantiles of \eqn{Y}
#'   given \eqn{R = 1, X = x_i} for **every** row, fitted out of fold, at the
#'   levels \eqn{(k - 1/2)/K}.
#' @param r Cross-fitted \eqn{\hat r(x_i)} for every row, in (0, 1].
#' @param response Response indicator for every row, 0/1 or logical.
#' @param outcome Outcome; ignored (and may be `NA`) where `response` is 0.
#' @param target `"all"` (default), `"unobserved"` (the mean of the missing
#'   outcomes, \eqn{E(Y \mid R = 0)}), or a logical vector selecting the rows
#'   of the portfolio by \eqn{X}.
#' @param segment Optional logical vector chosen by \eqn{X}, one per row, that
#'   restricts a target of `"unobserved"` (or `"all"`) to one segment, as in
#'   [sharp_scan_ci()].
#' @param lambda Grid of \eqn{\lambda} values.
#' @param weights Optional non-negative row weights, a function of \eqn{X}.
#' @param threshold Optional value of the mean at which the decision changes.
#' @param level Confidence level. Default 0.9.
#' @param trim Rows with \eqn{\hat r <} `trim` leave the target. Default 0.
#' @param r_bounds \eqn{\hat r} is clipped into this range before use.
#'   Default `c(0.01, 0.99)`. The number of clipped target rows is reported,
#'   below and above apart: clipping at the top only stabilises the cap
#'   \eqn{\pi \le 1} for rows that are almost always observed, clipping at the
#'   bottom limits inverse weights.
#' @return An object of class
#'   `c("sharp_scan_mean_ci", "sharp_scan_ci", "sharp_scan_mean",
#'   "sharp_scan")` with `curve` (columns `lambda`, `lower`, `upper`,
#'   `se_lower`, `se_upper`, `ci_lower`, `ci_upper`, `naive_lower`,
#'   `naive_upper`; the naive columns are the plug-in outer band), `point`
#'   (the estimate at \eqn{\lambda = 0}), `flip` (`sharp`, `ci`, `naive`),
#'   `trimmed`, `clipped` (rows entering the set, that is all rows of the
#'   segment for `"unobserved"`), `capped` (the weighted share of the rows
#'   entering the set in which the cap \eqn{\pi \le 1} binds at the largest \eqn{\lambda}, per
#'   bound), `tail`, `n` and `n_target`. Estimates and confidence set are kept
#'   monotone in \eqn{\lambda}, as in [sharp_scan_ci()]. `as.data.frame(x,
#'   lambda = )` gives the plug-in bounds per row, with \eqn{\hat r} clipped
#'   into `r_bounds`; they carry no correction and no confidence statement.
#' @references
#' Bücker, M. (2026). Sharp intervals for a mean under bounded selection
#' ratios when the response rate is known. Manuscript.
#' @seealso [sharp_scan_mean()], [crossfit_quantiles()], [sharp_scan_ci()]
#' @export
#' @examples
#' d <- simulate_shadow(1000, tilt = -1, outcome = "continuous", seed = 4)
#' xs <- c("x1", "x2", "z_shadow", "z_policy", "z_score")
#' fit_r <- function(train, newdata)
#'   predict(glm(reformulate(xs, "observed"), binomial, train), newdata,
#'           type = "response")
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   # location smooth in every covariate, scale linear; quantiles from the
#'   # standardised residuals
#'   sm <- paste0("s(", xs, ")")
#'   fit_q <- function(train, newdata, tau) {
#'     m <- mgcv::gam(list(reformulate(sm, "y"), reformulate(xs)),
#'                    family = mgcv::gaulss(), data = train)
#'     p_tr <- predict(m, train, type = "response")
#'     p_new <- predict(m, newdata, type = "response")
#'     e <- (train$y - p_tr[, 1]) * p_tr[, 2]
#'     outer(p_new[, 1], rep(1, length(tau))) +
#'       outer(1 / p_new[, 2], quantile(e, tau, names = FALSE))
#'   }
#'   cf <- crossfit_quantiles(d, "observed", fit_q, fit_r, K = 100, folds = 3,
#'                            seed = 1)
#'   k <- sharp_scan_mean_ci(cf$q, cf$r, d$observed, d$y, threshold = 0.5)
#'   print(k)
#'   plot(k)
#' }
sharp_scan_mean_ci <- function(q, r, response, outcome, target = "all",
                                lambda = seq(0, 2, length.out = 81),
                                weights = NULL, threshold = NULL, level = 0.9,
                                trim = 0, r_bounds = c(0.01, 0.99),
                                segment = NULL) {
  if (is.null(dim(q))) q <- matrix(q, ncol = 1L)
  q <- as.matrix(q)
  n <- nrow(q)
  if (!is.numeric(q) || anyNA(q) || any(!is.finite(q)))
    stop("`q` must be a finite numeric matrix.", call. = FALSE)
  if (!is.numeric(r) || length(r) != n)
    stop("`q` must have one row per value of `r`.", call. = FALSE)
  if (anyNA(r) || any(r <= 0 | r > 1))
    stop("`r` must lie in (0, 1], without NA.", call. = FALSE)
  if (is.logical(response)) response <- as.integer(response)
  if (!is.numeric(response) || length(response) != n || anyNA(response) ||
      !all(response %in% c(0, 1)))
    stop("`response` must be 0/1 or logical, one per row, without NA.",
         call. = FALSE)
  if (length(outcome) != n)
    stop("`outcome` must have one value per row.", call. = FALSE)
  if (!is.numeric(outcome))
    stop("`outcome` must be numeric. Convert a logical or factor outcome to ",
         "numbers first; for a binary outcome, sharp_scan_ci() is the direct ",
         "route.", call. = FALSE)
  R <- response
  if (!any(R == 1)) stop("No observed rows.", call. = FALSE)
  if (!is.numeric(outcome[R == 1]) || anyNA(outcome[R == 1]) ||
      any(!is.finite(outcome[R == 1])))
    stop("`outcome` must be finite where `response` is 1.", call. = FALSE)
  Y <- ifelse(R == 1, outcome, 0)
  if (!is.numeric(lambda) || length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda < 0))
    stop("`lambda` must be finite and non-negative.", call. = FALSE)
  lambda <- sort(unique(lambda))
  if (!is.null(threshold) && (!is.numeric(threshold) || length(threshold) != 1L ||
                              !is.finite(threshold)))
    stop("`threshold` must be a single finite number.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1)
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(trim) || length(trim) != 1L || is.na(trim) ||
      trim < 0 || trim >= 1)
    stop("`trim` must be a single number in [0, 1).", call. = FALSE)
  if (!is.numeric(r_bounds) || length(r_bounds) != 2L || anyNA(r_bounds) ||
      r_bounds[1] <= 0 || r_bounds[2] > 1 || r_bounds[1] >= r_bounds[2])
    stop("`r_bounds` must be two numbers with 0 < lower < upper <= 1.",
         call. = FALSE)
  w <- check_weights(weights, n)
  if (is.null(w)) w <- rep(1, n)

  # The target: the mean over all rows or over rows selected by X (weight
  # v = w s), or E(Y | R = 0) over the rows of a segment, whose sharp set is
  # {psi - E(RY)} / pr(R = 0) with psi the set of the mean over all rows of
  # the segment (mnar-sharp, Corollary 1).
  seg <- rep(1, n)
  if (is.character(target)) {
    target <- match.arg(target, c("all", "unobserved"))
    target_label <- target
  } else {
    if (!is.logical(target) || length(target) != n || anyNA(target))
      stop("`target` must be \"all\", \"unobserved\" or a logical vector, ",
           "one per row.", call. = FALSE)
    seg <- as.numeric(target); target_label <- "rows selected by X"
  }
  if (!is.null(segment)) {
    if (target_label == "rows selected by X")
      stop("`segment` combines with target = \"unobserved\" or \"all\"; for ",
           "a segment of all rows pass it as `target`.", call. = FALSE)
    if (!is.logical(segment) || length(segment) != n || anyNA(segment))
      stop("`segment` must be a logical vector chosen by X, one per row, ",
           "without NA.", call. = FALSE)
    seg <- as.numeric(segment)
    target_label <- paste(target_label, "rows of a segment chosen by X")
  }
  unobs <- startsWith(target_label, "unobserved")
  keep <- r >= trim
  in_target <- seg > 0 & (if (unobs) R == 0 else TRUE)
  n_trimmed <- sum(in_target & !keep)
  seg <- seg * keep; in_target <- in_target & keep
  ws <- w * seg
  den <- if (unobs) mean(ws * (1 - R)) else mean(ws)
  if (den <= 0) stop("The target is empty.", call. = FALSE)
  rc <- pmin(pmax(r, r_bounds[1]), r_bounds[2])
  n_clipped <- c(low = sum(seg > 0 & r < r_bounds[1]),
                 high = sum(seg > 0 & r > r_bounds[2]))
  z <- stats::qnorm(1 - (1 - level) / 2)
  prep <- mean_prep(q)
  mean_wRY <- mean(ws * R * Y)

  endpoint <- function(lam, side) {
    phi <- if (side == "upper") mean_score(prep$hi, prep$m, rc, Y, R, lam)$phi
           else -mean_score(prep$lo, -prep$m, rc, -Y, R, lam)$phi
    if (unobs) {
      est <- (mean(ws * phi) - mean_wRY) / den
      infl <- ws * (phi - R * Y - est * (1 - R)) / den
    } else {
      est <- mean(ws * phi) / den
      infl <- ws * (phi - est) / den
    }
    c(est = est, se = stats::sd(infl) / sqrt(n))
  }
  # plug-in outer band (without the cap pi <= 1) of the target
  naive <- function(lam, side) {
    o <- mean_bounds(prep, rc, lam, cap = FALSE)[[side]]
    if (unobs) (mean(ws * o) - mean_wRY) / den else mean(ws * o) / den
  }

  out <- data.frame(lambda = lambda, lower = NA_real_, upper = NA_real_,
                    se_lower = NA_real_, se_upper = NA_real_,
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    naive_lower = NA_real_, naive_upper = NA_real_)
  for (k in seq_along(lambda)) {
    lo <- endpoint(lambda[k], "lower"); up <- endpoint(lambda[k], "upper")
    out$lower[k] <- lo[["est"]]; out$se_lower[k] <- lo[["se"]]
    out$upper[k] <- up[["est"]]; out$se_upper[k] <- up[["se"]]
    out$naive_lower[k] <- naive(lambda[k], "lower")
    out$naive_upper[k] <- naive(lambda[k], "upper")
  }
  out$ci_lower <- out$lower - z * out$se_lower
  out$ci_upper <- out$upper + z * out$se_upper
  out <- monotone_curve(out)
  # The bounds lie in the range of the target rows' quantiles; far outside it
  # the inverse weights have blown up.
  rng <- range(q[seg > 0, , drop = FALSE])
  off <- max(rng[1] - out$lower, out$upper - rng[2], 0)
  if (off > 0.1 * max(diff(rng), .Machine$double.eps)) {
    obs <- R == 1 & seg > 0
    small_r <- any(obs & r < 0.05)
    outside <- mean(outcome[obs] < apply(q[obs, , drop = FALSE], 1, min) |
                    outcome[obs] > apply(q[obs, , drop = FALSE], 1, max))
    cause <- c(
      if (small_r) paste0("observed target rows with r-hat below 0.05 get large ",
                          "inverse weights 1/r-hat; check the response model for ",
                          "separation, or use `trim`"),
      if (outside > 0.1) sprintf(paste0("the quantiles do not cover the observed ",
                                        "outcomes (%.0f%% of them lie outside their ",
                                        "row's quantile range); check the quantile ",
                                        "model"), 100 * outside))
    if (is.null(cause))
      cause <- paste0("check the response model for separation and the quantile ",
                      "model for fit")
    warning(sprintf("The debiased estimate leaves the range of the quantiles by up to %.3g: %s.",
                    off, paste(cause, collapse = "; and ")), call. = FALSE)
  }

  point <- endpoint(0, "upper")[["est"]]
  lmax <- max(lambda)
  capped <- c(lower = sum(ws * mean_score(prep$lo, -prep$m, rc, -Y, R, lmax)$cap) / sum(ws),
              upper = sum(ws * mean_score(prep$hi, prep$m, rc, Y, R, lmax)$cap) / sum(ws))
  flip <- c(sharp = NA_real_, naive = NA_real_, ci = NA_real_)
  if (!is.null(threshold)) {
    side <- if (threshold >= point) "upper" else "lower"
    sgn <- if (side == "upper") 1 else -1
    est_fun <- function(l) endpoint(l, side)[["est"]]
    ci_fun <- function(l) { e <- endpoint(l, side); e[["est"]] + sgn * z * e[["se"]] }
    naive_fun <- function(l) naive(l, side)
    lam_max <- max(lambda)
    flip[["sharp"]] <- first_crossing(est_fun, threshold, lam_max, side,
                                      grid = out$lambda, values = out[[side]])
    flip[["naive"]] <- first_crossing(naive_fun, threshold, lam_max, side,
                                      grid = out$lambda, values = out[[paste0("naive_", side)]])
    flip[["ci"]] <- first_crossing(ci_fun, threshold, lam_max, side,
                                   grid = out$lambda, values = out[[paste0("ci_", side)]])
  }

  structure(list(
    curve = out, point = point, threshold = threshold, flip = flip,
    aggregate = "union", binding = NULL, level = level, target = target_label,
    trim = trim, trimmed = n_trimmed, clipped = n_clipped, r_bounds = r_bounds,
    unobserved = unobs,
    capped = capped,
    tail = c(target = sum(in_target & r < r_bounds[1]),
             observed = sum(R == 1 & r < r_bounds[1]),
             min_r_observed = min(r[R == 1])),
    n = n, n_target = sum(in_target),
    # for per-row plug-in bounds via as.data.frame(x, lambda = )
    q = q, r = rc, weights = weights
  ), class = c("sharp_scan_mean_ci", "sharp_scan_ci", "sharp_scan_mean",
               "sharp_scan"))
}

#' @export
print.sharp_scan_mean_ci <- function(x, digits = 3, ...) {
  cv <- x$curve; last <- nrow(cv)
  f <- function(v) format(v, digits = digits)
  cat(sprintf("<sharp_scan_mean_ci> sharp set for a mean, debiased, with %g%% confidence set\n",
              100 * x$level))
  cat(sprintf("  rows: %d   target: %s (%d rows)\n", x$n, x$target, x$n_target))
  if (x$trim > 0)
    cat(sprintf("  trimmed: %d rows with r-hat < %g, out of the target\n", x$trimmed, x$trim))
  if (sum(x$clipped) > 0) {
    part <- c(if (x$clipped[["low"]] > 0)
                sprintf("%d rows entering the set with r-hat below %g (raised to it)", x$clipped[["low"]], x$r_bounds[1]),
              if (x$clipped[["high"]] > 0)
                sprintf("%d rows entering the set with r-hat above %g (lowered to it)", x$clipped[["high"]], x$r_bounds[2]))
    cat(strwrap(paste0("clipped: ", paste(part, collapse = ", ")), indent = 2, exdent = 4), sep = "\n")
  }
  cat(sprintf("  r-hat < %g: %d target rows, %d observed rows; smallest r-hat observed %s\n",
              x$r_bounds[1], as.integer(x$tail[["target"]]), as.integer(x$tail[["observed"]]),
              format(x$tail[["min_r_observed"]], digits = 2)))
  i0 <- which(cv$lambda == 0)
  if (length(i0))
    cat(sprintf("  lambda = 0: %s [%s, %s]\n", f(x$point), f(cv$ci_lower[i0]), f(cv$ci_upper[i0])))
  cat(sprintf("  lambda = %s (propensity ratio %s): set [%s, %s], confidence set [%s, %s]\n",
              f(cv$lambda[last]), f(exp(cv$lambda[last])), f(cv$lower[last]),
              f(cv$upper[last]), f(cv$ci_lower[last]), f(cv$ci_upper[last])))
  if (!is.null(x$threshold)) {
    fl <- function(v) if (is.na(v)) "not within the scanned range" else
      sprintf("lambda = %s (propensity ratio %s)", f(v), f(exp(v)))
    cat("  threshold", f(x$threshold), "reached by\n")
    cat("    estimated set:       ", fl(x$flip[["sharp"]]), "\n")
    cat("    confidence set:      ", fl(x$flip[["ci"]]), "\n")
    cat("    outer band (plug-in):", fl(x$flip[["naive"]]), "\n")
  }
  cap_lines(x, digits)
  invisible(x)
}
