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
#' `shadow_scan_mean_ci()` estimates the sharp identified set of
#' [shadow_scan_mean()] -- the mean of an outcome that may be missing not at
#' random, when the ratio of the response propensities of any two outcome
#' values is at most \eqn{e^\lambda} -- with a debiased estimator and a
#' confidence set. It is the counterpart of [shadow_scan_ci()] for an outcome
#' that is not binary.
#'
#' @section The estimator:
#' Each endpoint of the set is a mean over rows of \eqn{\mu^\pm(x)}, which
#' depends on \eqn{r(x)} and on the quantile function of \eqn{Y} given
#' \eqn{R = 1, X = x}. The plug-in, [shadow_scan_mean()], inherits the
#' first-order error of both estimates and is biased towards the inside of
#' the set. The estimator here adds to each row a correction term that
#' removes that first-order error: an inverse-probability weighted residual
#' of the observed outcome, and, where the constraint \eqn{\pi \le 1} binds,
#' a term in \eqn{R - \hat r}. The resulting score is Neyman-orthogonal in
#' both regimes, so the estimate is insensitive to small errors in the
#' nuisances, and the correction in \eqn{r} vanishes continuously where the
#' regime changes, so no margin condition is needed. The confidence set is
#' the Wald interval from the empirical standard deviation of the score.
#'
#' The quantiles enter as a grid of \eqn{K} levels \eqn{(k - 1/2)/K}, as in
#' [shadow_scan_mean()]. An estimator that uses the full empirical
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
#' the quantiles of the standardised residuals, is a good default. A linear
#' model for the location only cannot follow heteroskedasticity and makes the
#' set markedly too wide; quantile forests are valid but conservative.
#'
#' The target is the mean over all rows (`target = "all"`) or over a segment
#' chosen by \eqn{X} (a logical vector). The mean over the unobserved rows
#' alone is not covered by this theory.
#'
#' @param q An \eqn{n \times K} matrix of conditional quantiles of \eqn{Y}
#'   given \eqn{R = 1, X = x_i} for **every** row, fitted out of fold, at the
#'   levels \eqn{(k - 1/2)/K}.
#' @param r Cross-fitted \eqn{\hat r(x_i)} for every row, in (0, 1].
#' @param response Response indicator for every row, 0/1 or logical.
#' @param outcome Outcome; ignored (and may be `NA`) where `response` is 0.
#' @param target `"all"` (default) or a logical vector selecting the rows of
#'   the portfolio by \eqn{X}.
#' @param lambda Grid of \eqn{\lambda} values.
#' @param weights Optional non-negative row weights, a function of \eqn{X}.
#' @param threshold Optional value of the mean at which the decision changes.
#' @param level Confidence level. Default 0.9.
#' @param trim Rows with \eqn{\hat r <} `trim` leave the target. Default 0.
#' @param r_bounds \eqn{\hat r} is clipped into this range before use.
#'   Default `c(0.01, 0.99)`. The number of clipped rows is reported.
#' @return An object of class
#'   `c("shadow_scan_mean_ci", "shadow_scan_ci", "shadow_scan_mean",
#'   "shadow_scan")` with `curve` (columns `lambda`, `lower`, `upper`,
#'   `se_lower`, `se_upper`, `ci_lower`, `ci_upper`, `naive_lower`,
#'   `naive_upper`; the naive columns are the plug-in outer band), `point`
#'   (the estimate at \eqn{\lambda = 0}), `flip` (`sharp`, `ci`, `naive`),
#'   `trimmed`, `clipped`, `tail`, `n` and `n_target`. `as.data.frame(x,
#'   lambda = )` gives the plug-in bounds per row, with \eqn{\hat r} clipped
#'   into `r_bounds`; they carry no correction and no confidence statement.
#' @seealso [shadow_scan_mean()], [crossfit_quantiles()], [shadow_scan_ci()]
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
#'   k <- shadow_scan_mean_ci(cf$q, cf$r, d$observed, d$y, threshold = 0.5)
#'   print(k)
#'   plot(k)
#' }
shadow_scan_mean_ci <- function(q, r, response, outcome, target = "all",
                                lambda = seq(0, 2, length.out = 81),
                                weights = NULL, threshold = NULL, level = 0.9,
                                trim = 0, r_bounds = c(0.01, 0.99)) {
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
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 1)
    stop("`trim` must be a single number in [0, 1).", call. = FALSE)
  if (!is.numeric(r_bounds) || length(r_bounds) != 2L || anyNA(r_bounds) ||
      r_bounds[1] <= 0 || r_bounds[2] > 1 || r_bounds[1] >= r_bounds[2])
    stop("`r_bounds` must be two numbers with 0 < lower < upper <= 1.",
         call. = FALSE)
  w <- check_weights(weights, n)
  if (is.null(w)) w <- rep(1, n)

  if (is.character(target)) {
    if (!identical(target, "all"))
      stop("Only target = \"all\" or a logical vector chosen by X is covered; ",
           "the mean over the unobserved rows alone is not covered by this ",
           "theory. See ?shadow_scan_mean_ci.", call. = FALSE)
    s <- rep(1, n); target_label <- "all"
  } else {
    if (!is.logical(target) || length(target) != n || anyNA(target))
      stop("`target` must be \"all\" or a logical vector, one per row.",
           call. = FALSE)
    s <- as.numeric(target); target_label <- "rows selected by X"
  }
  keep <- r >= trim
  n_trimmed <- sum(s > 0 & !keep)
  s <- s * keep
  if (sum(w * s) <= 0) stop("The target is empty.", call. = FALSE)
  rc <- pmin(pmax(r, r_bounds[1]), r_bounds[2])
  n_clipped <- sum(rc != r)
  z <- stats::qnorm(1 - (1 - level) / 2)
  prep <- mean_prep(q)
  ws <- w * s
  den <- mean(ws)

  endpoint <- function(lam, side) {
    phi <- if (side == "upper") mean_score(prep$hi, prep$m, rc, Y, R, lam)$phi
           else -mean_score(prep$lo, -prep$m, rc, -Y, R, lam)$phi
    est <- mean(ws * phi) / den
    infl <- (ws * phi - est * ws) / den
    c(est = est, se = stats::sd(infl) / sqrt(n))
  }

  out <- data.frame(lambda = lambda, lower = NA_real_, upper = NA_real_,
                    se_lower = NA_real_, se_upper = NA_real_,
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    naive_lower = NA_real_, naive_upper = NA_real_)
  for (k in seq_along(lambda)) {
    lo <- endpoint(lambda[k], "lower"); up <- endpoint(lambda[k], "upper")
    out$lower[k] <- lo[["est"]]; out$se_lower[k] <- lo[["se"]]
    out$upper[k] <- up[["est"]]; out$se_upper[k] <- up[["se"]]
    o <- mean_bounds(prep, rc, lambda[k], cap = FALSE)
    out$naive_lower[k] <- sum(ws * o$lower) / sum(ws)
    out$naive_upper[k] <- sum(ws * o$upper) / sum(ws)
  }
  out$ci_lower <- out$lower - z * out$se_lower
  out$ci_upper <- out$upper + z * out$se_upper
  # The bounds lie in the range of the target rows' quantiles; far outside it
  # the inverse weights have blown up.
  rng <- range(q[s > 0, , drop = FALSE])
  off <- max(rng[1] - out$lower, out$upper - rng[2], 0)
  if (off > 0.1 * max(diff(rng), .Machine$double.eps))
    warning(sprintf(paste0(
      "The debiased estimate leaves the range of the quantiles by up to %.3g: ",
      "the inverse weights 1/r-hat have blown up. Check the response model for ",
      "separation, or use `trim`."), off), call. = FALSE)

  point <- endpoint(0, "upper")[["est"]]
  flip <- c(sharp = NA_real_, naive = NA_real_, ci = NA_real_)
  if (!is.null(threshold)) {
    side <- if (threshold >= point) "upper" else "lower"
    sgn <- if (side == "upper") 1 else -1
    est_fun <- function(l) endpoint(l, side)[["est"]]
    ci_fun <- function(l) { e <- endpoint(l, side); e[["est"]] + sgn * z * e[["se"]] }
    naive_fun <- function(l) {
      o <- mean_bounds(prep, rc, l, cap = FALSE)
      sum(ws * o[[side]]) / sum(ws)
    }
    lam_max <- max(lambda)
    flip[["sharp"]] <- first_crossing(est_fun, threshold, lam_max, side)
    flip[["naive"]] <- first_crossing(naive_fun, threshold, lam_max, side)
    flip[["ci"]] <- first_crossing(ci_fun, threshold, lam_max, side)
  }

  structure(list(
    curve = out, point = point, threshold = threshold, flip = flip,
    aggregate = "union", binding = NULL, level = level, target = target_label,
    trim = trim, trimmed = n_trimmed, clipped = n_clipped, r_bounds = r_bounds,
    tail = c(target = sum(s > 0 & r < 0.01), observed = sum(R == 1 & r < 0.01),
             min_r_observed = min(r[R == 1])),
    n = n, n_target = sum(s > 0),
    # for per-row plug-in bounds via as.data.frame(x, lambda = )
    q = q, r = rc, weights = weights
  ), class = c("shadow_scan_mean_ci", "shadow_scan_ci", "shadow_scan_mean",
               "shadow_scan"))
}

#' @export
print.shadow_scan_mean_ci <- function(x, digits = 3, ...) {
  cv <- x$curve; last <- nrow(cv)
  f <- function(v) format(v, digits = digits)
  cat(sprintf("<shadow_scan_mean_ci> sharp set for a mean, debiased, with %g%% confidence set\n",
              100 * x$level))
  cat(sprintf("  rows: %d   target: %s (%d rows)\n", x$n, x$target, x$n_target))
  if (x$trim > 0)
    cat(sprintf("  trimmed: %d rows with r-hat < %g, out of the target\n", x$trimmed, x$trim))
  if (x$clipped > 0)
    cat(sprintf("  clipped: %d rows with r-hat outside [%g, %g]\n", x$clipped,
                x$r_bounds[1], x$r_bounds[2]))
  cat(sprintf("  r-hat < 0.01: %d target rows, %d observed rows; smallest r-hat observed %s\n",
              as.integer(x$tail[["target"]]), as.integer(x$tail[["observed"]]),
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
    cat("    estimated set:  ", fl(x$flip[["sharp"]]), "\n")
    cat("    confidence set: ", fl(x$flip[["ci"]]), "\n")
    cat("    outer band:     ", fl(x$flip[["naive"]]), "\n")
  }
  invisible(x)
}
