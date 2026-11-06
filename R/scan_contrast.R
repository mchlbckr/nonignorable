#' Debiased scan for the difference between two portfolios
#'
#' `shadow_scan_contrast()` compares two disjoint portfolios chosen by
#' \eqn{X} -- treatment arms of a trial, or two segments such as two regions --
#' under the same assumption as [shadow_scan()]: in every row the ratio of the
#' response propensities is at most \eqn{e^\lambda}. It answers the question
#' a comparison is run for: up to which strength of selection does the
#' finding "A differs from B" survive?
#'
#' Under that assumption each row may take its own extreme, so the two
#' portfolios may sit at opposite ends at once, and the sharp set for the
#' difference of rates \eqn{A - B} is
#' \eqn{[\mathrm{lower}_A - \mathrm{upper}_B,\ \mathrm{upper}_A - \mathrm{lower}_B]}.
#' Each endpoint is estimated by [shadow_scan_ci()]. Because the two targets
#' are disjoint, the influence functions of the two endpoints are non-zero on
#' different rows and have mean zero, so they are uncorrelated and the
#' standard error of the difference is the root of the sum of the squared
#' standard errors. Overlapping targets would need their covariance and are
#' not supported. By default each target is a whole segment of rows, observed
#' and unobserved alike; `among = "unobserved"` compares the unobserved rows
#' of the two segments instead, such as the rejected applicants of each.
#'
#' The flip points say from which \eqn{\lambda} on the set, and the confidence
#' set, reach `threshold` (by default 0: no difference). They are found by
#' linear interpolation on the `lambda` grid, so a fine grid gives a precise
#' flip.
#'
#' @inheritParams shadow_scan_ci
#' @param target_a,target_b Logical vectors, one per row, selecting the two
#'   portfolios by \eqn{X}. They must be disjoint and non-empty.
#' @param lambda Grid of \eqn{\lambda} values; must include 0.
#' @param threshold The difference at which the decision changes. Default 0.
#' @param labels Names of the two portfolios, for printing and plotting.
#' @param among `"all"` (default): each target is a whole segment, observed
#'   and unobserved rows alike. `"unobserved"`: only the unobserved rows of
#'   each segment, such as the rejected applications of two regions; see the
#'   `segment` argument of [shadow_scan_ci()]. The targets stay disjoint, so
#'   the standard errors still add in squares.
#' @return An object of class
#'   `c("shadow_scan_contrast", "shadow_scan_ci", "shadow_scan")` whose
#'   `curve` has the columns of [shadow_scan_ci()] for the difference
#'   \eqn{A - B}, with `point` the difference at \eqn{\lambda = 0} and `flip`
#'   the flip points of the set, the outer band and the confidence set.
#' @seealso [shadow_scan_ci()]
#' @export
#' @examples
#' d <- simulate_shadow(3000, tilt = -1, seed = 3)
#' f <- ~ x1 + x2 + z_shadow + z_policy + z_score
#' inp <- crossfit_inputs(d, "observed",
#'   fit_eta = function(train, new) predict(glm(update(f, y ~ .), binomial, train),
#'                                          new, type = "response"),
#'   fit_r   = function(train, new) predict(glm(update(f, observed ~ .), binomial, train),
#'                                          new, type = "response"),
#'   seed = 1)
#' k <- shadow_scan_contrast(inp$eta_obs, inp$r, d$observed, d$y,
#'                           target_a = d$x1 > 0, target_b = d$x1 <= 0,
#'                           labels = c("x1 > 0", "x1 <= 0"))
#' k
#' plot(k)
shadow_scan_contrast <- function(eta_obs, r, response, outcome,
                                 target_a, target_b,
                                 lambda = seq(0, 2, length.out = 81),
                                 weights = NULL, threshold = 0, level = 0.9,
                                 trim = 0,
                                 se_type = c("max", "conditional", "empirical"),
                                 labels = c("A", "B"), among = c("all", "unobserved")) {
  se_type <- match.arg(se_type)
  among <- match.arg(among)
  n <- length(eta_obs)
  check_contrast_args(n, target_a, target_b, lambda, labels)
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) ||
      abs(threshold) >= 1)
    stop("`threshold` must be a single number in (-1, 1).", call. = FALSE)
  arm <- function(tg, label) relabel_warnings(label,
    if (among == "all")
      shadow_scan_ci(eta_obs, r, response, outcome, target = tg,
                     lambda = lambda, weights = weights,
                     level = level, trim = trim, se_type = se_type)
    else
      shadow_scan_ci(eta_obs, r, response, outcome, target = "unobserved",
                     segment = tg, lambda = lambda, weights = weights,
                     level = level, trim = trim, se_type = se_type))
  a <- arm(target_a, labels[1]); b <- arm(target_b, labels[2])
  k <- contrast_from_arms(a, b, threshold, level, labels, n, range = c(-1, 1))
  k$se_type <- se_type
  k$among <- among
  k
}

#' Debiased scan for the difference of two portfolio means
#'
#' `shadow_scan_mean_contrast()` is [shadow_scan_contrast()] for an outcome
#' that is not binary: it compares the means of two disjoint portfolios chosen
#' by \eqn{X} under the assumption of [shadow_scan_mean()], with each
#' endpoint estimated by [shadow_scan_mean_ci()]. The sharp set for
#' \eqn{A - B} is \eqn{[\mathrm{lower}_A - \mathrm{upper}_B,\
#' \mathrm{upper}_A - \mathrm{lower}_B]}; because the targets are disjoint,
#' the standard error of each endpoint is the root of the sum of the squared
#' standard errors of the arms. Each target is a whole segment of rows,
#' observed and unobserved alike.
#'
#' @inheritParams shadow_scan_mean_ci
#' @param target_a,target_b Logical vectors, one per row, selecting the two
#'   portfolios by \eqn{X}. They must be disjoint and non-empty.
#' @param lambda Grid of \eqn{\lambda} values; must include 0.
#' @param threshold The difference at which the decision changes. Default 0.
#' @param labels Names of the two portfolios, for printing and plotting.
#' @return An object of class
#'   `c("shadow_scan_contrast", "shadow_scan_ci", "shadow_scan")` as from
#'   [shadow_scan_contrast()], with `of_means = TRUE`.
#' @seealso [shadow_scan_contrast()], [shadow_scan_mean_ci()],
#'   [crossfit_quantiles()]
#' @export
#' @examples
#' d <- simulate_shadow(3000, tilt = -1, outcome = "continuous", seed = 3)
#' q <- obs_quantiles(d, K = 50)   # the true quantile grids of the simulation
#' k <- shadow_scan_mean_contrast(q, d$obs_rate, d$observed, d$y,
#'                                target_a = d$x1 > 0, target_b = d$x1 <= 0,
#'                                lambda = seq(0, 2, length.out = 21),
#'                                labels = c("x1 > 0", "x1 <= 0"))
#' k
#' plot(k)
shadow_scan_mean_contrast <- function(q, r, response, outcome, target_a, target_b,
                                      lambda = seq(0, 2, length.out = 81),
                                      weights = NULL, threshold = 0, level = 0.9,
                                      trim = 0, r_bounds = c(0.01, 0.99),
                                      labels = c("A", "B")) {
  n <- if (is.null(dim(q))) length(q) else nrow(q)
  check_contrast_args(n, target_a, target_b, lambda, labels)
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("`threshold` must be a single finite number.", call. = FALSE)
  arm <- function(tg, label) relabel_warnings(label,
    shadow_scan_mean_ci(q, r, response, outcome, target = tg, lambda = lambda,
                        weights = weights, level = level, trim = trim,
                        r_bounds = r_bounds))
  a <- arm(target_a, labels[1]); b <- arm(target_b, labels[2])
  k <- contrast_from_arms(a, b, threshold, level, labels, n, range = c(-Inf, Inf))
  k$of_means <- TRUE
  k$clipped <- a$clipped + b$clipped
  k$r_bounds <- r_bounds
  k
}

check_contrast_args <- function(n, target_a, target_b, lambda, labels) {
  for (tg in list(target_a, target_b))
    if (!is.logical(tg) || length(tg) != n || anyNA(tg))
      stop("`target_a` and `target_b` must be logical vectors, one per row, ",
           "without NA.", call. = FALSE)
  if (!any(target_a) || !any(target_b))
    stop("`target_a` and `target_b` must not be empty.", call. = FALSE)
  if (any(target_a & target_b))
    stop("`target_a` and `target_b` must be disjoint.", call. = FALSE)
  if (!is.numeric(lambda) || !any(lambda == 0))
    stop("`lambda` must include lambda = 0.", call. = FALSE)
  if (length(unique(lambda)) < 2L)
    stop("`lambda` needs at least two distinct values, so that flip points can ",
         "be interpolated.", call. = FALSE)
  if (!is.character(labels) || length(labels) != 2L)
    stop("`labels` must be two strings.", call. = FALSE)
}

relabel_warnings <- function(label, expr) withCallingHandlers(expr,
  warning = function(cnd) {
    warning(label, ": ", conditionMessage(cnd), call. = FALSE)
    invokeRestart("muffleWarning")
  })

# The difference A - B from two arm scans on disjoint targets. `range` clamps
# the confidence set (a difference of rates lies in [-1, 1]).
contrast_from_arms <- function(a, b, threshold, level, labels, n, range) {
  ca <- a$curve; cb <- b$curve
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(lambda = ca$lambda,
                    lower = ca$lower - cb$upper, upper = ca$upper - cb$lower,
                    se_lower = sqrt(ca$se_lower^2 + cb$se_upper^2),
                    se_upper = sqrt(ca$se_upper^2 + cb$se_lower^2))
  out$ci_lower <- pmax(out$lower - z * out$se_lower, range[1])
  out$ci_upper <- pmin(out$upper + z * out$se_upper, range[2])
  out$naive_lower <- ca$naive_lower - cb$naive_upper
  out$naive_upper <- ca$naive_upper - cb$naive_lower
  out <- monotone_curve(out)
  point <- out$upper[out$lambda == 0]
  side <- if (threshold >= point) "upper" else "lower"
  lam_max <- max(out$lambda)
  on_grid <- function(v) stats::approxfun(out$lambda, v)
  pick <- function(lo, up) on_grid(if (side == "upper") up else lo)
  flip <- c(sharp = first_crossing(pick(out$lower, out$upper), threshold, lam_max, side, out$lambda),
            naive = first_crossing(pick(out$naive_lower, out$naive_upper), threshold, lam_max, side, out$lambda),
            ci = first_crossing(pick(out$ci_lower, out$ci_upper), threshold, lam_max, side, out$lambda))
  structure(list(
    curve = out, point = point, threshold = threshold, flip = flip,
    aggregate = "union", binding = NULL, level = level,
    trim = a$trim, trimmed = a$trimmed + b$trimmed, labels = labels,
    n = n, n_a = a$n_target, n_b = b$n_target, of_means = FALSE,
    tail = c(target = a$tail[["target"]] + b$tail[["target"]],
             observed = a$tail[["observed"]],
             min_r_observed = a$tail[["min_r_observed"]]),
    arms = list(a = a$curve, b = b$curve)
  ), class = c("shadow_scan_contrast", "shadow_scan_ci", "shadow_scan"))
}

#' @export
print.shadow_scan_contrast <- function(x, digits = 3, ...) {
  cv <- x$curve; last <- nrow(cv); i0 <- which(cv$lambda == 0)
  f <- function(v) format(v, digits = digits)
  cat(sprintf("<shadow_scan_contrast> difference of %s %s - %s, debiased, with %g%% confidence set\n",
              if (isTRUE(x$of_means)) "means" else "rates",
              x$labels[1], x$labels[2], 100 * x$level))
  cat(sprintf("  rows: %d   %s: %d %srows   %s: %d %srows\n", x$n, x$labels[1], x$n_a,
              if (identical(x$among, "unobserved")) "unobserved " else "",
              x$labels[2], x$n_b,
              if (identical(x$among, "unobserved")) "unobserved " else ""))
  if (x$trim > 0)
    cat(sprintf("  trimmed: %d rows with r-hat < %g, out of the targets\n", x$trimmed, x$trim))
  if (sum(x$clipped) > 0)
    cat(sprintf("  clipped: %d target rows with r-hat below %g, %d above %g\n",
                as.integer(x$clipped[["low"]]), x$r_bounds[1],
                as.integer(x$clipped[["high"]]), x$r_bounds[2]))
  cat(sprintf("  r-hat < 0.01: %d target rows, %d observed rows; smallest r-hat observed %s\n",
              as.integer(x$tail[["target"]]), as.integer(x$tail[["observed"]]),
              format(x$tail[["min_r_observed"]], digits = 2)))
  cat(sprintf("  lambda = 0: %s [%s, %s]\n", f(x$point), f(cv$ci_lower[i0]), f(cv$ci_upper[i0])))
  cat(sprintf("  lambda = %s (propensity ratio %s): set [%s, %s], confidence set [%s, %s]\n",
              f(cv$lambda[last]), f(exp(cv$lambda[last])), f(cv$lower[last]),
              f(cv$upper[last]), f(cv$ci_lower[last]), f(cv$ci_upper[last])))
  fl <- function(v) if (is.na(v)) "not within the scanned range" else
    sprintf("lambda = %s (propensity ratio %s)", f(v), f(exp(v)))
  cat(sprintf("  difference %s reached by\n", f(x$threshold)))
  cat("    estimated set:       ", fl(x$flip[["sharp"]]), "\n")
  cat("    confidence set:      ", fl(x$flip[["ci"]]), "\n")
  cat("    outer band (plug-in):", fl(x$flip[["naive"]]), "\n")
  invisible(x)
}
