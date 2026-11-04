#' Confidence set for the binary scan, by debiased estimation
#'
#' `shadow_scan_ci()` estimates the sharp set of [shadow_scan()] with its
#' sampling error and returns, at each \eqn{\lambda}, a confidence set that
#' covers the identified set. It follows Semenova (2025): each endpoint is a
#' mean over rows of the minimum (or maximum) of two pieces, the band and the
#' worst-case bound; the piece that binds in a row is chosen from nuisance
#' estimates fitted without that row, and what is averaged is an orthogonal
#' signal of the chosen piece, not its plug-in value. This removes the bias of
#' the plug-in, which takes in every row the smaller of two noisy pieces and
#' so lies inside the set, most of all where many rows are near the switch
#' between band and bound. The confidence set uses the normal approximation
#' that Semenova's theory gives.
#'
#' @section Inputs:
#' Unlike [shadow_scan()], the inputs are given for **all** rows, observed and
#' unobserved, together with the response indicator and the outcome:
#' the signals use \eqn{RY} and \eqn{R(Y - \hat\eta_{\mathrm{obs}})/\hat r}.
#' `eta_obs` and `r` must be **cross-fitted**: the prediction for a row comes
#' from models fitted without it. [crossfit_inputs()] does that for any
#' learner. In-sample predictions make the confidence set too narrow.
#'
#' `target` says which rows the portfolio is: `"unobserved"` (the rejected
#' applications, the usual case), `"all"`, or a logical vector of rows chosen
#' by \eqn{X} alone (a segment). For `"unobserved"`, the share of a row in the
#' target, \eqn{1 - r(x)}, is itself estimated, and the score accounts for it.
#'
#' @section What is covered:
#' At each \eqn{\lambda} separately, the interval
#' \eqn{[\hat\psi^- - z\,\widehat{se}^-,\ \hat\psi^+ + z\,\widehat{se}^+]} with
#' \eqn{z} the \eqn{1 - (1 - \mathrm{level})/2} normal quantile covers the
#' identified set \eqn{[\psi^-, \psi^+]} with probability at least `level`
#' asymptotically. It is a confidence set for the set, not for the parameter
#' inside it, and it holds pointwise in \eqn{\lambda}, not uniformly over the
#' curve. The theory requires that rows where band and bound are exactly
#' equal carry no probability mass, which holds when \eqn{X} has a continuous
#' component, and overlap: \eqn{r(x)} bounded away from zero.
#'
#' Overlap is the binding condition in practice. The correction term weights
#' each observed row by \eqn{1/\hat r}, so a few observed rows in regions
#' where almost nothing is observed can carry most of the variance of the
#' score. The usual standard error, the sample standard deviation of the
#' influence values (`se_type = "empirical"`), then depends on whether the
#' sample happens to contain such rows, and it is too small on average: in
#' simulations for the rejected rows it covered 0.84 to 0.88 at nominal 0.90,
#' without improving from 1,000 to 16,000 rows. The default
#' `se_type = "conditional"` computes the variance of the score given \eqn{X}
#' in closed form from the nuisances, for instance
#' \eqn{(1-\hat r)^2 f'(\hat\eta)^2 \hat\eta(1-\hat\eta)/\hat r} for the
#' band piece, and averages it over **all** rows, including the many
#' unobserved ones in the same regions. In the same simulations it covered
#' 0.91 to 0.97 from 4,000 rows on, with intervals about 2% wider; at 1,000
#' rows and small \eqn{\lambda}, 0.85. It relies on \eqn{\hat\eta_{\mathrm{obs}}}
#' for the conditional variance of \eqn{Y}.
#'
#' `trim` removes rows with \eqn{\hat r <} `trim` from the target, together
#' with their correction terms. That changes the target to the rows with
#' \eqn{\hat r \ge} `trim`, and the output says how many rows left it. Rows
#' with \eqn{\hat r} at zero, where the selection was deterministic, have to
#' leave the target in any case, as for [shadow_scan()].
#'
#' For the whole sample (`target = "all"`), where the weight does not depend
#' on \eqn{r}, both standard errors covered 0.88 to 0.92 from 1,000 rows on.
#' Percentile bootstrap intervals of the plug-in covered as little as 0.72 at
#' \eqn{\lambda = 2}.
#'
#' Only the union aggregation is covered: the common-tilt set is bounded by an
#' extreme over rows, for which this theory does not hold. For a non-binary
#' outcome, see [shadow_scan_mean()], for which no inference is offered.
#'
#' @param eta_obs Cross-fitted \eqn{\hat\eta_{\mathrm{obs}}(x_i)} for every row,
#'   in (0, 1).
#' @param r Cross-fitted \eqn{\hat r(x_i)} for every row, in (0, 1].
#' @param response Response indicator for every row, 0/1 or logical.
#' @param outcome Outcome, 0/1; ignored (and may be `NA`) where `response` is 0.
#' @param target `"unobserved"` (default), `"all"`, or a logical vector
#'   selecting the rows of the portfolio by \eqn{X}.
#' @param lambda,threshold As in [shadow_scan()].
#' @param weights Optional non-negative row weights, a function of \eqn{X}.
#' @param level Confidence level. Default 0.9.
#' @param trim Rows with \eqn{\hat r <} `trim` leave the target. Default 0.
#'   See the section on what is covered.
#' @param se_type `"conditional"` (default) or `"empirical"`; see the section
#'   on what is covered.
#' @return An object of class `c("shadow_scan_ci", "shadow_scan")`, as from
#'   [shadow_scan()], whose `curve` holds the debiased estimates `lower`,
#'   `upper`, their standard errors `se_lower`, `se_upper`, the confidence set
#'   `ci_lower`, `ci_upper`, and the plug-in naive band of the target rows.
#'   `flip` adds `ci`: the smallest \eqn{\lambda} at which the confidence set
#'   reaches the threshold, the point from which the data no longer rule out
#'   the other decision.
#' @references
#' Semenova, V. (2025). Debiased machine learning of aggregated intersection
#' bounds and other causal parameters. Working paper, SSRN 5134514.
#' \doi{10.2139/ssrn.5134514}
#' @seealso [crossfit_inputs()], [shadow_scan()]
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
#' s <- shadow_scan_ci(inp$eta_obs, inp$r, d$observed, d$y, threshold = 0.35)
#' s
#' plot(s)
shadow_scan_ci <- function(eta_obs, r, response, outcome,
                           target = "unobserved",
                           lambda = seq(0, 2, length.out = 81),
                           weights = NULL, threshold = NULL, level = 0.9,
                           trim = 0, se_type = c("conditional", "empirical")) {
  se_type <- match.arg(se_type)
  check_inputs(eta_obs, r)
  n <- length(eta_obs)
  if (is.logical(response)) response <- as.integer(response)
  if (!is.numeric(response) || length(response) != n || anyNA(response) ||
      !all(response %in% c(0, 1)))
    stop("`response` must be 0/1 or logical, one per row, without NA.",
         call. = FALSE)
  if (length(outcome) != n)
    stop("`outcome` must have one value per row.", call. = FALSE)
  y_obs <- outcome[response == 1]
  if (anyNA(y_obs) || !all(y_obs %in% c(0, 1)))
    stop("`outcome` must be 0/1 where `response` is 1.", call. = FALSE)
  if (!any(response == 1))
    stop("No observed rows.", call. = FALSE)
  RY <- ifelse(response == 1, outcome, 0)
  R <- response
  if (!is.numeric(lambda) || length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda < 0))
    stop("`lambda` must be finite and non-negative.", call. = FALSE)
  lambda <- sort(unique(lambda))
  w <- check_weights(weights, n)
  if (is.null(w)) w <- rep(1, n)
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L ||
       !is.finite(threshold) || threshold <= 0 || threshold >= 1))
    stop("`threshold` must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a single number in (0, 1).", call. = FALSE)

  # target membership s and its conditional mean p(x) = E[s | x]
  if (is.character(target)) {
    target <- match.arg(target, c("unobserved", "all"))
    if (target == "unobserved") { s <- 1 - R; p <- 1 - r }
    else { s <- rep(1, n); p <- s }
    target_label <- target
  } else {
    if (!is.logical(target) || length(target) != n || anyNA(target))
      stop("`target` must be \"unobserved\", \"all\" or a logical vector, ",
           "one per row.", call. = FALSE)
    s <- as.numeric(target); p <- s
    target_label <- "rows selected by X"
  }
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 1)
    stop("`trim` must be a single number in [0, 1).", call. = FALSE)
  # rows with r-hat below trim leave the target, with their correction terms
  keep <- r >= trim
  s <- s * keep; p <- p * keep
  if (sum(w * s) <= 0) stop("The target is empty.", call. = FALSE)
  z <- stats::qnorm(1 - (1 - level) / 2)

  ell <- stats::qlogis(eta_obs)
  m_lo <- r * eta_obs; m_up <- r * eta_obs + 1 - r
  # derivative of a band piece f(eta) = expit(logit eta +- lambda) in eta
  resid <- R * (outcome_or0(outcome, R) - eta_obs) / r
  band_signal <- function(f) f + f * (1 - f) / (eta_obs * (1 - eta_obs)) * resid
  den <- mean(w * s)

  # target membership: random (1 - R) for "unobserved", else a function of X
  s_random <- identical(target_label, "unobserved")
  eta_c <- eta_obs

  # debiased estimate and standard error of one endpoint
  endpoint <- function(lam, side) {
    if (side == "upper") {
      f <- stats::plogis(ell + lam)
      band <- f <= m_up
      h <- ifelse(band, f, m_up)
      rho <- ifelse(band, band_signal(f), RY + 1 - R)
    } else {
      f <- stats::plogis(ell - lam)
      band <- f >= m_lo
      h <- ifelse(band, f, m_lo)
      rho <- ifelse(band, band_signal(f), RY)
    }
    score <- w * (s * h + p * (rho - h))
    est <- mean(score) / den
    if (se_type == "empirical") {
      infl <- (score - est * w * s) / den
      return(c(est = est, se = stats::sd(infl) / sqrt(n)))
    }
    # Conditional variance given X of D = rho - h, and its covariance with the
    # random target indicator 1 - R, in closed form from the nuisances. They
    # are evaluated at every row, so the variance does not hinge on the few
    # observed rows with small r-hat that dominate the empirical one.
    fp <- f * (1 - f) / (eta_c * (1 - eta_c))
    if (side == "upper") {
      vD <- ifelse(band, fp^2 * eta_c * (1 - eta_c) / r,
                   r * (1 - eta_c) * (1 - r * (1 - eta_c)))
      cD <- ifelse(band, 0, (1 - r) * r * (1 - eta_c))
    } else {
      vD <- ifelse(band, fp^2 * eta_c * (1 - eta_c) / r,
                   r * eta_c * (1 - r * eta_c))
      cD <- ifelse(band, 0, -(1 - r) * r * eta_c)
    }
    a <- h - est
    if (s_random) {                 # IF * den = (1 - R) a + (1 - r) D, times w
      # p = (1 - r) * keep; trimmed rows contribute nothing
      m1 <- w * p * a
      v1 <- w^2 * keep * (a^2 * r * (1 - r) + p^2 * vD + 2 * a * p * cD)
    } else {                        # IF * den = s (a + D), s a function of X
      m1 <- w * s * a
      v1 <- w^2 * s^2 * vD
    }
    v <- (mean((m1 - mean(m1))^2) + mean(v1)) / den^2
    c(est = est, se = sqrt(v / n))
  }

  K <- length(lambda)
  out <- data.frame(lambda = lambda, lower = NA_real_, upper = NA_real_,
                    se_lower = NA_real_, se_upper = NA_real_,
                    ci_lower = NA_real_, ci_upper = NA_real_,
                    naive_lower = NA_real_, naive_upper = NA_real_)
  tw <- w * s
  for (k in seq_len(K)) {
    lo <- endpoint(lambda[k], "lower"); up <- endpoint(lambda[k], "upper")
    out$lower[k] <- lo[["est"]]; out$se_lower[k] <- lo[["se"]]
    out$upper[k] <- up[["est"]]; out$se_upper[k] <- up[["se"]]
    out$naive_lower[k] <- sum(tw * stats::plogis(ell - lambda[k])) / sum(tw)
    out$naive_upper[k] <- sum(tw * stats::plogis(ell + lambda[k])) / sum(tw)
  }
  out$ci_lower <- pmax(out$lower - z * out$se_lower, 0)
  out$ci_upper <- pmin(out$upper + z * out$se_upper, 1)
  p0 <- endpoint(0, "upper")
  point <- p0[["est"]]
  flip <- c(sharp = NA_real_, naive = NA_real_, ci = NA_real_)
  if (!is.null(threshold)) {
    side <- if (threshold >= point) "upper" else "lower"
    sgn <- if (side == "upper") 1 else -1
    est_fun <- function(l) endpoint(l, side)[["est"]]
    ci_fun <- function(l) {
      e <- endpoint(l, side)
      e[["est"]] + sgn * z * e[["se"]]
    }
    naive_fun <- function(l)
      sum(tw * stats::plogis(ell + sgn * l)) / sum(tw)
    flip[["sharp"]] <- first_crossing(est_fun, threshold, max(lambda), side)
    flip[["naive"]] <- first_crossing(naive_fun, threshold, max(lambda), side)
    flip[["ci"]] <- first_crossing(ci_fun, threshold, max(lambda), side)
  }

  structure(list(
    curve = out, point = point, threshold = threshold, flip = flip,
    aggregate = "union", binding = NULL, level = level, target = target_label,
    trim = trim, trimmed = sum(!keep), se_type = se_type,
    eta_obs = eta_obs, r = r, weights = weights, n = n,
    n_target = sum(s > 0)
  ), class = c("shadow_scan_ci", "shadow_scan"))
}

outcome_or0 <- function(outcome, R) ifelse(R == 1, outcome, 0)

#' @export
print.shadow_scan_ci <- function(x, digits = 3, ...) {
  cv <- x$curve
  last <- nrow(cv)
  cat("<shadow_scan_ci> sharp set for a rate, debiased, with",
      sprintf("%g%% confidence set\n", 100 * x$level))
  cat("  rows:", x$n, "  target:", x$target, sprintf("(%d rows)\n", x$n_target))
  if (x$trim > 0)
    cat(sprintf("  trimmed: %d rows with r-hat < %g, out of the target\n",
                x$trimmed, x$trim))

  cat("  lambda = 0:", format(x$point, digits = digits),
      sprintf("[%s, %s]\n", format(cv$ci_lower[1], digits = digits),
              format(cv$ci_upper[1], digits = digits)))
  cat(sprintf("  lambda = %s (propensity ratio %s): set [%s, %s], confidence set [%s, %s]\n",
              format(cv$lambda[last], digits = digits),
              format(exp(cv$lambda[last]), digits = digits),
              format(cv$lower[last], digits = digits),
              format(cv$upper[last], digits = digits),
              format(cv$ci_lower[last], digits = digits),
              format(cv$ci_upper[last], digits = digits)))
  if (!is.null(x$threshold)) {
    fl <- function(v) if (is.na(v)) "not within the scanned range" else
      sprintf("lambda = %s (propensity ratio %s)", format(v, digits = digits),
              format(exp(v), digits = digits))
    cat("  threshold", format(x$threshold, digits = digits), "reached by\n")
    cat("    estimated set:  ", fl(x$flip[["sharp"]]), "\n")
    cat("    confidence set: ", fl(x$flip[["ci"]]), "\n")
    cat("    naive band:     ", fl(x$flip[["naive"]]), "\n")
  }
  invisible(x)
}

#' Cross-fitted scan inputs for any learner
#'
#' Splits the rows into folds and predicts \eqn{\hat\eta_{\mathrm{obs}}} and
#' \eqn{\hat r} for each fold from models fitted on the other folds, as
#' [shadow_scan_ci()] requires. The learner is yours: `fit_eta` and `fit_r`
#' are functions that fit on a training frame and return probabilities for
#' new rows. `fit_eta` receives only the observed rows of the training folds.
#'
#' @param data A data frame with all rows, observed and unobserved.
#' @param response Name of the response indicator column (0/1 or logical).
#' @param fit_eta `function(train, newdata)` returning
#'   \eqn{P(Y = 1 \mid X, R = 1)} for `newdata`; `train` holds observed rows.
#' @param fit_r `function(train, newdata)` returning \eqn{P(R = 1 \mid X)}.
#' @param folds Number of folds. Default 5.
#' @param seed Optional seed for the fold assignment.
#' @param clip Predictions are clipped to `[clip, 1 - clip]`. Default 1e-6.
#' @return A list with `eta_obs`, `r` (one per row of `data`) and `fold`.
#' @seealso [shadow_scan_ci()]
#' @export
crossfit_inputs <- function(data, response, fit_eta, fit_r, folds = 5L,
                            seed = NULL, clip = 1e-6) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!response %in% names(data))
    stop("Column `", response, "` not found.", call. = FALSE)
  if (!is.function(fit_eta) || !is.function(fit_r))
    stop("`fit_eta` and `fit_r` must be functions(train, newdata).", call. = FALSE)
  n <- nrow(data)
  folds <- as.integer(folds)
  if (folds < 2L || folds > n) stop("`folds` must be between 2 and the number of rows.",
                                    call. = FALSE)
  R <- as.integer(as.logical(data[[response]]))
  draw <- function() sample(rep_len(seq_len(folds), n))
  fold <- if (is.null(seed)) draw() else with_seed(seed, draw())
  eta <- r <- rep(NA_real_, n)
  for (k in seq_len(folds)) {
    tr <- data[fold != k, , drop = FALSE]
    te <- data[fold == k, , drop = FALSE]
    pe <- as.numeric(fit_eta(tr[R[fold != k] == 1, , drop = FALSE], te))
    pr <- as.numeric(fit_r(tr, te))
    if (length(pe) != nrow(te) || length(pr) != nrow(te))
      stop("A fit function returned too few or too many predictions: it must ",
           "return one per row of `newdata`.", call. = FALSE)
    eta[fold == k] <- pe; r[fold == k] <- pr
  }
  if (anyNA(eta) || anyNA(r))
    stop("A fit function returned missing or too few predictions.", call. = FALSE)
  list(eta_obs = pmin(pmax(eta, clip), 1 - clip),
       r = pmin(pmax(r, clip), 1), fold = fold)
}
