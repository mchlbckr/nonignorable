#' Sharp identified set for a binary outcome under a bounded ratio of response propensities
#'
#' `sharp_scan()` answers *how far the answer can move* if the outcome is
#' missing not at random. It assumes nothing but a bound \eqn{|b(x)| \le
#' \lambda} on the log ratio of the response propensities of the two outcome
#' values and needs **no shadow variable**. It is a different tool from [shadow_screen()], which asks
#' whether a variable is fit to serve as a shadow variable for a
#' point-identified model; the two rest on different assumptions.
#'
#' @section The set:
#' Let \eqn{\eta(x) = P(Y=1 \mid X=x)} be the target,
#' \eqn{\eta_{\mathrm{obs}}(x) = P(Y=1 \mid X=x, R=1)} its observed-data
#' counterpart, \eqn{r(x) = P(R=1 \mid X=x)} the conditional observation rate
#' and \eqn{b(x) = \log\{\pi(1,x)/\pi(0,x)\}} with \eqn{\pi(y,x) = P(R=1 \mid
#' Y=y, X=x)}. By Bayes' rule \eqn{\mathrm{logit}\,\eta_{\mathrm{obs}} =
#' \mathrm{logit}\,\eta + b}.
#'
#' Note the scale: \eqn{b} is the log of a **risk ratio** of response
#' propensities. That it shifts the logit of \eqn{\eta} does not make it an
#' odds ratio. It is not the odds-ratio parameter of Vansteelandt et al.
#' (2006), whose set neither depends on \eqn{r} nor saturates, and not the
#' \eqn{\Lambda} of the marginal sensitivity model. The assumption is a member
#' of the class of Manski (2016, Sect. 4.2), with the ratio taken between
#' outcome values rather than relative to the response rate; for a
#' non-binary outcome see [sharp_scan_mean()]. With \eqn{\ell = \mathrm{logit}\,
#' \eta_{\mathrm{obs}}} and \eqn{\eta_t = \mathrm{expit}(\ell - t)}, the
#' identified set is
#' \deqn{\mathcal H_\lambda(x) = \{\eta_t : |t| \le \lambda,\; c_t \in (0,1],\;
#'   c_t e^{t} \le 1\}, \qquad c_t = \frac{r}{e^{t}\eta_t + 1 - \eta_t},}
#' where \eqn{c_t = \pi(0,x)} and \eqn{c_t e^t = \pi(1,x)} are the two
#' propensities that reproduce the observed \eqn{r}. Solving the two
#' constraints gives the closed form
#' \deqn{\mathcal H_\lambda(x) = \big[\max\{\mathrm{expit}(\ell-\lambda),\;
#'   r\eta_{\mathrm{obs}}\},\; \min\{\mathrm{expit}(\ell+\lambda),\;
#'   r\eta_{\mathrm{obs}} + 1 - r\}\big],}
#' the outer band \eqn{[\mathrm{expit}(\ell-\lambda),
#' \mathrm{expit}(\ell+\lambda)]} intersected with the worst-case (Manski)
#' bounds of Manski (1989). Both endpoints are attained, and at \eqn{\lambda = 0} the set is
#' \eqn{\{\eta_{\mathrm{obs}}\}}: the value any model fitted on the observed
#' rows asserts.
#'
#' @section When the sharp set is narrower, and when it is not:
#' It is never wider than the outer band. It is strictly narrower exactly where
#' a worst-case bound cuts the band: the lower end moves when \eqn{r >
#' \mathrm{expit}(\ell-\lambda)/\eta_{\mathrm{obs}}}, the upper end when
#' \eqn{r > \{1-\mathrm{expit}(\ell+\lambda)\}/(1-\eta_{\mathrm{obs}})}. For
#' \eqn{r \le e^{-\lambda}} neither moves. In reject inference the rows of
#' interest are the *rejected* applications, whose observation rate is
#' typically small, so the **upper** end there often coincides with the outer
#' one. Ignoring \eqn{r} is never wrong, only wasteful where \eqn{r} is
#' large: it gives an outer set, not a sharp one, as Dorn and Guo (2023)
#' observe for the marginal sensitivity model.
#'
#' @section The two inputs:
#' The scan takes two numeric vectors, not one, evaluated at the rows for which
#' the set is wanted: \eqn{\hat\eta_{\mathrm{obs}}(x)} and \eqn{\hat r(x)}.
#' `r` is the *conditional observation rate* \eqn{P(R=1 \mid X=x)}, not the row
#' indicator and not the global response rate. For rejected applications
#' \eqn{\hat\eta_{\mathrm{obs}}} is an extrapolation of a model fitted on the
#' accepted ones; that is a property of the task, not of this package. The
#' package never sees the learner: `glm`, a scorecard, a random forest or a
#' model from another language serve equally. [scan_inputs()] extracts both
#' vectors from fitted models with a `predict()` method. The set requires
#' \eqn{\hat r(x) > 0}: where acceptance was deterministic, nothing is
#' observed and \eqn{\eta_{\mathrm{obs}}} is undefined.
#'
#' @section Estimated inputs, no inference:
#' The set is an *identified* set: it is what the data would permit if the
#' observed-data law were known exactly. The scan plugs in
#' \eqn{\hat\eta_{\mathrm{obs}}} and \eqn{\hat r} as if they were the truth.
#' Their estimation error is **not** included, so the reported set is
#' narrower than an honest statement about the population would be. For a
#' confidence set use [sharp_scan_ci()], which needs cross-fitted inputs for
#' all rows and the observed outcomes. Inference for such sets is not routine: the
#' endpoints are maxima and minima of estimated quantities, where the usual
#' bootstrap can fail; a flexible learner needs more than a naive plug-in;
#' and covering the set differs from covering the parameter inside it. Nor
#' does resampling capture misspecification of the extrapolation to rows that
#' are rarely observed.
#'
#' If you bootstrap the set yourself, re-estimating both inputs in every
#' draw, do not take percentile intervals of the endpoints. Each endpoint is a
#' mean over rows of a minimum (or maximum) of the band and the worst-case
#' bound; its plug-in is biased towards the inside of the set, and the
#' percentile bootstrap inherits that bias. In simulations at
#' \eqn{\lambda = 2}, nominal 90% percentile intervals covered 72%
#' (\eqn{n = 1000}) and 80% (\eqn{n = 4000}). The basic bootstrap,
#' \eqn{2\hat\theta - Q(\hat\theta^*)}, came close to nominal. The bias is
#' largest where many rows lie near the switch between band and bound, which
#' is where the scan looks for its flip.
#'
#' @section Aggregation:
#' The set is pointwise. For a portfolio rate (a weighted mean of
#' \eqn{\eta(x)}) there are two ways to aggregate:
#' * `"union"` (default): every row at its own extreme, averaged. Under the
#'   bound \eqn{|b(x)| \le \lambda} alone, with the tilt free to vary with
#'   \eqn{x}, this is the sharp set for the portfolio rate. Each row
#'   contributes only its own bounds, so no single row dominates.
#' * `"common"`: one tilt \eqn{t} for all rows -- "had selection the strength
#'   \eqn{t} everywhere, the rate would be \eqn{R(t)}". This adds an
#'   assumption, a tilt constant in \eqn{x}, and gives a narrower set. Then
#'   \eqn{t} must be admissible in *every* row; the admissible range is the
#'   intersection of the row ranges \eqn{[\ell - \mathrm{logit}\,U,\; \ell -
#'   \mathrm{logit}\,L]}. It always contains \eqn{t = 0}. A single row with
#'   \eqn{\hat r} near one can narrow it a lot; the rows that bind are
#'   reported in `$binding`. Because the range is set by the most extreme
#'   rows, with *estimated* inputs it inherits their estimation error, at
#'   exactly the rows where \eqn{\hat r} is least reliable. Inspect the
#'   binding rows before reading a plateau as information. This is why it is
#'   not the default. It is a sensitivity calculation: no confidence set is
#'   offered for it, since its bound is an extremum over rows of an estimated
#'   function, set by the rows where the estimate is least reliable.
#'
#' @param eta_obs Numeric vector of \eqn{\hat\eta_{\mathrm{obs}}(x)}, in (0, 1).
#' @param r Numeric vector of \eqn{\hat r(x)}, in (0, 1], same length.
#' @param lambda Non-negative bounds on \eqn{|b|}. \eqn{e^\lambda} is the
#'   largest admissible ratio of the two response propensities (a risk ratio,
#'   not an odds ratio). Sorted and de-duplicated.
#' @param aggregate `"union"` (default) or `"common"`, see Aggregation.
#' @param weights Optional non-negative row weights for the portfolio mean.
#' @param threshold Optional decision threshold on the portfolio rate. If
#'   given, the smallest \eqn{\lambda} at which the sharp set (and the outer
#'   band) reaches it is computed exactly, by root finding on the closed form.
#'
#' @return An object of class `"sharp_scan"`, a list with
#'   * `curve`: data frame with `lambda`, `lower`, `upper` (sharp set of the
#'     portfolio rate) and `naive_lower`, `naive_upper`;
#'   * `point`: the portfolio rate at \eqn{\lambda = 0};
#'   * `threshold`, `flip`: the threshold and the \eqn{\lambda} at which the
#'     sharp set and the outer band first reach it (`NA` if not within the
#'     scanned range);
#'   * `binding`: for `"common"`, the row indices that bind the lower and the
#'     upper end at each \eqn{\lambda} (`NA` where the outer bound binds);
#'   * `capped`: for `"union"`, the weighted share of rows in which the cap
#'     \eqn{\pi \le 1} (the worst-case bound) binds at the largest
#'     \eqn{\lambda}, for the `lower` and the `upper` bound;
#'   * `held`: for `"union"`, the \eqn{\lambda} from which the cap binds in
#'     every row, so that this bound stays flat (`NA` if not within the
#'     scanned range). With \eqn{\hat r} near one this comes early: the lower
#'     bound \eqn{r\eta_{\mathrm{obs}}} binds from
#'     \eqn{\lambda = \mathrm{logit}\,\eta_{\mathrm{obs}} -
#'     \mathrm{logit}(r\eta_{\mathrm{obs}})}, about \eqn{-\log r} for rare
#'     events. That flat edge is the identified set, not a defect;
#'   * the inputs, so that `as.data.frame(x, lambda = )` can return the
#'     pointwise sets.
#' @references
#' Dorn, J. and Guo, K. (2023). Sharp sensitivity analysis for inverse
#' propensity weighting via quantile balancing. *Journal of the American
#' Statistical Association*, 118(544), 2645--2657.
#' \doi{10.1080/01621459.2022.2069572}
#'
#' Manski, C. F. (1989). Anatomy of the selection problem. *Journal of Human
#' Resources*, 24(3), 343--360. \doi{10.2307/145818}
#'
#' Manski, C. F. (2016). Credible interval estimates for official statistics
#' with survey nonresponse. *Journal of Econometrics*, 191(2), 293--301.
#' \doi{10.1016/j.jeconom.2015.12.002}
#'
#' Vansteelandt, S., Goetghebeur, E., Kenward, M. G. and Molenberghs, G.
#' (2006). Ignorance and uncertainty regions as inferential tools in a
#' sensitivity analysis. *Statistica Sinica*, 16(3), 953--979.
#' @seealso [plot.sharp_scan()], [sharp_scan_mean()] for a non-binary
#'   outcome, [shadow_screen()] for the other question.
#' @export
#' @examples
#' d <- simulate_shadow(5000, tilt = -1, seed = 1)
#' rej <- d[d$observed == 0, ]
#' s <- sharp_scan(rej$eta_obs, rej$obs_rate, lambda = seq(0, 2, by = 0.05),
#'                  threshold = 0.35)
#' s
#' mean(rej$eta)            # the truth, inside the set at lambda = 1
#' plot(s)
#' head(as.data.frame(s, lambda = 1))
sharp_scan <- function(eta_obs, r, lambda = seq(0, 2, length.out = 81),
                        aggregate = c("union", "common"), weights = NULL,
                        threshold = NULL) {
  aggregate <- match.arg(aggregate)
  check_inputs(eta_obs, r)
  if (!is.numeric(lambda) || length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda < 0))
    stop("`lambda` must be finite and non-negative.", call. = FALSE)
  lambda <- sort(unique(lambda))
  w <- check_weights(weights, length(eta_obs))
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L ||
       !is.finite(threshold) || threshold <= 0 || threshold >= 1))
    stop("`threshold` must be a single number in (0, 1).", call. = FALSE)

  ell <- stats::qlogis(eta_obs)
  # Manski bounds do not depend on lambda; their tilts give each row's
  # admissible t-range before the |t| <= lambda cap.
  m_lo <- r * eta_obs
  m_up <- r * eta_obs + 1 - r
  t_min_row <- ell - stats::qlogis(m_up)   # most negative admissible tilt
  t_max_row <- ell - stats::qlogis(m_lo)   # most positive admissible tilt

  K <- length(lambda)
  out <- data.frame(lambda = lambda, lower = NA_real_, upper = NA_real_,
                    naive_lower = NA_real_, naive_upper = NA_real_)
  binding <- data.frame(lambda = lambda, lower_row = NA_integer_,
                        upper_row = NA_integer_)
  # rows that can bind the common tilt: the most restrictive ones
  i_min <- which.max(t_min_row)
  i_max <- which.min(t_max_row)
  T_min <- t_min_row[i_min]
  T_max <- t_max_row[i_max]

  for (k in seq_len(K)) {
    lam <- lambda[k]
    out$naive_lower[k] <- wmean(stats::plogis(ell - lam), w)
    out$naive_upper[k] <- wmean(stats::plogis(ell + lam), w)
    if (aggregate == "common") {
      t_up <- max(-lam, T_min)
      t_lo <- min(lam, T_max)
      out$upper[k] <- wmean(stats::plogis(ell - t_up), w)
      out$lower[k] <- wmean(stats::plogis(ell - t_lo), w)
      if (T_min > -lam) binding$upper_row[k] <- i_min
      if (T_max < lam) binding$lower_row[k] <- i_max
    } else {
      out$lower[k] <- wmean(pmax(stats::plogis(ell - lam), m_lo), w)
      out$upper[k] <- wmean(pmin(stats::plogis(ell + lam), m_up), w)
    }
  }
  point <- wmean(eta_obs, w)
  # at lambda = 0 every expression above reduces to eta_obs; make it exact
  z <- lambda == 0
  out$lower[z] <- out$upper[z] <- out$naive_lower[z] <- out$naive_upper[z] <- point

  bound_fun <- function(side, sharp) {
    function(lam) {
      if (!sharp) {
        return(wmean(stats::plogis(ell + if (side == "upper") lam else -lam), w))
      }
      if (aggregate == "common") {
        t <- if (side == "upper") max(-lam, T_min) else min(lam, T_max)
        wmean(stats::plogis(ell - t), w)
      } else if (side == "upper") {
        wmean(pmin(stats::plogis(ell + lam), m_up), w)
      } else {
        wmean(pmax(stats::plogis(ell - lam), m_lo), w)
      }
    }
  }
  flip <- c(sharp = NA_real_, naive = NA_real_)
  if (!is.null(threshold)) {
    side <- if (threshold >= point) "upper" else "lower"
    for (nm in names(flip)) {
      flip[[nm]] <- first_crossing(bound_fun(side, nm == "sharp"), threshold,
                                   max(lambda), side)
    }
  }

  cap <- if (aggregate == "union") cap_summary(eta_obs, r, w, lambda) else NULL

  structure(list(
    curve = out, point = point, threshold = threshold, flip = flip,
    aggregate = aggregate, held = cap$held, capped = cap$capped,
    binding = if (aggregate == "common") binding else NULL,
    eta_obs = eta_obs, r = r, weights = weights, n = length(eta_obs)
  ), class = "sharp_scan")
}

# Where the cap pi <= 1 binds. Row i's lower bound is r * eta_obs, not the
# band, once lambda exceeds logit(eta_obs) - logit(r * eta_obs); its upper
# bound once lambda exceeds logit(r * eta_obs + 1 - r) - logit(eta_obs).
# `capped` is the weighted share of target rows (positive weight) where the
# cap binds at the largest lambda; `held` is the lambda from which it binds
# in every target row, so that the bound no longer moves (NA if not within
# the scanned range).
cap_summary <- function(eta_obs, r, w, lambda, target = NULL) {
  if (is.null(w)) w <- rep(1, length(eta_obs))
  if (!is.null(target)) w <- w * target
  keep <- w > 0
  ell <- stats::qlogis(eta_obs[keep]); rr <- r[keep]; e <- eta_obs[keep]; ww <- w[keep]
  from <- list(lower = ell - stats::qlogis(rr * e),
               upper = stats::qlogis(rr * e + 1 - rr) - ell)
  lam <- max(lambda)
  capped <- vapply(from, function(f) sum(ww[f < lam]) / sum(ww), 0)
  held <- vapply(from, function(f) { h <- max(max(f), 0); if (h <= lam) h else NA_real_ }, 0)
  list(capped = capped, held = held)
}

cap_lines <- function(x, digits = 3) {
  if (is.null(x$capped)) return(invisible())
  lam <- max(x$curve$lambda)
  pc <- function(v) sprintf("%.0f%%", 100 * v)
  cat(sprintf("  cap pi <= 1 at lambda = %s binds in %s of target rows for the lower bound, %s for the upper\n",
              format(lam, digits = digits), pc(x$capped[["lower"]]), pc(x$capped[["upper"]])))
  for (side in names(x$held)[!is.na(x$held)])
    cat(sprintf("  %s bound held by the cap from lambda = %s (propensity ratio %s): flat from there\n",
                side, format(x$held[[side]], digits = digits),
                format(exp(x$held[[side]]), digits = digits)))
}

#' Pointwise sharp bounds
#'
#' Internal workhorse, exposed through `as.data.frame.sharp_scan()`.
#' @noRd
sharp_bounds <- function(eta_obs, r, lambda) {
  ell <- stats::qlogis(eta_obs)
  nlo <- stats::plogis(ell - lambda)
  nup <- stats::plogis(ell + lambda)
  data.frame(
    eta_obs = eta_obs, r = r,
    lower = pmax(nlo, r * eta_obs),
    upper = pmin(nup, r * eta_obs + 1 - r),
    naive_lower = nlo, naive_upper = nup
  )
}

# The smallest lambda at which f reaches threshold. The plug-in bounds are
# monotone in lambda, so a root on [0, lam_max] is the first one. Estimated
# bounds need not be: with `grid`, the first grid interval where the curve
# reaches the threshold is found first and the root is sought only there.
# `values`, the curve already computed on `grid`, saves evaluating f there.
first_crossing <- function(f, threshold, lam_max, side, grid = NULL, values = NULL) {
  g <- if (side == "upper") function(l) f(l) - threshold
       else function(l) threshold - f(l)
  if (g(0) >= 0) return(0)
  if (is.null(grid)) {
    if (g(lam_max) < 0) return(NA_real_)
    return(stats::uniroot(g, c(0, lam_max), tol = 1e-12)$root)
  }
  if (is.null(values)) values <- vapply(grid, f, 0)
  o <- order(grid); grid <- grid[o]; values <- values[o]
  gv <- if (side == "upper") values - threshold else threshold - values
  i <- which(gv >= 0)[1]
  if (is.na(i)) return(NA_real_)
  lo <- if (i == 1L) 0 else grid[i - 1L]        # g(0) < 0 is known
  if (grid[i] <= lo) return(grid[i])
  stats::uniroot(g, c(lo, grid[i]), tol = 1e-12)$root
}

check_inputs <- function(eta_obs, r) {
  if (!is.numeric(eta_obs) || !is.numeric(r))
    stop("`eta_obs` and `r` must be numeric vectors.", call. = FALSE)
  if (length(eta_obs) != length(r) || length(r) == 0L)
    stop("`eta_obs` and `r` must have the same, non-zero length.",
         call. = FALSE)
  if (anyNA(eta_obs) || anyNA(r))
    stop("`eta_obs` and `r` must not contain missing values.", call. = FALSE)
  if (any(eta_obs <= 0 | eta_obs >= 1))
    stop("`eta_obs` must lie strictly inside (0, 1).", call. = FALSE)
  if (any(r <= 0 | r > 1))
    stop("`r` must lie in (0, 1]. r = 0 means no overlap: nothing is ",
         "observed there and eta_obs is undefined.", call. = FALSE)
  invisible(TRUE)
}

check_weights <- function(weights, n) {
  if (is.null(weights)) return(NULL)
  if (!is.numeric(weights) || length(weights) != n || anyNA(weights) ||
      any(weights < 0) || sum(weights) <= 0)
    stop("`weights` must be non-negative, of the same length as `r`, ",
         "and not all zero.", call. = FALSE)
  weights
}

wmean <- function(x, w) if (is.null(w)) mean(x) else sum(w * x) / sum(w)

#' @rdname sharp_scan
#' @param x A `"sharp_scan"` object.
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#' @return `as.data.frame()`: the portfolio curve, or, if `lambda` is a single
#'   value, the pointwise sharp set and outer band for every row.
#' @export
as.data.frame.sharp_scan <- function(x, row.names = NULL, optional = FALSE,
                                      ..., lambda = NULL) {
  if (is.null(lambda)) return(x$curve)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0)
    stop("`lambda` must be a single non-negative number.", call. = FALSE)
  sharp_bounds(x$eta_obs, x$r, lambda)
}

#' @export
print.sharp_scan <- function(x, digits = 3, ...) {
  cv <- x$curve
  last <- nrow(cv)
  is_mean <- inherits(x, "sharp_scan_mean")
  cat(if (is_mean) "<sharp_scan_mean> sharp identified set for a mean"
      else "<sharp_scan> sharp identified set for a rate",
      "under a propensity ratio <= exp(lambda)\n")
  outer <- "outer"
  cat("  rows:", x$n, "  aggregation:", x$aggregate, "\n")
  cat("  lambda = 0 (model on observed rows):", format(x$point, digits = digits),
      "\n")
  cat(sprintf("  lambda = %s (propensity ratio %s): sharp [%s, %s], %s [%s, %s]\n",
              format(cv$lambda[last], digits = digits),
              format(exp(cv$lambda[last]), digits = digits),
              format(cv$lower[last], digits = digits),
              format(cv$upper[last], digits = digits), outer,
              format(cv$naive_lower[last], digits = digits),
              format(cv$naive_upper[last], digits = digits)))
  if (!is.null(x$threshold)) {
    fl <- function(v) if (is.na(v)) "not within the scanned range" else
      sprintf("lambda = %s (propensity ratio %s)", format(v, digits = digits),
              format(exp(v), digits = digits))
    cat("  threshold", format(x$threshold, digits = digits), "reached by\n")
    cat("    sharp set:  ", fl(x$flip[["sharp"]]), "\n")
    cat(sprintf("    %s band: ", outer), fl(x$flip[["naive"]]), "\n")
  }
  cap_lines(x, digits)
  invisible(x)
}

#' Extract the two scan inputs from fitted models
#'
#' A convenience, not an architectural commitment: calls `predict()` with
#' `type = "response"` on each model at `newdata`. For a classification
#' forest that returns class probabilities (e.g. `ranger` with
#' `probability = TRUE`) the column for class `"1"` (or the last column) is
#' used. Anything else, compute the two vectors yourself.
#'
#' @param eta_model Model of \eqn{Y} fitted on the observed rows.
#' @param r_model Model of \eqn{R} fitted on all rows.
#' @param newdata Rows at which the set is wanted.
#' @return A list with numeric vectors `eta_obs` and `r`.
#' @export
#' @examples
#' d <- simulate_shadow(3000, seed = 2)
#' f <- ~ x1 + x2 + z_shadow + z_policy + z_score
#' m_eta <- glm(update(f, y ~ .), binomial, data = d[d$observed == 1, ])
#' m_r <- glm(update(f, observed ~ .), binomial, data = d)
#' inp <- scan_inputs(m_eta, m_r, d[d$observed == 0, ])
#' sharp_scan(inp$eta_obs, inp$r)
scan_inputs <- function(eta_model, r_model, newdata) {
  list(eta_obs = response_prob(eta_model, newdata),
       r = response_prob(r_model, newdata))
}

response_prob <- function(m, newdata) {
  p <- stats::predict(m, newdata = newdata, type = "response")
  if (is.list(p) && !is.null(p$predictions)) p <- p$predictions
  if (is.matrix(p) || is.data.frame(p)) {
    p <- as.matrix(p)
    p <- if ("1" %in% colnames(p)) p[, "1"] else p[, ncol(p)]
  }
  as.numeric(p)
}
