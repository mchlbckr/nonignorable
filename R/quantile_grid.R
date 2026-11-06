#' Quantile grid for the mean scan
#'
#' `quantile_grid()` produces the \eqn{n \times K} matrix of conditional
#' quantiles of \eqn{Y} given \eqn{R = 1, X = x_i} at the levels
#' \eqn{(k - 1/2)/K} that [sharp_scan_mean()] takes as `q`. It is a
#' convenience like [scan_inputs()], not a commitment to a learner, and
#' accepts three kinds of `fit`:
#'
#' * a quantile regression from \pkg{quantreg} (`rq` or `rqs`), fitted on the
#'   observed rows. It is refitted at the `K` levels with `update()`, so the
#'   data it was fitted on must still be reachable, and predicted at
#'   `newdata`;
#' * a function `function(tau, newdata)` returning, for a vector of levels
#'   `tau`, a matrix with one row per row of `newdata` and one column per
#'   level -- for any other conditional quantile model;
#' * a numeric matrix of draws from the conditional distribution, one row per
#'   row of interest and any number of draws per row -- for a distributional
#'   or generative model. `newdata` is then not used.
#'
#' Quantile curves estimated level by level can cross. Each row is sorted
#' (rearranged), which yields a proper quantile function and never moves it
#' further from the truth (Chernozhukov, Fernández-Val and Galichon, 2010).
#' With many levels nearly every row crosses somewhere by sampling noise
#' alone, so the share of rows that needed it (attribute `"crossing"`) says
#' little. Attribute `"shift"` says how much was moved: the median over rows of
#' the largest change by sorting, relative to the row's range. A few percent is
#' noise; a large shift says the quantile model is poor.
#'
#' @param fit An `rq`/`rqs` object, a `function(tau, newdata)`, or a numeric
#'   matrix of draws.
#' @param newdata Rows at which the quantiles are wanted (for `rq` and
#'   function fits).
#' @param K Number of levels. Default 200.
#' @return A numeric matrix, `nrow(newdata)` (or `nrow(fit)` for draws) by
#'   `K`, rows sorted, with attributes `"levels"`, `"crossing"` and `"shift"`.
#' @references
#' Chernozhukov, V., Fernández-Val, I. and Galichon, A. (2010). Quantile and
#' probability curves without crossing. *Econometrica*, 78(3), 1093--1125.
#' \doi{10.3982/ECTA7880}
#' @seealso [sharp_scan_mean()]
#' @export
#' @examples
#' dc <- simulate_shadow(3000, tilt = -1, outcome = "continuous", seed = 1)
#' obs <- dc[dc$observed == 1, ]
#' rej <- dc[dc$observed == 0, ]
#' # any conditional quantile model, here a location model fitted by lm()
#' m <- lm(y ~ x1 + x2 + z_shadow + z_policy + z_score, data = obs)
#' res <- residuals(m)
#' q <- quantile_grid(function(tau, newdata)
#'   outer(predict(m, newdata), quantile(res, tau, names = FALSE), "+"),
#'   rej, K = 100)
#' sm <- sharp_scan_mean(q, rej$obs_rate, threshold = 0.5)
#' sm
#' if (requireNamespace("quantreg", quietly = TRUE)) {
#'   fq <- quantreg::rq(y ~ x1 + x2 + z_shadow + z_policy + z_score,
#'                      tau = 0.5, data = obs)
#'   q2 <- quantile_grid(fq, rej, K = 50)
#'   attr(q2, "crossing")
#' }
quantile_grid <- function(fit, newdata = NULL, K = 200L) {
  if (!is.numeric(K) || length(K) != 1L || K < 2)
    stop("`K` must be a single number >= 2.", call. = FALSE)
  K <- as.integer(K)
  tau <- (seq_len(K) - 0.5) / K
  if (is.matrix(fit) && is.numeric(fit)) {
    if (anyNA(fit)) stop("The draws must not contain missing values.", call. = FALSE)
    q <- t(apply(fit, 1L, stats::quantile, probs = tau, names = FALSE, type = 8))
    if (K == 1L) q <- t(q)
  } else {
    if (is.null(newdata))
      stop("`newdata` is required unless `fit` is a matrix of draws.", call. = FALSE)
    if (is.function(fit)) {
      q <- fit(tau, newdata)
    } else if (inherits(fit, c("rq", "rqs"))) {
      if (!requireNamespace("quantreg", quietly = TRUE))
        stop("Package 'quantreg' is needed for an rq fit.", call. = FALSE)
      # re-evaluate where the model was fitted, not here
      call <- fit$call
      call$tau <- tau
      env <- environment(stats::formula(fit))
      if (is.null(env)) env <- parent.frame()
      refit <- eval(call, env)
      q <- stats::predict(refit, newdata = newdata)
    } else {
      stop("`fit` must be an rq/rqs object, a function(tau, newdata) or a ",
           "numeric matrix of draws.", call. = FALSE)
    }
    q <- as.matrix(q)
    if (nrow(q) != nrow(newdata) || ncol(q) != K)
      stop("The quantile model must return one row per row of `newdata` and ",
           "one column per level.", call. = FALSE)
    if (anyNA(q) || any(!is.finite(q)))
      stop("The quantile model returned missing or infinite values.", call. = FALSE)
  }
  q <- unname(q)
  crossed <- apply(q, 1L, is.unsorted)
  shift <- numeric(nrow(q))
  if (any(crossed)) {
    raw <- q[crossed, , drop = FALSE]
    srt <- t(apply(raw, 1L, sort))
    span <- srt[, K] - srt[, 1L]
    shift[crossed] <- ifelse(span > 0, apply(abs(srt - raw), 1L, max) / span, 0)
    q[crossed, ] <- srt
  }
  attr(q, "levels") <- tau
  attr(q, "crossing") <- mean(crossed)
  attr(q, "shift") <- stats::median(shift)
  q
}
