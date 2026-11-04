#' Simulate a nonignorably selected binary outcome with known truth
#'
#' Draws covariates, a binary outcome \eqn{Y} and an observation indicator
#' \eqn{R} whose propensity depends on \eqn{Y} through a constant logit tilt
#' \eqn{b = \log\{\pi(1,x)/\pi(0,x)\}}. Because the mechanism is known, the
#' target \eqn{\eta(x) = P(Y = 1 \mid X = x)}, the observed-data probability
#' \eqn{\eta_{\mathrm{obs}}(x) = P(Y = 1 \mid X = x, R = 1)} and the conditional
#' observation rate \eqn{r(x) = P(R = 1 \mid X = x)} are returned exactly. The
#' simulator exists so that examples and tests depend on no outside data.
#'
#' The covariates play fixed roles, so that [shadow_screen()] has something to
#' find:
#' * `x1`, `x2`: ordinary covariates entering both outcome and selection.
#' * `z_shadow`: enters the outcome only, a valid shadow variable.
#' * `z_policy`: enters the selection only, high exclusion violation and no
#'   relevance.
#' * `z_score`: enters both strongly, like a bureau score that is also used in
#'   the acceptance decision.
#'
#' Outcome: \eqn{\mathrm{logit}\,\eta = -1.5 + 0.8 x_1 - 0.5 x_2 + 0.6 z_{\mathrm{shadow}} + 1.0 z_{\mathrm{score}}}.
#' Selection: \eqn{\pi(0,x) = \mathrm{expit}(g(x)) \cdot \min(1, e^{-b})} and
#' \eqn{\pi(1,x) = \pi(0,x)e^{b}}, with
#' \eqn{g = 0.8 - 0.8 x_1 + 0.4 x_2 + 1.0 z_{\mathrm{policy}} - 1.0 z_{\mathrm{score}}}.
#' All truth columns condition on all five covariates.
#'
#' @param n Number of rows.
#' @param tilt The log ratio \eqn{b} of the observation propensities of
#'   \eqn{Y = 1} and \eqn{Y = 0}. Negative values mean that a bad outcome
#'   (\eqn{Y = 1}) lowers the chance of being observed, as when defaults are
#'   more likely to have been rejected. `0` is ignorable selection.
#' @param seed Seed for the random number generator. The caller's RNG state is
#'   restored on exit.
#'
#' @return A data frame with columns `x1`, `x2`, `z_shadow`, `z_policy`,
#'   `z_score`, `y` (the outcome, `NA` where unobserved), `observed` (0/1),
#'   and the exact truth `eta`, `eta_obs` and `obs_rate`. The tilt is stored
#'   in attribute `"tilt"`.
#' @export
#' @examples
#' d <- simulate_shadow(2000, tilt = -1, seed = 1)
#' head(d)
#' mean(d$observed)
simulate_shadow <- function(n = 5000, tilt = -1, seed = 1L) {
  stopifnot(length(n) == 1L, n >= 1, length(tilt) == 1L, is.finite(tilt))
  with_seed(seed, {
    x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    z_shadow <- stats::rnorm(n); z_policy <- stats::rnorm(n)
    z_score <- stats::rnorm(n)
    eta <- stats::plogis(-1.5 + 0.8 * x1 - 0.5 * x2 + 0.6 * z_shadow + 1.0 * z_score)
    g <- 0.8 - 0.8 * x1 + 0.4 * x2 + 1.0 * z_policy - 1.0 * z_score
    pi0 <- stats::plogis(g) * min(1, exp(-tilt))
    pi1 <- pi0 * exp(tilt)
    y_full <- stats::rbinom(n, 1L, eta)
    observed <- stats::rbinom(n, 1L, ifelse(y_full == 1L, pi1, pi0))
  })
  obs_rate <- pi1 * eta + pi0 * (1 - eta)
  d <- data.frame(
    x1 = x1, x2 = x2, z_shadow = z_shadow, z_policy = z_policy,
    z_score = z_score,
    y = ifelse(observed == 1L, y_full, NA_integer_),
    observed = observed,
    eta = eta, eta_obs = pi1 * eta / obs_rate, obs_rate = obs_rate
  )
  attr(d, "tilt") <- tilt
  d
}

with_seed <- function(seed, code) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit(
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv),
    add = TRUE
  )
  set.seed(seed)
  eval.parent(substitute(code))
}
