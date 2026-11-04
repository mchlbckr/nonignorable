#' Simulate a nonignorably selected outcome with known truth
#'
#' Draws covariates, an outcome \eqn{Y} and an observation indicator \eqn{R}
#' whose propensity depends on \eqn{Y}, with the ratio of response
#' propensities between outcome values bounded by \eqn{e^{|\mathrm{tilt}|}}.
#' Because the mechanism is known, the target, the observed-data law and the
#' conditional observation rate \eqn{r(x) = P(R = 1 \mid X = x)} are returned
#' exactly. The simulator exists so that examples and tests depend on no
#' outside data.
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
#' Linear index of the outcome: \eqn{\mu = 0.8 x_1 - 0.5 x_2 + 0.6
#' z_{\mathrm{shadow}} + 1.0 z_{\mathrm{score}}}; of selection: \eqn{g = 0.8 -
#' 0.8 x_1 + 0.4 x_2 + 1.0 z_{\mathrm{policy}} - 1.0 z_{\mathrm{score}}}. All
#' truth columns condition on all five covariates.
#'
#' **Binary** (`outcome = "binary"`): \eqn{\eta = \mathrm{expit}(\mu - 1.5)},
#' \eqn{\pi(0,x) = \mathrm{expit}(g)\min(1, e^{-b})} and \eqn{\pi(1,x) =
#' \pi(0,x)e^{b}} with \eqn{b} = `tilt`. Truth: `eta`, `eta_obs`, `obs_rate`.
#'
#' **Continuous** (`outcome = "continuous"`): \eqn{Y = \mu + \varepsilon},
#' \eqn{\varepsilon \sim N(0,1)}, and \eqn{\pi(y,x) = c(x)\exp\{b\,
#' \Phi(\varepsilon)\}} with \eqn{c(x) = \mathrm{expit}(g)\min(1, e^{-b})}.
#' Since \eqn{\Phi(\varepsilon)} is uniform, the ratio of propensities is at
#' most \eqn{e^{|b|}}, \eqn{r(x) = c(x)(e^{b} - 1)/b}, and the quantile
#' function of \eqn{Y} given \eqn{R = 1, X = x} is \eqn{\mu(x) +
#' \Phi^{-1}[\log\{1 + p(e^{b} - 1)\}/b]}, all in closed form. Truth: `mu`
#' \eqn{= E(Y \mid x)}, `mean_obs` \eqn{= E(Y \mid R = 1, x)}, `obs_rate`;
#' [obs_quantiles()] returns the quantile grid that [shadow_scan_mean()]
#' takes.
#'
#' @param n Number of rows.
#' @param tilt The log of the largest ratio \eqn{b} of response propensities.
#'   Negative values mean that a high outcome (a default, \eqn{Y = 1}, or a
#'   large loss) lowers the chance of being observed, as when such cases are
#'   more likely to have been rejected. `0` is ignorable selection.
#' @param outcome `"binary"` (default) or `"continuous"`.
#' @param seed Seed for the random number generator. The caller's RNG state is
#'   restored on exit.
#'
#' @return A data frame with columns `x1`, `x2`, `z_shadow`, `z_policy`,
#'   `z_score`, `y` (the outcome, `NA` where unobserved), `observed` (0/1)
#'   and the exact truth columns described above. Attributes `"tilt"` and
#'   `"outcome"`.
#' @export
#' @examples
#' d <- simulate_shadow(2000, tilt = -1, seed = 1)
#' head(d)
#' dc <- simulate_shadow(2000, tilt = -1, outcome = "continuous", seed = 1)
#' head(dc)
#' dim(obs_quantiles(dc, K = 50))
simulate_shadow <- function(n = 5000, tilt = -1,
                            outcome = c("binary", "continuous"), seed = 1L) {
  outcome <- match.arg(outcome)
  stopifnot(length(n) == 1L, n >= 1, length(tilt) == 1L, is.finite(tilt))
  c_max <- min(1, exp(-tilt))
  with_seed(seed, {
    # Keep the order of the draws: standalone replications of this simulator
    # reproduce it draw for draw from the same seed.
    x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    z_shadow <- stats::rnorm(n); z_policy <- stats::rnorm(n)
    z_score <- stats::rnorm(n)
    mu <- 0.8 * x1 - 0.5 * x2 + 0.6 * z_shadow + 1.0 * z_score
    g <- 0.8 - 0.8 * x1 + 0.4 * x2 + 1.0 * z_policy - 1.0 * z_score
    if (outcome == "binary") {
      eta <- stats::plogis(mu - 1.5)
      pi0 <- stats::plogis(g) * c_max
      pi1 <- pi0 * exp(tilt)
      y_full <- stats::rbinom(n, 1L, eta)
      observed <- stats::rbinom(n, 1L, ifelse(y_full == 1L, pi1, pi0))
    } else {
      eps <- stats::rnorm(n)
      cx <- stats::plogis(g) * c_max
      y_full <- mu + eps
      observed <- stats::rbinom(n, 1L, cx * exp(tilt * stats::pnorm(eps)))
    }
  })
  d <- data.frame(
    x1 = x1, x2 = x2, z_shadow = z_shadow, z_policy = z_policy,
    z_score = z_score,
    y = ifelse(observed == 1L, y_full, NA),
    observed = observed
  )
  if (outcome == "binary") {
    d$y <- as.integer(d$y)
    obs_rate <- pi1 * eta + pi0 * (1 - eta)
    d$eta <- eta
    d$eta_obs <- pi1 * eta / obs_rate
    d$obs_rate <- obs_rate
  } else {
    d$mu <- mu
    d$mean_obs <- mu + tilted_eps_mean(tilt)
    d$obs_rate <- cx * uniform_mgf(tilt)
  }
  attr(d, "tilt") <- tilt
  attr(d, "outcome") <- outcome
  d
}

#' @rdname simulate_shadow
#' @param sim A data frame from `simulate_shadow(outcome = "continuous")`.
#' @param K Number of quantile levels, \eqn{(k - 1/2)/K}.
#' @return `obs_quantiles()`: the exact \eqn{n \times K} matrix of
#'   conditional quantiles of \eqn{Y} given \eqn{R = 1, X = x_i}.
#' @export
obs_quantiles <- function(sim, K = 200L) {
  if (!identical(attr(sim, "outcome"), "continuous"))
    stop("`sim` must come from simulate_shadow(outcome = \"continuous\").",
         call. = FALSE)
  qe <- tilted_eps_quantile(((seq_len(K)) - 0.5) / K, attr(sim, "tilt"))
  outer(sim$mu, qe, "+")
}

# E exp(b U), U uniform on [0, 1]
uniform_mgf <- function(b) if (b == 0) 1 else expm1(b) / b

# quantile of eps under density proportional to phi(e) exp(b Phi(e))
tilted_eps_quantile <- function(p, b) {
  if (b == 0) return(stats::qnorm(p))
  stats::qnorm(log1p(p * expm1(b)) / b)
}

tilted_eps_mean <- function(b) {
  if (b == 0) return(0)
  num <- stats::integrate(function(e) e * stats::dnorm(e) *
                            exp(b * stats::pnorm(e)), -Inf, Inf,
                          rel.tol = 1e-12)$value
  num / uniform_mgf(b)
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
