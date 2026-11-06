# true eta_t(x) = P(Y <= t | R = 1, x) of the continuous simulator: given
# R = 1, u = Phi(eps) has density proportional to exp(b u) on (0, 1)
true_eta_t <- function(d, t) {
  b <- attr(d, "tilt")
  sapply(t, function(tj) (exp(b * pnorm(tj - d$mu)) - 1) / (exp(b) - 1))
}

test_that("each t reproduces the rate scan for 1(Y <= t)", {
  d <- simulate_shadow(1500, tilt = -1, outcome = "continuous", seed = 2)
  tg <- c(-1, 0, 1)
  et <- true_eta_t(d, tg)
  for (tgt in c("all", "unobserved")) {
    k <- sharp_scan_cdf(et, d$obs_rate, d$observed, d$y, t = tg, lambda = 0.7,
                        target = tgt)
    for (j in seq_along(tg)) {
      z <- ifelse(d$observed == 1, as.integer(d$y <= tg[j]), NA)
      s <- sharp_scan_ci(et[, j], d$obs_rate, d$observed, z, target = tgt,
                         lambda = c(0, 0.7))
      expect_equal(k$curve$lower[j], s$curve$lower[2])
      expect_equal(k$curve$upper[j], s$curve$upper[2])
      expect_equal(k$curve$ci_upper[j], s$curve$ci_upper[2])
      expect_equal(k$curve$point[j], s$point)
    }
  }
})

test_that("the confidence sets cover the true distribution function", {
  d <- simulate_shadow(20000, tilt = -1, outcome = "continuous", seed = 3)
  tg <- c(-2, -1, 0, 1, 2)
  k <- sharp_scan_cdf(true_eta_t(d, tg), d$obs_rate, d$observed, d$y, t = tg,
                      lambda = 1)
  truth <- sapply(tg, function(tj) mean(pnorm(tj - d$mu)))
  expect_true(all(k$curve$ci_lower <= truth & truth <= k$curve$ci_upper))
  ku <- sharp_scan_cdf(true_eta_t(d, tg), d$obs_rate, d$observed, d$y, t = tg,
                       lambda = 1, target = "unobserved")
  mis <- d$observed == 0
  truth_u <- sapply(tg, function(tj) mean(pnorm(tj - d$mu[mis])))
  expect_true(all(ku$curve$ci_lower <= truth_u & truth_u <= ku$curve$ci_upper))
  # high outcomes are less often observed: the missing rows lie higher, so
  # MAR overstates their distribution function
  expect_true(all(ku$curve$point > truth_u))
})

test_that("quantile() inverts the bounds on the grid", {
  d <- simulate_shadow(3000, tilt = -1, outcome = "continuous", seed = 5)
  tg <- seq(-4, 3, by = 0.25)
  k <- sharp_scan_cdf(true_eta_t(d, tg), d$obs_rate, d$observed, d$y, t = tg,
                      lambda = 0.5)
  q <- quantile(k, c(0.1, 0.5, 0.9))
  expect_named(q, c("prob", "point", "lower", "upper"))
  expect_true(all(q$lower <= q$point & q$point <= q$upper))
  cv <- k$curve
  # widened outwards on the grid: F+ is still below 0.5 at the reported Q-,
  # F- has reached it at the reported Q+, and not one grid step earlier
  j <- match(q$lower[2], tg)
  expect_lt(sort(cv$upper)[j], 0.5)
  expect_gte(sort(cv$upper)[j + 1], 0.5)
  j <- match(q$upper[2], tg)
  expect_gte(sort(cv$lower)[j], 0.5)
  expect_lt(sort(cv$lower)[j - 1], 0.5)
  # the point is interpolated, so it need not lie on the grid
  expect_equal(approx(sort(cv$point), tg, 0.5)$y, q$point[2])
  # a level the grid does not reach
  expect_true(is.na(quantile(k, 0.999)$upper))
  expect_error(quantile(k, 1), "probs")
})

test_that("crossfit_cdf fits one classification per t", {
  d <- simulate_shadow(600, tilt = -1, outcome = "continuous", seed = 6)
  calls <- 0L
  fe <- function(train, newdata) {
    calls <<- calls + 1L
    expect_true(all(train$z %in% 0:1))
    predict(glm(z ~ x1 + x2, binomial, train), newdata, type = "response")
  }
  fr <- function(train, newdata)
    predict(glm(observed ~ x1 + x2, binomial, train), newdata, type = "response")
  tg <- c(-100, -1, 0, 1)
  cf <- crossfit_cdf(d, "observed", "y", tg, fe, fr, folds = 3, seed = 1)
  expect_equal(dim(cf$eta_t), c(600L, 4L))
  expect_equal(calls, 3L * 3L)          # the t below all outcomes is not fitted
  expect_true(all(cf$eta_t[, 1] == 1e-6))
  expect_true(all(diff(colMeans(cf$eta_t)) > 0))
  expect_error(crossfit_cdf(d, "observed", "y", c(1, 0), fe, fr), "increasing")
})

test_that("bad input is refused", {
  d <- simulate_shadow(300, tilt = -1, outcome = "continuous", seed = 7)
  et <- true_eta_t(d, c(0, 1))
  expect_error(sharp_scan_cdf(et, d$obs_rate, d$observed, d$y, t = 0), "one per column")
  expect_error(sharp_scan_cdf(et, d$obs_rate, d$observed, d$y, t = c(1, 0)), "increasing")
  expect_error(sharp_scan_cdf(et, d$obs_rate, d$observed, d$y, t = c(0, 1),
                              lambda = c(0, 1)), "single")
  y <- d$y; y[which(d$observed == 1)[1]] <- NA
  expect_error(sharp_scan_cdf(et, d$obs_rate, d$observed, y, t = c(0, 1)), "NA")
})

test_that("print and plot work, and the kink share is reported", {
  d <- simulate_shadow(800, tilt = -1, outcome = "continuous", seed = 8)
  tg <- seq(-2, 2, by = 1)
  k <- sharp_scan_cdf(true_eta_t(d, tg), d$obs_rate, d$observed, d$y, t = tg,
                      lambda = 1, target = "unobserved")
  expect_true(all(k$curve$near_lower >= 0 & k$curve$near_upper <= 1))
  expect_output(print(k), "pointwise")
  f <- tempfile(fileext = ".png"); grDevices::png(f)
  on.exit({grDevices::dev.off(); unlink(f)})
  expect_silent(plot(k))
})
