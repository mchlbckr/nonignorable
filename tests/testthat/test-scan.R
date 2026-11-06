brute_set <- function(eta, r, lam, m = 200001) {
  t <- seq(-lam, lam, length.out = m)
  et <- plogis(qlogis(eta) - t)
  ct <- r / (exp(t) * et + 1 - et)
  ok <- ct > 0 & ct <= 1 & ct * exp(t) <= 1
  range(et[ok])
}

test_that("closed form agrees with the definition of the set", {
  set.seed(11)
  for (i in 1:300) {
    e <- runif(1, 0.01, 0.99); r <- runif(1, 0.01, 0.999); l <- runif(1, 0, 3)
    cf <- unlist(nonignorable:::sharp_bounds(e, r, l)[c("lower", "upper")])
    expect_equal(unname(cf), brute_set(e, r, l), tolerance = 1e-4)
  }
})

test_that("endpoints are attained by admissible propensities", {
  set.seed(12)
  e <- runif(500, 0.01, 0.99); r <- runif(500, 0.01, 1); l <- 1.3
  b <- nonignorable:::sharp_bounds(e, r, l)
  for (end in c("lower", "upper")) {
    v <- b[[end]]
    t <- qlogis(e) - qlogis(v)
    pi0 <- r / (exp(t) * v + 1 - v)
    pi1 <- pi0 * exp(t)
    tol <- 1e-10
    expect_true(all(abs(t) <= l + tol))
    expect_true(all(pi0 > 0 & pi0 <= 1 + tol))
    expect_true(all(pi1 > 0 & pi1 <= 1 + tol))
    # the propensities reproduce the observed law
    expect_equal(pi1 * v + pi0 * (1 - v), r, tolerance = 1e-10)
    expect_equal(pi1 * v / r, e, tolerance = 1e-10)
  }
})

test_that("set collapses to eta_obs at lambda = 0", {
  e <- c(0.1, 0.5, 0.9); r <- c(0.2, 0.6, 1)
  b <- nonignorable:::sharp_bounds(e, r, 0)
  expect_equal(b$lower, e)
  expect_equal(b$upper, e)
  s <- sharp_scan(e, r, lambda = c(0, 1))
  expect_identical(s$curve$lower[1], mean(e))
  expect_identical(s$curve$upper[1], mean(e))
  expect_identical(s$point, mean(e))
})

test_that("never wider than naive; strictly narrower exactly where Manski cuts", {
  set.seed(13)
  e <- runif(2000, 0.01, 0.99); r <- runif(2000, 0.01, 1); l <- 0.8
  b <- nonignorable:::sharp_bounds(e, r, l)
  expect_true(all(b$lower >= b$naive_lower))
  expect_true(all(b$upper <= b$naive_upper))
  cut_lo <- r * e > plogis(qlogis(e) - l)
  cut_up <- r * e + 1 - r < plogis(qlogis(e) + l)
  expect_identical(b$lower > b$naive_lower, cut_lo)
  expect_identical(b$upper < b$naive_upper, cut_up)
  # no cut at all where r <= exp(-lambda)
  small <- r <= exp(-l)
  expect_false(any(cut_lo[small] | cut_up[small]))
})

test_that("scan contains the known truth on the simulator", {
  for (tilt in c(-1.5, -0.5, 0, 0.7)) {
    d <- simulate_shadow(4000, tilt = tilt, seed = 3)
    lam <- abs(tilt)
    b <- nonignorable:::sharp_bounds(d$eta_obs, d$obs_rate, lam)
    tol <- 1e-12
    expect_true(all(d$eta >= b$lower - tol & d$eta <= b$upper + tol))
    rej <- d[d$observed == 0, ]
    for (agg in c("common", "union")) {
      s <- sharp_scan(rej$eta_obs, rej$obs_rate, lambda = lam, aggregate = agg)
      truth <- mean(rej$eta)
      expect_gte(truth, s$curve$lower - tol)
      expect_lte(truth, s$curve$upper + tol)
    }
  }
})

test_that("common tilt is inside union and inside naive", {
  d <- simulate_shadow(3000, tilt = -1, seed = 4)
  lam <- seq(0, 2, by = 0.25)
  sc <- sharp_scan(d$eta_obs, d$obs_rate, lam, aggregate = "common")$curve
  su <- sharp_scan(d$eta_obs, d$obs_rate, lam, aggregate = "union")$curve
  tol <- 1e-12
  expect_true(all(sc$lower >= su$lower - tol & sc$upper <= su$upper + tol))
  expect_true(all(su$lower >= su$naive_lower - tol &
                  su$upper <= su$naive_upper + tol))
  expect_equal(sc$naive_lower, su$naive_lower)
})

test_that("common tilt reports the binding row and respects it", {
  e <- c(0.2, 0.2, 0.2); r <- c(0.1, 0.1, 0.99)
  s <- sharp_scan(aggregate = "common", e, r, lambda = c(0.5, 2))
  # row 3 has r near one and pins the common tilt
  expect_equal(s$binding$upper_row, c(3L, 3L))
  expect_equal(s$binding$lower_row, c(3L, 3L))
  # at the binding tilt row 3 sits exactly at its Manski bound
  t_up <- qlogis(0.2) - qlogis(0.99 * 0.2 + 0.01)
  expect_equal(s$curve$upper, rep(plogis(qlogis(0.2) - t_up), 2))
  # non-binding: the naive bound binds
  s2 <- sharp_scan(aggregate = "common", c(0.2, 0.3), c(0.1, 0.2), lambda = 0.5)
  expect_true(is.na(s2$binding$upper_row))
  expect_equal(s2$curve$upper, s2$curve$naive_upper)
})

test_that("union is the default aggregation", {
  expect_identical(sharp_scan(0.2, 0.5)$aggregate, "union")
  expect_null(sharp_scan(0.2, 0.5)$binding)
})

test_that("weights enter the portfolio mean", {
  e <- c(0.1, 0.4); r <- c(0.5, 0.5)
  s <- sharp_scan(e, r, lambda = 0, weights = c(3, 1))
  expect_equal(s$point, (3 * 0.1 + 0.4) / 4)
})

test_that("flip point is exact and sharp flips no earlier than naive", {
  d <- simulate_shadow(5000, tilt = -1, seed = 1)
  rej <- d[d$observed == 0, ]
  s <- sharp_scan(aggregate = "common", rej$eta_obs, rej$obs_rate, lambda = seq(0, 2, by = 0.1),
                   threshold = 0.35)
  expect_false(is.na(s$flip[["sharp"]]))
  up <- sharp_scan(aggregate = "common", rej$eta_obs, rej$obs_rate, lambda = s$flip[["sharp"]])
  expect_equal(up$curve$upper, 0.35, tolerance = 1e-9)
  expect_gte(s$flip[["sharp"]], s$flip[["naive"]] - 1e-9)
  # a threshold the naive band reaches but the common-tilt set never does
  s2 <- sharp_scan(aggregate = "common", rej$eta_obs, rej$obs_rate, lambda = seq(0, 2, by = 0.1),
                    threshold = 0.4)
  expect_true(is.na(s2$flip[["sharp"]]))
  expect_false(is.na(s2$flip[["naive"]]))
  # unreachable threshold
  s3 <- sharp_scan(rej$eta_obs, rej$obs_rate, lambda = c(0, 0.1),
                    threshold = 0.99)
  expect_true(is.na(s3$flip[["sharp"]]))
  # threshold below the point uses the lower end
  s4 <- sharp_scan(rej$eta_obs, rej$obs_rate, lambda = c(0, 2),
                    aggregate = "union", threshold = s$point / 2)
  lo <- sharp_scan(rej$eta_obs, rej$obs_rate, lambda = s4$flip[["sharp"]],
                    aggregate = "union")
  expect_equal(lo$curve$lower, s$point / 2, tolerance = 1e-9)
})

test_that("r = 1 pins the set to a point", {
  b <- nonignorable:::sharp_bounds(0.3, 1, 2)
  expect_equal(c(b$lower, b$upper), c(0.3, 0.3))
})

test_that("inputs are validated", {
  expect_error(sharp_scan(c(0.1, 0.2), 0.5), "same")
  expect_error(sharp_scan(0, 0.5), "strictly inside")
  expect_error(sharp_scan(0.2, 0), "overlap")
  expect_error(sharp_scan(0.2, 0.5, lambda = -1), "non-negative")
  expect_error(sharp_scan(0.2, NA_real_), "missing")
  expect_error(sharp_scan(0.2, 0.5, threshold = 2), "threshold")
})

test_that("as.data.frame returns curve or pointwise sets", {
  s <- sharp_scan(c(0.1, 0.5), c(0.3, 0.9), lambda = c(0, 1))
  expect_identical(as.data.frame(s), s$curve)
  rows <- as.data.frame(s, lambda = 1)
  expect_equal(nrow(rows), 2L)
  expect_true(all(c("lower", "upper", "naive_lower", "naive_upper") %in%
                  names(rows)))
})

test_that("scan_inputs extracts predictions", {
  d <- simulate_shadow(1500, seed = 5)
  m_eta <- glm(y ~ x1 + x2, binomial, data = d[d$observed == 1, ])
  m_r <- glm(observed ~ x1 + x2, binomial, data = d)
  inp <- scan_inputs(m_eta, m_r, d[1:10, ])
  expect_equal(inp$eta_obs, unname(predict(m_eta, d[1:10, ], type = "response")))
  expect_length(inp$r, 10L)
})

test_that("first_crossing finds the first crossing of a curve that is not monotone", {
  f <- function(l) (l - 0.2) * (l - 0.4) * (l - 1.8)      # crosses 0 at 0.2, 0.4, 1.8
  expect_equal(nonignorable:::first_crossing(f, 0, 2, "upper",
                                           grid = seq(0, 2, length.out = 81)),
               0.2, tolerance = 1e-8)
})
