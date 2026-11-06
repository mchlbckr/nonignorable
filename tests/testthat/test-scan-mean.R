bounds1 <- function(q, r, lam, cap = TRUE) {
  q <- matrix(q, nrow = 1)
  unlist(nonignorable:::mean_bounds(nonignorable:::mean_prep(q), r, lam, cap))
}
qnorm_grid <- function(K) qnorm(((1:K) - 0.5) / K)

test_that("1. binary quantile rows reproduce the closed form of sharp_scan", {
  K <- 1000
  worst <- 0
  for (eta in c(0.05, 0.2, 0.5, 0.8)) for (r in c(0.1, 0.4, 0.7, 0.95))
    for (lam in c(0.1, 0.5, 1, log(3), 2, 5)) {
      n1 <- round(K * eta)
      q <- c(rep(0, K - n1), rep(1, n1))
      got <- bounds1(q, r, lam)
      cf <- nonignorable:::sharp_bounds(eta, r, lam)
      worst <- max(worst, abs(got - c(cf$lower, cf$upper)))
    }
  expect_lt(worst, 1e-12)
  # and the outer band without the cap is the naive band
  got <- bounds1(c(rep(0, 800), rep(1, 200)), 0.6, 1, cap = FALSE)
  expect_equal(unname(got), plogis(qlogis(0.2) + c(-1, 1)), tolerance = 1e-12)
})

test_that("reference values for N(0,1), K = 20000", {
  q <- qnorm_grid(20000)
  ref <- rbind(c(0.9, 0.5, 0.1485, 0.1993), c(0.7, 1.0, 0.3817, 0.3975),
               c(0.5, 1.0, 0.3975, 0.3975), c(0.7, 2.0, 0.5811, 0.7865))
  for (i in seq_len(nrow(ref))) {
    expect_lt(abs(bounds1(q, ref[i, 1], ref[i, 2])[["upper"]] - ref[i, 3]),
              6e-5)
    expect_lt(abs(bounds1(q, ref[i, 1], ref[i, 2], cap = FALSE)[["upper"]] -
                    ref[i, 4]), 6e-5)
  }
})

test_that("2. lambda = 0 gives the observed mean; monotone in lambda", {
  set.seed(41)
  q <- t(apply(matrix(rexp(50 * 100), 50), 1, sort))
  r <- runif(50, 0.05, 1)
  s <- sharp_scan_mean(q, r, lambda = seq(0, 3, by = 0.1))
  expect_equal(s$curve$lower[1], mean(q))
  expect_equal(s$curve$upper[1], mean(q))
  expect_equal(s$point, mean(q))
  expect_true(all(diff(s$curve$upper) >= -1e-12))
  expect_true(all(diff(s$curve$lower) <= 1e-12))
  # pointwise too
  rows <- sapply(c(0.5, 1, 2), function(l) as.data.frame(s, lambda = l)$upper)
  expect_true(all(apply(rows, 1, diff) >= -1e-12))
})

test_that("3. never wider than the outer band; equal when the cap is slack", {
  set.seed(42)
  q <- t(apply(matrix(rnorm(40 * 200) * runif(40, .5, 3), 40), 1, sort))
  r <- runif(40, 0.02, 1)
  for (lam in c(0.3, 1, 2.5)) {
    b <- nonignorable:::mean_bounds(nonignorable:::mean_prep(q), r, lam, TRUE)
    o <- nonignorable:::mean_bounds(nonignorable:::mean_prep(q), r, lam, FALSE)
    expect_true(all(b$upper <= o$upper + 1e-12 & b$lower >= o$lower - 1e-12))
    # s_circ: the maximiser over [0, 1]; cap slack iff sbar >= s_circ
    K <- ncol(q); kap <- expm1(lam)
    for (i in seq_len(nrow(q))) {
      ys <- sort(q[i, ], decreasing = TRUE)
      f <- (mean(ys) + kap * c(0, cumsum(ys)) / K) / (1 + kap * (0:K) / K)
      s_circ <- ((0:K) / K)[which.max(f)]
      sbar <- min(1, (1 - r[i]) / (r[i] * kap))
      if (sbar >= s_circ) expect_equal(b$upper[i], o$upper[i], tolerance = 1e-12)
      else expect_lt(b$upper[i], o$upper[i])
    }
  }
})

test_that("4. symmetric quantiles give a symmetric set", {
  q <- qnorm_grid(2000)
  for (r in c(0.2, 0.7, 0.95)) for (lam in c(0.4, 1.5)) {
    b <- bounds1(q, r, lam)
    expect_equal(b[["lower"]], -b[["upper"]], tolerance = 1e-12)
  }
})

test_that("5. bounded Y and large lambda approach the Manski bounds", {
  q <- qbeta(((1:500) - 0.5) / 500, 2, 3)
  m <- mean(q)
  for (r in c(0.3, 0.6, 0.9)) {
    b <- bounds1(q, r, 40)
    manski <- c(r * m + (1 - r) * min(q), r * m + (1 - r) * max(q))
    expect_equal(unname(b), manski, tolerance = 1e-10)
    # and never beyond them at moderate lambda
    b2 <- bounds1(q, r, 2)
    expect_true(b2[["lower"]] >= manski[1] - 1e-12 &&
                b2[["upper"]] <= manski[2] + 1e-12)
  }
})

test_that("6. random admissible weights never beat the bounds", {
  set.seed(20261104)
  for (trial in 1:6) {
    K <- 7
    y <- sort(rnorm(K) * 2)
    r <- runif(1, 0.2, 0.9); lam <- runif(1, 0.2, 2.5); kap <- expm1(lam)
    b <- bounds1(y, r, lam)
    n_try <- 30000
    u <- matrix(runif(n_try * K), n_try)^runif(n_try, 0.2, 5)
    a <- 1 / (r * (1 + kap * rowMeans(u)))
    ok <- a >= 1
    w <- a[ok] * (1 + kap * u[ok, , drop = FALSE])  # w = 1/pi >= 1, ratio <= e^lam
    expect_lt(max(abs(rowMeans(w) - 1 / r)), 1e-10)
    expect_true(all(apply(w, 1, max) / apply(w, 1, min) <= exp(lam) + 1e-12))
    mu <- r * as.vector(w %*% y) / K
    best <- max(mu); worst <- min(mu)
    expect_lte(best, b[["upper"]] + 1e-12)
    expect_gte(worst, b[["lower"]] - 1e-12)
    expect_gt(best, b[["lower"]])        # the search does get somewhere
  }
})

test_that("extremal weights attain the upper bound", {
  y <- sort(c(-1.2, 0.3, 0.5, 2.0, 3.1)); K <- 5
  r <- 0.6; lam <- 1.1; kap <- expm1(lam)
  ub <- bounds1(y, r, lam)[["upper"]]
  ys <- sort(y, decreasing = TRUE)
  sbar <- min(1, (1 - r) / (r * kap))
  ss <- sort(unique(c((0:K) / K, sbar)))
  ss <- ss[ss <= sbar]
  U <- function(s) { u <- pmin(pmax(s * K - (0:(K - 1)), 0), 1); sum(u * ys) / K }
  f <- sapply(ss, function(s) (mean(y) + kap * U(s)) / (1 + kap * s))
  s_star <- ss[which.max(f)]
  u <- pmin(pmax(s_star * K - (0:(K - 1)), 0), 1)
  a <- 1 / (r * (1 + kap * s_star))
  w <- a * (1 + kap * u)
  expect_true(min(w) >= 1 - 1e-12 && max(w) / min(w) <= exp(lam) + 1e-12)
  expect_equal(mean(w), 1 / r, tolerance = 1e-12)
  expect_equal(r * mean(w * ys), ub, tolerance = 1e-12)
})

test_that("7. truth on the simulator is inside from lambda = |tilt| on", {
  for (tilt in c(-1.2, -0.4, 0, 0.8)) {
    d <- simulate_shadow(1500, tilt = tilt, outcome = "continuous", seed = 43)
    q <- obs_quantiles(d, K = 1000)
    expect_equal(rowMeans(q), d$mean_obs, tolerance = 1e-3)
    rej <- d$observed == 0
    s <- sharp_scan_mean(q[rej, ], d$obs_rate[rej],
                          lambda = abs(tilt) + c(0, 0.5))
    truth <- mean(d$mu[rej])
    expect_true(all(truth >= s$curve$lower & truth <= s$curve$upper))
    rows <- as.data.frame(s, lambda = abs(tilt))
    expect_true(all(d$mu[rej] >= rows$lower - 1e-9 &
                    d$mu[rej] <= rows$upper + 1e-9))
  }
})

test_that("simulator: continuous truth is internally consistent", {
  d <- simulate_shadow(200000, tilt = -1, outcome = "continuous", seed = 44)
  # observed rate and observed mean agree with the closed forms
  expect_equal(mean(d$observed), mean(d$obs_rate), tolerance = 0.01)
  expect_equal(mean(d$y[d$observed == 1] - d$mu[d$observed == 1]),
               mean(d$mean_obs - d$mu), tolerance = 0.01)
  q <- obs_quantiles(d[1:3, ], K = 4)
  expect_equal(dim(q), c(3L, 4L))
  expect_error(obs_quantiles(simulate_shadow(10)), "continuous")
})

test_that("flip, print, as.data.frame and validation", {
  d <- simulate_shadow(2000, tilt = -1, outcome = "continuous", seed = 45)
  q <- obs_quantiles(d, K = 200)
  rej <- d$observed == 0
  s <- sharp_scan_mean(q[rej, ], d$obs_rate[rej], threshold = 0.5)
  expect_s3_class(s, c("sharp_scan_mean", "sharp_scan"))
  expect_false(is.na(s$flip[["sharp"]]))
  at <- sharp_scan_mean(q[rej, ], d$obs_rate[rej], lambda = s$flip[["sharp"]])
  expect_equal(at$curve$upper, 0.5, tolerance = 1e-8)
  expect_gte(s$flip[["sharp"]], s$flip[["naive"]] - 1e-9)
  expect_output(print(s), "outer band")
  expect_equal(nrow(as.data.frame(s, lambda = 1)), sum(rej))
  w <- runif(sum(rej))
  sw <- sharp_scan_mean(q[rej, ], d$obs_rate[rej], lambda = 0, weights = w)
  expect_equal(sw$point, sum(w * rowMeans(q[rej, ])) / sum(w))
  expect_error(sharp_scan_mean(q, d$obs_rate[1:2]), "one value per row")
  expect_error(sharp_scan_mean(q[1:2, ], c(0, 1)), "overlap")
  expect_error(sharp_scan_mean(q[1:2, ], c(1, 1), lambda = -1), "non-negative")
  qq <- q[1:2, ]; qq[1, 1] <- NA
  expect_error(sharp_scan_mean(qq, c(1, 1)), "finite")
  # r = 1 pins the set to the observed mean
  b <- bounds1(qnorm_grid(100), 1, 3)
  expect_equal(unname(b), rep(mean(qnorm_grid(100)), 2), tolerance = 1e-12)
  # a vector is a single row
  expect_equal(sharp_scan_mean(qnorm_grid(10), 0.5, lambda = 1)$n, 1L)
  pdf(NULL); on.exit(dev.off())
  expect_identical(plot(s), s)
})

test_that("unobserved = TRUE bounds the mean of the rows' own missing outcomes", {
  d <- simulate_shadow(20000, tilt = -1, outcome = "continuous", seed = 82)
  mis <- d[d$observed == 0, ]
  q <- obs_quantiles(mis, K = 200); r <- mis$obs_rate
  s <- sharp_scan_mean(q, r, lambda = c(0, 1), unobserved = TRUE)
  b <- sharp_scan_mean(q, r, lambda = c(0, 1))
  m <- rowMeans(q)
  pr <- as.data.frame(s, lambda = 1); pb <- as.data.frame(b, lambda = 1)
  expect_equal(pr$lower, (pb$lower - r * m) / (1 - r))
  expect_equal(s$point, b$point)
  truth <- mean((mis$mu - r * mis$mean_obs) / (1 - r))
  expect_lte(s$curve$lower[2], truth); expect_gte(s$curve$upper[2], truth)
  # binary input reproduces the binary scan here too
  e <- c(0.1, 0.4, 0.7); rr <- c(0.3, 0.6, 0.8); K <- 1000
  qb <- t(sapply(e, function(p) as.numeric(((seq_len(K) - 0.5) / K) > 1 - p)))
  expect_equal(sharp_scan_mean(qb, rr, lambda = 1, unobserved = TRUE)$curve$upper,
               sharp_scan(e, rr, lambda = 1, unobserved = TRUE)$curve$upper, tolerance = 2e-3)
})
