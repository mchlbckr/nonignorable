test_that("draws reproduce known quantiles and the scan input shape", {
  set.seed(71)
  draws <- matrix(rnorm(50 * 4000), 50)
  q <- quantile_grid(draws, K = 20)
  expect_equal(dim(q), c(50, 20))
  expect_equal(attr(q, "levels"), ((1:20) - 0.5) / 20)
  expect_lt(max(abs(colMeans(q) - qnorm(attr(q, "levels")))), 0.06)
  expect_equal(attr(q, "crossing"), 0)
})

test_that("a function fit is used as given and crossings are rearranged", {
  nd <- data.frame(a = 1:5)
  f <- function(tau, newdata) outer(newdata$a, tau, function(a, t) a + qnorm(t))
  q <- quantile_grid(f, nd, K = 10)
  expect_equal(q[3, ], 3 + qnorm(((1:10) - 0.5) / 10), ignore_attr = TRUE)
  g <- function(tau, newdata) {
    m <- outer(newdata$a, tau, function(a, t) a + qnorm(t))
    m[2, ] <- rev(m[2, ])                         # a crossing row
    m
  }
  q2 <- quantile_grid(g, nd, K = 10)
  expect_equal(attr(q2, "crossing"), 0.2)
  expect_equal(attr(q2, "shift"), 0)            # median over rows: 4 of 5 untouched
  expect_gt(attr(quantile_grid(function(tau, nd) g(tau, nd)[c(2, 2, 2), ], nd[1:3, , drop = FALSE], K = 10), "shift"), 0.5)
  expect_false(is.unsorted(q2[2, ]))
  expect_error(quantile_grid(function(tau, newdata) 1, nd, K = 10), "one row per row")
  expect_error(quantile_grid(f, NULL), "newdata")
  expect_error(quantile_grid("x", nd), "must be")
})

test_that("an rq fit is refitted at the grid and feeds the mean scan", {
  skip_if_not_installed("quantreg")
  dc <- simulate_shadow(2000, tilt = -1, outcome = "continuous", seed = 72)
  obs <- dc[dc$observed == 1, ]; rej <- dc[dc$observed == 0, ]
  fq <- quantreg::rq(y ~ x1 + x2 + z_shadow + z_policy + z_score, tau = 0.5, data = obs)
  q <- quantile_grid(fq, rej, K = 40)
  expect_equal(dim(q), c(nrow(rej), 40))
  # close to the true observed-data quantiles of the simulator
  qt <- obs_quantiles(dc, K = 40)[dc$observed == 0, ]
  expect_lt(abs(mean(q) - mean(qt)), 0.05)
  sm <- sharp_scan_mean(q, rej$obs_rate, lambda = c(0, 1))
  expect_s3_class(sm, "sharp_scan_mean")
})
