mfx <- local({
  d <- simulate_shadow(4000, tilt = -1, outcome = "continuous", seed = 93)
  list(d = d, q = obs_quantiles(d, K = 50), A = d$x1 > 0, B = d$x1 <= 0)
})
mlam <- seq(0, 2, length.out = 21)
mk <- function(a = mfx$A, b = mfx$B, ...)
  shadow_scan_mean_contrast(mfx$q, mfx$d$obs_rate, mfx$d$observed, mfx$d$y,
                            target_a = a, target_b = b, lambda = mlam, ...)

test_that("mean contrast endpoints are differences of the arm scans", {
  k <- mk()
  arm <- function(tg) shadow_scan_mean_ci(mfx$q, mfx$d$obs_rate, mfx$d$observed, mfx$d$y,
                                          target = tg, lambda = mlam)$curve
  ca <- arm(mfx$A); cb <- arm(mfx$B)
  expect_equal(k$curve$upper, ca$upper - cb$lower)
  expect_equal(k$curve$lower, ca$lower - cb$upper)
  expect_equal(k$curve$se_upper, sqrt(ca$se_upper^2 + cb$se_lower^2))
  expect_equal(k$curve$naive_lower, ca$naive_lower - cb$naive_upper)
  z <- qnorm(0.95)
  expect_equal(k$curve$ci_upper, cummax(k$curve$upper + z * k$curve$se_upper))
  expect_equal(k$point, k$curve$upper[1])
  expect_s3_class(k, c("shadow_scan_contrast", "shadow_scan_ci", "shadow_scan"))
  expect_true(k$of_means)
})

test_that("mean contrast flips, prints and plots as a difference of means", {
  k <- mk(threshold = 0)
  expect_true(all(c("sharp", "naive", "ci") %in% names(k$flip)))
  expect_true(is.na(k$flip[["ci"]]) || k$flip[["ci"]] <= k$flip[["sharp"]] || is.na(k$flip[["sharp"]]))
  out <- paste(capture.output(print(k)), collapse = "\n")
  expect_match(out, "difference of means A - B")
  expect_match(out, "outer band \\(plug-in\\)")
  pdf(NULL); on.exit(dev.off())
  expect_identical(plot(k), k)
})

test_that("mean contrast checks its targets and threshold", {
  expect_error(mk(mfx$A, mfx$A), "disjoint")
  expect_error(mk(mfx$A, rep(FALSE, nrow(mfx$d))), "empty")
  expect_error(mk(threshold = NA_real_), "threshold")
  expect_error(shadow_scan_mean_contrast(mfx$q, mfx$d$obs_rate, mfx$d$observed, mfx$d$y,
                                         mfx$A, mfx$B, lambda = c(0.5, 1)), "lambda = 0")
})
