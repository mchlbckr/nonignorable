contrast_fixture <- function() {
  d <- simulate_shadow(4000, tilt = -1, seed = 91)
  f <- ~ x1 + x2 + z_shadow + z_policy + z_score
  inp <- crossfit_inputs(d, "observed",
    fit_eta = function(train, new) predict(glm(update(f, y ~ .), binomial, train),
                                           new, type = "response"),
    fit_r   = function(train, new) predict(glm(update(f, observed ~ .), binomial, train),
                                           new, type = "response"),
    seed = 1)
  list(d = d, inp = inp, A = d$x1 > 0, B = d$x1 <= 0)
}
fx <- contrast_fixture()
lam <- seq(0, 2, length.out = 41)
k_ab <- function(a = fx$A, b = fx$B, ...)
  shadow_scan_contrast(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y,
                       target_a = a, target_b = b, lambda = lam, ...)

test_that("contrast endpoints are differences of the arm scans", {
  k <- k_ab()
  ca <- shadow_scan_ci(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y, target = fx$A, lambda = lam)
  cb <- shadow_scan_ci(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y, target = fx$B, lambda = lam)
  expect_equal(k$curve$upper, ca$curve$upper - cb$curve$lower)
  expect_equal(k$curve$lower, ca$curve$lower - cb$curve$upper)
  expect_equal(k$curve$naive_upper, ca$curve$naive_upper - cb$curve$naive_lower)
  expect_equal(k$curve$naive_lower, ca$curve$naive_lower - cb$curve$naive_upper)
})

test_that("contrast se is the root sum of squares", {
  k <- k_ab()
  ca <- shadow_scan_ci(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y, target = fx$A, lambda = lam)
  cb <- shadow_scan_ci(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y, target = fx$B, lambda = lam)
  expect_equal(k$curve$se_upper, sqrt(ca$curve$se_upper^2 + cb$curve$se_lower^2))
  expect_equal(k$curve$se_lower, sqrt(ca$curve$se_lower^2 + cb$curve$se_upper^2))
  z <- qnorm(0.95)
  expect_equal(k$curve$ci_upper, pmin(k$curve$upper + z * k$curve$se_upper, 1))
  expect_equal(k$curve$ci_lower, pmax(k$curve$lower - z * k$curve$se_lower, -1))
})

test_that("contrast flips on the right side", {
  k1 <- k_ab(); k2 <- k_ab(fx$B, fx$A)
  expect_equal(k1$point, -k2$point)
  expect_equal(k1$flip, k2$flip, tolerance = 1e-8)
  expect_false(is.na(k1$flip[["sharp"]]))
  expect_lte(k1$flip[["ci"]], k1$flip[["sharp"]])
})

test_that("contrast validates targets", {
  expect_error(k_ab(fx$A, fx$A | fx$B), "disjoint")
  expect_error(k_ab(fx$A, rep(FALSE, length(fx$A))), "empty")
  expect_error(k_ab(as.numeric(fx$A), fx$B), "logical")
  expect_error(shadow_scan_contrast(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y,
                                    fx$A, fx$B, lambda = c(0.5, 1)), "lambda = 0")
})

test_that("contrast prints and plots", {
  k <- k_ab()
  expect_s3_class(k, c("shadow_scan_contrast", "shadow_scan_ci", "shadow_scan"))
  expect_output(print(k), "difference")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(k))
})

test_that("contrast needs a grid of at least two lambdas", {
  expect_error(shadow_scan_contrast(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y,
                                    fx$A, fx$B, lambda = 0), "at least two")
})

test_that("contrast labels the warnings of each arm", {
  eta <- fx$inp$eta_obs; r <- fx$inp$r
  bad <- which(fx$d$observed == 1 & fx$d$y %in% 1 & fx$A)[1:3]
  r[bad] <- 1e-5; eta[bad] <- 0.01
  w <- character()
  withCallingHandlers(
    shadow_scan_contrast(eta, r, fx$d$observed, fx$d$y, fx$A, fx$B,
                         lambda = c(0, 1), labels = c("arm one", "arm two")),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("^arm one: ", w)))
  expect_false(any(grepl("^arm two: ", w)))
})

test_that("the contrast print names the outer band", {
  d <- simulate_shadow(1500, tilt = -1, seed = 71)
  arm <- d$x1 > 0
  k <- shadow_scan_contrast(rep(0.3, nrow(d)), rep(0.6, nrow(d)), d$observed, d$y,
                            target_a = arm, target_b = !arm)
  out <- paste(capture.output(print(k)), collapse = "\n")
  expect_match(out, "outer band \\(plug-in\\)")
  expect_no_match(out, "naive")
})

test_that("the contrast can compare the unobserved rows of two segments", {
  k <- k_ab(among = "unobserved")
  arm <- function(tg) shadow_scan_ci(fx$inp$eta_obs, fx$inp$r, fx$d$observed, fx$d$y,
                                     target = "unobserved", segment = tg, lambda = lam)$curve
  ca <- arm(fx$A); cb <- arm(fx$B)
  expect_equal(k$curve$upper, cummax(ca$upper - cb$lower))
  expect_equal(k$curve$se_lower, sqrt(ca$se_lower^2 + cb$se_upper^2))
  expect_identical(k$n_a, sum(fx$A & fx$d$observed == 0))
  expect_identical(k$among, "unobserved")
  expect_output(print(k), "unobserved rows")
  expect_error(k_ab(among = "observed"), "should be one of")
})
