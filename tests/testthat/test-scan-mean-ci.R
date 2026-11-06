# Reference score for the upper bound, after ~/Projects/mnar-sharp/spike/10_orthogonal_score.R
# (functions prep, at, side_parts, score), generalised from a location model with
# one residual grid to one quantile grid per row. One deliberate change: the
# spike finds the kink index as floor(s* K) in floating point, which at a few
# kinks lands one grid value too high; here it is the exact integer index of
# the maximum (MEAN-CI-DESIGN.md, decision 9).
ref_score <- function(q, r, Y, R, lam) {
  kap <- expm1(lam)
  vapply(seq_len(nrow(q)), function(i) {
    qs <- sort(q[i, ], decreasing = TRUE); K <- length(qs)
    V <- c(0, cumsum(qs)) / K; s <- (0:K) / K; m <- mean(qs)
    g <- (m + kap * V) / (1 + kap * s)
    jstar <- which.max(g) - 1L
    sbar <- min(1, (1 - r[i]) / (r[i] * kap))
    cap <- sbar < jstar / K
    if (cap) { so <- sbar; j <- min(floor(sbar * K), K - 1L) }
    else { so <- jstar / K; j <- min(jstar, K - 1L) }
    qa <- qs[j + 1]; Vc <- V[j + 1] + (so - s[j + 1]) * qa
    theta <- (m + kap * Vc) / (1 + kap * so); D <- 1 + kap * so
    theta + R[i] * (Y[i] + kap * max(Y[i] - qa, 0) + kap * qa * so - theta * D) / (r[i] * D) -
      cap * (qa - theta) * (R[i] - r[i]) / (r[i]^2 * D)
  }, numeric(1))
}

mean_fixture <- local({
  d <- simulate_shadow(3000, tilt = -1, outcome = "continuous", seed = 7)
  list(d = d, q = obs_quantiles(d, K = 50),
       r = pmin(pmax(d$obs_rate, 0.01), 0.99),
       Y = ifelse(d$observed == 1, d$y, 0), R = d$observed)
})

test_that("the score equals the reference from the paper's spike", {
  f <- mean_fixture
  pr <- shadowscan:::mean_prep(f$q)
  for (lam in c(0, 0.5, 2)) {
    up <- shadowscan:::mean_score(pr$hi, pr$m, f$r, f$Y, f$R, lam)
    lo <- shadowscan:::mean_score(pr$lo, -pr$m, f$r, -f$Y, f$R, lam)
    expect_equal(up$phi, ref_score(f$q, f$r, f$Y, f$R, lam), tolerance = 1e-10)
    expect_equal(-lo$phi, -ref_score(-f$q, f$r, -f$Y, f$R, lam), tolerance = 1e-10)
  }
  cap2 <- shadowscan:::mean_score(pr$hi, pr$m, f$r, f$Y, f$R, 2)$cap
  expect_true(any(cap2))
  expect_true(any(!cap2))
})

test_that("theta is the plug-in bound of shadow_scan_mean", {
  f <- mean_fixture
  pr <- shadowscan:::mean_prep(f$q)
  for (lam in c(0.5, 2)) {
    up <- shadowscan:::mean_score(pr$hi, pr$m, f$r, f$Y, f$R, lam)
    expect_equal(up$theta, shadowscan:::mean_upper(pr$hi, pr$m, f$r, lam, TRUE),
                 tolerance = 1e-12)
  }
})

test_that("constant grids, observed Y beyond the grid and K = 1 stay finite", {
  q <- rbind(c(1, 2, 3, 4), c(5, 5, 5, 5), c(-1, 0, 1, 8))
  r <- c(0.5, 0.3, 0.9); R <- c(1L, 1L, 1L); Y <- c(100, 5, -50)
  pr <- shadowscan:::mean_prep(q)
  for (lam in c(0.5, 2)) {
    expect_no_warning(up <- shadowscan:::mean_score(pr$hi, pr$m, r, Y, R, lam))
    expect_true(all(is.finite(up$phi)))
    expect_equal(up$phi, ref_score(q, r, Y, R, lam), tolerance = 1e-12)
  }
  q1 <- matrix(c(1, 2, 3), ncol = 1)
  pr1 <- shadowscan:::mean_prep(q1)
  up1 <- shadowscan:::mean_score(pr1$hi, pr1$m, c(0.4, 0.5, 0.95), c(2, 0, 9), c(1L, 0L, 1L), 1)
  expect_true(all(is.finite(up1$phi)))
  expect_equal(up1$theta, c(1, 2, 3))
})

ci_fixture <- local({
  d <- simulate_shadow(20000, tilt = -1, outcome = "continuous", seed = 12)
  list(d = d, q = obs_quantiles(d, K = 100))
})

test_that("with the true inputs the debiased bounds hit the truth", {
  big <- simulate_shadow(100000, tilt = -1, outcome = "continuous", seed = 11)
  truth <- shadow_scan_mean(obs_quantiles(big, K = 100), big$obs_rate, lambda = c(0, 1))
  f <- ci_fixture
  k <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y, lambda = c(0, 1))
  # the truth is itself computed on another sample, so compare within 3 SE
  expect_true(all(abs(k$curve$upper - truth$curve$upper) < 3 * k$curve$se_upper))
  expect_true(all(abs(k$curve$lower - truth$curve$lower) < 3 * k$curve$se_lower))
  expect_true(all(k$curve$ci_lower <= k$curve$lower & k$curve$ci_upper >= k$curve$upper))
  expect_true(all(k$curve$se_upper > 0 & k$curve$se_lower > 0))
  expect_s3_class(k, c("shadow_scan_mean_ci", "shadow_scan_ci", "shadow_scan_mean", "shadow_scan"))
  expect_equal(k$point, k$curve$upper[1])
})

test_that("targets: all, a segment, and unobserved is refused", {
  f <- ci_fixture
  seg <- f$d$x1 > 0
  k <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y, target = seg,
                           lambda = c(0, 1))
  expect_equal(k$n_target, sum(seg))
  kall <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y, lambda = c(0, 1))
  expect_equal(kall$n_target, nrow(f$d))
  expect_false(isTRUE(all.equal(k$curve$upper, kall$curve$upper)))
  expect_error(shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y,
                                   target = "unobserved"), "not covered")
})

test_that("r-hat is clipped into r_bounds and the rows are counted", {
  f <- ci_fixture
  r <- f$d$obs_rate; r[1:5] <- 0.999; r[6:8] <- 0.001
  k <- shadow_scan_mean_ci(f$q, r, f$d$observed, f$d$y, lambda = c(0, 1))
  expect_equal(k$clipped, c(low = sum(r < 0.01), high = sum(r > 0.99)))
  expect_gte(k$clipped[["low"]], 3); expect_gte(k$clipped[["high"]], 5)
  expect_output(print(k), "clipped.*below 0.01.*above 0.99")
  for (bad in list(replace(r, 1, 0), replace(r, 1, 1.2), replace(r, 1, NA)))
    expect_error(shadow_scan_mean_ci(f$q, bad, f$d$observed, f$d$y), "`r`")
})

test_that("trim removes target rows and says how many", {
  f <- ci_fixture
  k <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y, lambda = c(0, 1),
                           trim = 0.2)
  expect_equal(k$trimmed, sum(f$d$obs_rate < 0.2))
  expect_equal(k$n_target, sum(f$d$obs_rate >= 0.2))
  expect_output(print(k), "trimmed")
})

test_that("inputs are validated", {
  f <- ci_fixture
  y <- f$d$y
  expect_true(anyNA(y[f$d$observed == 0]))   # NA on unobserved rows is fine
  y[which(f$d$observed == 1)[1]] <- NA
  expect_error(shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, y), "outcome")
  expect_error(shadow_scan_mean_ci(f$q[-1, ], f$d$obs_rate, f$d$observed, f$d$y), "one row")
  expect_error(shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y,
                                   lambda = c(-1, 0)), "lambda")
  expect_error(shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y,
                                   threshold = Inf), "threshold")
  expect_error(shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y,
                                   r_bounds = c(0.5, 0.1)), "r_bounds")
  expect_no_error(shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y, lambda = 0))
})

test_that("flip points, print and plot", {
  f <- ci_fixture
  k0 <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y, lambda = c(0, 2))
  thr <- k0$point + 0.5 * (k0$curve$upper[2] - k0$point)
  k <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y,
                           lambda = seq(0, 2, length.out = 21), threshold = thr)
  fl <- k$flip
  expect_false(is.na(fl[["sharp"]])); expect_false(is.na(fl[["ci"]]))
  expect_true(fl[["ci"]] <= fl[["sharp"]])
  expect_output(print(k), "shadow_scan_mean_ci")
  pdf(NULL); on.exit(dev.off())
  expect_identical(plot(k), k)
})

test_that("crossfit_quantiles fits out of fold on observed rows", {
  d <- simulate_shadow(600, tilt = -1, outcome = "continuous", seed = 21)
  d$id <- seq_len(nrow(d))
  seen <- list(); taus <- list()
  fq <- function(train, newdata, tau) {
    seen[[length(seen) + 1]] <<- list(obs = all(train$observed == 1), ids = train$id,
                                      new = newdata$id)
    taus[[length(taus) + 1]] <<- tau
    outer(rep(mean(train$y), nrow(newdata)), tau, "+")
  }
  fr <- function(train, newdata) rep(mean(train$observed), nrow(newdata))
  cf <- crossfit_quantiles(d, "observed", fq, fr, K = 20, folds = 4, seed = 3)
  expect_equal(dim(cf$q), c(600, 20))
  expect_length(cf$r, 600)
  expect_true(all(vapply(seen, `[[`, logical(1), "obs")))
  for (s in seen) expect_length(intersect(s$ids, s$new), 0)       # out of fold
  expect_setequal(unlist(lapply(seen, `[[`, "new")), d$id)        # each row once
  expect_equal(taus[[1]], (seq_len(20) - 0.5) / 20)
  expect_equal(cf$fold, crossfit_quantiles(d, "observed", fq, fr, K = 20, folds = 4,
                                           seed = 3)$fold)
})

test_that("crossfit_quantiles checks the shape of fit_q and the response", {
  d <- simulate_shadow(300, tilt = -1, outcome = "continuous", seed = 22)
  fr <- function(train, newdata) rep(0.5, nrow(newdata))
  expect_error(crossfit_quantiles(d, "observed", function(tr, nd, tau) rep(0, nrow(nd)),
                                  fr, K = 5), "one row per")
  expect_error(crossfit_quantiles(d, "observed",
                                  function(tr, nd, tau) matrix(0, nrow(nd), 3), fr, K = 5),
               "one row per")
  d$observed[2] <- NA
  expect_error(crossfit_quantiles(d, "observed", function(tr, nd, tau)
    matrix(0, nrow(nd), length(tau)), fr, K = 5), "0/1 or logical")
})

test_that("a constant grid gives the same score whatever its value", {
  # g is flat; floating-point noise must not pick an arbitrary kink.
  # With s* = 0: theta = c, q = c, D = 1, phi = c + (1 + kappa) (Y - c) / r.
  for (cc in c(0.1, 3.3, -7.9)) {
    q <- matrix(cc, 2, 200)
    pr <- shadowscan:::mean_prep(q)
    sc <- shadowscan:::mean_score(pr$hi, pr$m, c(0.3, 0.3), c(cc + 1, 0), c(1L, 0L), 2)
    expect_false(any(sc$cap))
    expect_equal(sc$phi, c(cc + exp(2) / 0.3, cc), tolerance = 1e-10)
  }
})

test_that("per-row bounds come out of as.data.frame, as for shadow_scan_mean", {
  d <- simulate_shadow(1500, tilt = -1, outcome = "continuous", seed = 1)
  q <- obs_quantiles(d, K = 50)
  k <- shadow_scan_mean_ci(q, d$obs_rate, d$observed, d$y, lambda = c(0, 1))
  rows <- as.data.frame(k, lambda = 1)
  expect_equal(nrow(rows), nrow(d))
  ref <- as.data.frame(shadow_scan_mean(q, pmin(pmax(d$obs_rate, 0.01), 0.99)), lambda = 1)
  expect_equal(rows$upper, ref$upper)
  expect_equal(rows$lower, ref$lower)
})

test_that("a numeric outcome is required, with a message that says so", {
  f <- mean_fixture
  expect_error(shadow_scan_mean_ci(f$q, f$r, f$R, f$d$y > 0), "must be numeric")
  expect_error(shadow_scan_mean_ci(f$q, f$r, f$R, factor(f$d$y > 0)), "must be numeric")
})

test_that("NA level or trim are refused by the package, not by base R", {
  f <- mean_fixture
  expect_error(shadow_scan_mean_ci(f$q, f$r, f$R, f$d$y, level = NA_real_), "`level`")
  expect_error(shadow_scan_mean_ci(f$q, f$r, f$R, f$d$y, trim = NA_real_), "`trim`")
})

test_that("clipped and the tail count target rows, below the lower bound", {
  f <- mean_fixture
  r <- f$d$obs_rate
  tg <- f$d$x1 > 0
  k <- shadow_scan_mean_ci(f$q, r, f$R, f$d$y, target = tg, lambda = c(0, 1),
                           r_bounds = c(0.05, 0.99))
  expect_equal(k$clipped, c(low = sum(tg & r < 0.05), high = sum(tg & r > 0.99)))
  expect_equal(k$tail[["target"]], sum(tg & r < 0.05))
  expect_output(print(k), "r-hat < 0.05")
})

test_that("the print labels the outer band as plug-in", {
  f <- mean_fixture
  k <- shadow_scan_mean_ci(f$q, f$r, f$R, f$d$y, lambda = c(0, 1), threshold = 0)
  expect_output(print(k), "outer band \\(plug-in\\)")
})

test_that("the range warning names the quantiles when r-hat is harmless", {
  f <- mean_fixture
  # quantiles squeezed to a tiny range around the mean, r-hat well inside (0, 1)
  qq <- mean(f$d$y[f$R == 1]) + f$q * 1e-3
  w <- tryCatch(shadow_scan_mean_ci(qq, rep(0.5, nrow(qq)), f$R, f$d$y,
                                    lambda = c(0, 1)),
                warning = conditionMessage)
  expect_match(w, "do not cover the observed outcomes")
  expect_no_match(w, "blown up")
})

test_that("crossfit_quantiles refuses a fractional K and non-numeric fits", {
  d <- simulate_shadow(300, tilt = -1, outcome = "continuous", seed = 2)
  fr <- function(train, newdata) rep(0.5, nrow(newdata))
  fq <- function(train, newdata, tau) matrix(0, nrow(newdata), length(tau))
  expect_error(crossfit_quantiles(d, "observed", fq, fr, K = 1.5), "`K`")
  fq_chr <- function(train, newdata, tau) matrix("a", nrow(newdata), length(tau))
  expect_error(crossfit_quantiles(d, "observed", fq_chr, fr, K = 5, seed = 1),
               "numeric")
})

test_that("clipping only at the top says so, without alarm", {
  f <- ci_fixture
  r <- pmax(f$d$obs_rate, 0.02); r[1:30] <- 0.995
  k <- shadow_scan_mean_ci(f$q, r, f$d$observed, f$d$y, lambda = c(0, 1))
  expect_identical(k$clipped[["low"]], 0L)
  out <- paste(capture.output(print(k)), collapse = "\n")
  expect_match(out, "clipped: 30 target rows with r-hat above 0.99")
  expect_no_match(out, "below 0.01, ")
})

test_that("debiased bounds and confidence sets are monotone in lambda", {
  f <- ci_fixture
  k <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y,
                           lambda = seq(0, 2, by = 0.1))
  cv <- k$curve
  expect_true(all(diff(cv$lower) <= 0) && all(diff(cv$upper) >= 0))
  expect_true(all(diff(cv$ci_lower) <= 0) && all(diff(cv$ci_upper) >= 0))
})

test_that("the mean print reports where the cap binds", {
  f <- ci_fixture
  k <- shadow_scan_mean_ci(f$q, f$d$obs_rate, f$d$observed, f$d$y, lambda = c(0, 2))
  expect_named(k$capped, c("lower", "upper"))
  expect_true(all(k$capped >= 0 & k$capped <= 1))
  expect_output(print(k), "cap pi <= 1 at lambda = 2 binds in")
})
