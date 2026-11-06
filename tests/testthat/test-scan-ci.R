ci_sim <- function(n = 2000, seed = 61) {
  d <- simulate_shadow(n, tilt = -1, seed = seed)
  f <- ~ x1 + x2 + z_shadow + z_policy + z_score
  inp <- crossfit_inputs(d, "observed",
    fit_eta = function(tr, new) predict(glm(update(f, y ~ .), binomial, tr), new, type = "response"),
    fit_r = function(tr, new) predict(glm(update(f, observed ~ .), binomial, tr), new, type = "response"),
    seed = 1)
  list(d = d, inp = inp)
}

test_that("crossfit predictions come from other folds", {
  d <- simulate_shadow(500, seed = 62)
  seen <- integer()
  inp <- crossfit_inputs(d, "observed",
    fit_eta = function(tr, new) { seen <<- c(seen, nrow(tr)); rep(0.3, nrow(new)) },
    fit_r = function(tr, new) rep(0.5, nrow(new)), folds = 4, seed = 2)
  expect_length(seen, 4)
  expect_true(all(seen < sum(d$observed)))     # never all observed rows
  expect_identical(sort(unique(inp$fold)), 1:4)
  expect_equal(inp$eta_obs, rep(0.3, 500))
  inp2 <- crossfit_inputs(d, "observed", function(tr, new) rep(0.3, nrow(new)),
                          function(tr, new) rep(0.5, nrow(new)), folds = 4, seed = 2)
  expect_identical(inp$fold, inp2$fold)
  expect_error(crossfit_inputs(d, "observed", function(tr, new) 0.3,
                               function(tr, new) rep(0.5, nrow(new))), "too few|missing")
})

test_that("the confidence set contains the estimate and widens with lambda", {
  x <- ci_sim()
  s <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, target = "unobserved")
  cv <- s$curve
  expect_true(all(cv$ci_lower <= cv$lower & cv$upper <= cv$ci_upper))
  expect_true(all(cv$se_lower > 0 & cv$se_upper > 0))
  expect_gt(diff(range(cv$ci_upper - cv$ci_lower)), 0)
  expect_gt(cv$ci_upper[nrow(cv)] - cv$ci_lower[nrow(cv)], cv$ci_upper[1] - cv$ci_lower[1])
  expect_s3_class(s, "sharp_scan")
})

test_that("with in-sample truth the estimate is the plug-in when no Manski piece binds", {
  # perfect nuisances and no residual: the correction vanishes on average
  x <- ci_sim(4000, 63)
  d <- x$d
  s <- sharp_scan_ci(d$eta_obs, d$obs_rate, d$observed, d$y, target = "all",
                      lambda = 0)
  expect_lt(abs(s$point - mean(d$eta_obs)), 4 * s$curve$se_upper[1])
})

test_that("targets: unobserved, all and a logical segment", {
  x <- ci_sim()
  a <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, target = "all",
                      lambda = c(0, 1))
  u <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, lambda = c(0, 1), target = "unobserved")
  g <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y,
                      target = x$d$x1 > 0, lambda = c(0, 1))
  expect_false(isTRUE(all.equal(a$curve$upper, u$curve$upper)))
  expect_identical(u$target, "unobserved")
  expect_identical(g$n_target, sum(x$d$x1 > 0))
  expect_error(sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y,
                              target = "nope"), "should be one of")
})

test_that("threshold gives three crossings, the confidence set first", {
  x <- ci_sim()
  s <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, threshold = 0.35, target = "unobserved")
  expect_named(s$flip, c("sharp", "naive", "ci"))
  expect_lte(s$flip[["ci"]], s$flip[["sharp"]])
  expect_output(print(s), "confidence set")
  pdf(NULL); on.exit(dev.off())
  expect_identical(plot(s), s)
})

test_that("inputs are validated", {
  x <- ci_sim(500, 64)
  e <- x$inp$eta_obs; r <- x$inp$r; R <- x$d$observed; Y <- x$d$y
  expect_error(sharp_scan_ci(e, r, R[-1], Y), "one per row")
  expect_error(sharp_scan_ci(e, r, R, Y[-1]), "one value per row")
  Y2 <- Y; Y2[which(R == 1)[1]] <- NA
  expect_error(sharp_scan_ci(e, r, R, Y2), "0/1 where")
  expect_error(sharp_scan_ci(e, r, R, Y, level = 1), "level")
  expect_error(sharp_scan_ci(e, r, R, Y, target = rep(TRUE, 3)), "logical vector")
})

test_that("both standard errors are positive and of similar size; trim shrinks the target", {
  x <- ci_sim()
  a <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, lambda = c(0.5, 2), target = "unobserved")
  b <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, lambda = c(0.5, 2),
                      se_type = "empirical", target = "unobserved")
  expect_equal(a$curve$upper, b$curve$upper)
  ratio <- a$curve$se_upper / b$curve$se_upper
  expect_true(all(ratio > 0.5 & ratio < 2))
  t <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, lambda = 1, trim = 0.05, target = "unobserved")
  expect_gt(t$trimmed, 0)
  expect_lt(t$n_target, a$n_target)
  expect_output(print(t), "trimmed")
  expect_error(sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, trim = 1, target = "unobserved"), "trim")
})

test_that("the default standard error is the larger of the two", {
  x <- ci_sim()
  args <- list(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, lambda = c(0, 1, 2))
  m <- do.call(sharp_scan_ci, args)
  cnd <- do.call(sharp_scan_ci, c(args, se_type = "conditional"))
  emp <- do.call(sharp_scan_ci, c(args, se_type = "empirical"))
  expect_equal(m$curve$se_upper, pmax(cnd$curve$se_upper, emp$curve$se_upper))
  expect_equal(m$curve$se_lower, pmax(cnd$curve$se_lower, emp$curve$se_lower))
  expect_identical(m$se_type, "max")
})

test_that("blown-up inverse weights and separation are flagged", {
  d <- simulate_shadow(2000, tilt = -1, seed = 81)
  eta <- pmin(pmax(d$eta_obs, 0.01), 0.99); r <- d$obs_rate
  ok <- sharp_scan_ci(eta, r, d$observed, d$y, lambda = c(0, 1), target = "unobserved")
  expect_named(ok$tail, c("target", "observed", "min_r_observed"))
  expect_output(print(ok), "r-hat < 0.01")
  bad <- which(d$observed == 1 & d$y == 1)[1:3]
  r2 <- r; r2[bad] <- 1e-5; eta2 <- eta; eta2[bad] <- 0.01
  expect_warning(sharp_scan_ci(eta2, r2, d$observed, d$y, lambda = c(0, 1), target = "unobserved"),
                 "separation")
  expect_warning(
    crossfit_inputs(d, "observed",
                    fit_eta = function(train, new) rep(0.3, nrow(new)),
                    fit_r = function(train, new) ifelse(new$observed == 1 & new$y %in% 1, 1e-9, 0.5),
                    seed = 1),
    "separates")
})

test_that("trimmed counts only target rows that left the target", {
  d <- simulate_shadow(2000, tilt = -1, seed = 82)
  eta <- pmin(pmax(d$eta_obs, 0.01), 0.99)
  s <- sharp_scan_ci(eta, d$obs_rate, d$observed, d$y, lambda = c(0, 1), trim = 0.2, target = "unobserved")
  expect_equal(s$trimmed, sum(d$observed == 0 & d$obs_rate < 0.2))
  k <- sharp_scan_contrast(eta, d$obs_rate, d$observed, d$y, d$x1 > 0, d$x1 <= 0,
                            lambda = c(0, 1), trim = 0.2)
  expect_equal(k$trimmed, sum(d$obs_rate < 0.2))
})

test_that("crossfit_inputs rejects a response column with NA", {
  d <- simulate_shadow(500, tilt = -1, seed = 83)
  d$observed[3] <- NA
  expect_error(crossfit_inputs(d, "observed",
                               function(tr, nd) rep(0.3, nrow(nd)),
                               function(tr, nd) rep(0.5, nrow(nd)), seed = 1),
               "0/1 or logical")
})

test_that("the print names the outer band and reports where the cap binds", {
  x <- ci_sim(1000)
  s <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, threshold = 0.6, target = "unobserved")
  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "outer band \\(plug-in\\)")
  expect_no_match(out, "naive")
  expect_match(out, "cap pi <= 1")
  expect_length(s$capped, 2)
  expect_true(all(s$capped >= 0 & s$capped <= 1))
  # a toenail-like case: r near one, rare events, the lower bound flat early
  n <- 400
  set.seed(5)
  R <- rbinom(n, 1, 0.92); y <- ifelse(R == 1, rbinom(n, 1, 0.04), NA)
  t <- sharp_scan_ci(rep(0.04, n), rep(0.92, n), R, y, target = "unobserved")
  expect_equal(t$held[["lower"]], qlogis(0.04) - qlogis(0.92 * 0.04), tolerance = 1e-8)
  expect_output(print(t), "held by the cap")
})

test_that("the default target is all rows, as for the mean", {
  x <- ci_sim(1000)
  a <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, lambda = c(0, 1))
  expect_identical(a$target, "all")
  expect_identical(a$n_target, nrow(x$d))
  expect_identical(formals(sharp_scan_ci)$target, formals(sharp_scan_mean_ci)$target)
})

test_that("debiased bounds are monotone in lambda even where the cap holds them", {
  # toenail-like: r near one, rare events; the raw debiased lower bound can
  # rise by sampling error once the cap takes over
  set.seed(9)
  n <- 300
  R <- rbinom(n, 1, 0.92); y <- ifelse(R == 1, rbinom(n, 1, 0.04), NA)
  eta <- runif(n, 0.02, 0.06)
  for (tg in c("all", "unobserved")) {
    s <- sharp_scan_ci(eta, rep(0.92, n), R, y, target = tg, lambda = seq(0, 2, by = 0.05))
    cv <- s$curve
    expect_true(all(diff(cv$lower) <= 0) && all(diff(cv$upper) >= 0))
    expect_true(all(diff(cv$ci_lower) <= 0) && all(diff(cv$ci_upper) >= 0))
  }
})

test_that("a segment restricts the unobserved target to rows chosen by X", {
  x <- ci_sim(2000)
  seg <- x$d$x1 > 0
  lam <- c(0, 0.5, 1)
  s <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, target = "unobserved",
                      segment = seg, lambda = lam)
  expect_identical(s$n_target, sum(seg & x$d$observed == 0))
  expect_match(s$target, "unobserved.*segment")
  # the estimate is that of the scan on the segment's rows alone
  i <- which(seg)
  sub <- sharp_scan_ci(x$inp$eta_obs[i], x$inp$r[i], x$d$observed[i], x$d$y[i],
                        target = "unobserved", lambda = lam)
  expect_equal(s$curve$lower, sub$curve$lower)
  expect_equal(s$curve$upper, sub$curve$upper)
  expect_true(all(s$curve$se_upper > 0))
  # an all-TRUE segment changes nothing
  a <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, target = "unobserved",
                      segment = rep(TRUE, nrow(x$d)), lambda = lam)
  u <- sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, target = "unobserved",
                      lambda = lam)
  expect_equal(a$curve, u$curve)
  expect_error(sharp_scan_ci(x$inp$eta_obs, x$inp$r, x$d$observed, x$d$y, target = seg,
                              segment = seg), "segment")
})
