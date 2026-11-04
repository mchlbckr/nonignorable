screen_sim <- function(...) {
  d <- simulate_shadow(8000, tilt = -1, seed = 21)
  shadow_screen(d, "y", "observed", c("x1", "x2"),
                c("z_shadow", "z_policy", "z_score"), ...)
}

test_that("high exclusion and vanishing relevance is rejected", {
  sc <- screen_sim()
  p <- sc[sc$candidate == "z_policy", ]
  expect_identical(p$verdict, "reject")
  expect_gt(p$exclusion, 10 * p$relevance)
  # Under nonignorable selection z_policy picks up some relevance among the
  # observed through the collider R; under ignorable selection it has none
  # and is rejected on both counts.
  d0 <- simulate_shadow(8000, tilt = 0, seed = 21)
  p0 <- shadow_screen(d0, "y", "observed", c("x1", "x2"), "z_policy")
  expect_lt(p0$relevance_chisq_df, 10)
  expect_identical(p0$reason, "weak and violates exclusion")
})

test_that("valid shadow variable passes, the double-duty score does not", {
  sc <- screen_sim()
  expect_identical(sc$verdict[sc$candidate == "z_shadow"], "pass")
  sco <- sc[sc$candidate == "z_score", ]
  expect_gt(sco$relevance_chisq_df, 10)       # strongly relevant ...
  expect_identical(sco$verdict, "reject")     # ... but used in selection
  expect_identical(sco$reason, "violates exclusion")
})

test_that("measure is the incremental McFadden R^2 of the paper", {
  d <- simulate_shadow(3000, seed = 22)
  sc <- shadow_screen(d, "y", "observed", c("x1", "x2"), "z_shadow")
  obs <- d[d$observed == 1, ]
  m0 <- glm(y ~ x1 + x2, binomial, obs)
  m1 <- glm(y ~ x1 + x2 + z_shadow, binomial, obs)
  expect_equal(sc$relevance, (m0$deviance - m1$deviance) / m0$null.deviance)
  r0 <- glm(observed ~ x1 + x2, binomial, d)
  r1 <- glm(observed ~ x1 + x2 + z_shadow, binomial, d)
  expect_equal(sc$exclusion, (r0$deviance - r1$deviance) / r0$null.deviance)
  expect_equal(sc$relevance_chisq_df, m0$deviance - m1$deviance)
})

test_that("factor candidates, no covariates, logical response", {
  d <- simulate_shadow(3000, seed = 23)
  d$zf <- cut(d$z_shadow, c(-Inf, -1, 0, 1, Inf))
  d$obs_l <- d$observed == 1
  sc <- shadow_screen(d, "y", "obs_l", character(), "zf")
  expect_equal(sc$df, 3)
  expect_s3_class(sc, "shadow_screen")
})

test_that("inputs are validated", {
  d <- simulate_shadow(500, seed = 24)
  expect_error(shadow_screen(d, "y", "observed", "x1", "nope"), "not found")
  expect_error(shadow_screen(d, "y", "observed", "x1", "x1"), "also be")
  d$x0 <- 1
  expect_error(shadow_screen(d, "y", "observed", "x1", "x0"), "degrees")
})

test_that("too few minority cases per parameter are flagged, per side", {
  d <- simulate_shadow(400, tilt = -1, seed = 5)
  w <- character()
  sc <- withCallingHandlers(
    shadow_screen(d, "y", "observed", c("x1", "x2"), "z_shadow"),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_equal(sc$epp_relevance, 19 / 4)
  expect_equal(sc$epp_exclusion, 185 / 4)
  expect_true(any(grepl("outcome model", w)))
  expect_false(any(grepl("response model", w)))
})

test_that("ample minority cases stay quiet", {
  d <- simulate_shadow(5000, tilt = -1, seed = 1)
  expect_no_warning(shadow_screen(d, "y", "observed", c("x1", "x2"), "z_shadow"))
})
