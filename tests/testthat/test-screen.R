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

test_that("a ratio over an exclusion below its noise floor is marked unreadable", {
  d <- simulate_shadow(2000, tilt = -1, seed = 2)
  d$noise <- with_seed(105, rnorm(nrow(d)))
  sc <- shadow_screen(d, "y", "observed", c("x1", "x2"), c("z_shadow", "noise"))
  expect_equal(sc$ratio_readable, sc$exclusion_chisq_df >= 1)
  expect_identical(sc$ratio_readable, c(TRUE, FALSE))
  expect_output(print(sc), "violation at noise level.*noise")
  expect_false(any(grepl("noise level.*z_shadow", capture.output(print(sc)))))
  expect_false(any(grepl("not readable", capture.output(print(sc)))))
})

test_that("the screen print fits 80 columns and formats each number on its own", {
  d <- simulate_shadow(3000, tilt = -1, seed = 4)
  d$noise <- with_seed(9, rnorm(nrow(d)))
  sc <- shadow_screen(d, "y", "observed", c("x1", "x2"),
                      c("z_shadow", "z_policy", "z_score", "noise"))
  out <- capture.output(print(sc))
  expect_true(all(nchar(out) <= 80))
  expect_true(any(grepl("chi2/df", out)))
  # one line per candidate, with its verdict on the same line
  for (cand in sc$candidate) {
    ln <- grep(paste0("^ *", cand, " "), out, value = TRUE)
    expect_length(ln, 1)
    expect_match(ln, sc$verdict[sc$candidate == cand])
  }
  # no padding zeros from formatting a whole column alike
  expect_false(any(grepl("[0-9]\\.[0-9]*000\\b", out)))
})

test_that("internal helpers are not exported", {
  expect_false("verdict_text" %in% getNamespaceExports("shadowscan"))
})

# --- continuous outcome: ordinal relevance on outcome classes ---------------

cont <- local({
  d <- simulate_shadow(4000, tilt = -1, outcome = "continuous", seed = 21)
  d$y[d$observed == 0] <- NA
  d
})

test_that("a continuous outcome is screened by an ordinal model on outcome classes", {
  sc <- shadow_screen(cont, "y", "observed", c("x1", "x2"), c("z_shadow", "z_policy"))
  oc <- attr(sc, "outcome_coding")
  expect_identical(oc$model, "ordinal")
  expect_identical(oc$classes, 5L)
  expect_identical(sc$verdict[sc$candidate == "z_shadow"], "pass")
  # relevance by hand: polr on five classes of about equal mass of the
  # observed y, cut between distinct values as screen_classes() does
  obs <- cont[cont$observed == 1, ]
  yc <- screen_classes(obs$y, k = 5, support = range(obs$y))
  yc <- factor(yc, levels = levels(yc), ordered = TRUE)
  expect_equal(as.vector(table(yc)), c(323, 323, 323, 323, 321))
  f0 <- MASS::polr(yc ~ x1 + x2, data = obs)
  f1 <- MASS::polr(yc ~ x1 + x2 + z_shadow, data = obs)
  d0 <- -2 * sum(table(yc) * log(prop.table(table(yc))))
  expect_equal(sc$relevance[1], (deviance(f0) - deviance(f1)) / d0, tolerance = 1e-5)
  expect_equal(sc$relevance_chisq_df[1], deviance(f0) - deviance(f1), tolerance = 1e-5)
  expect_true(is.na(sc$epp_relevance[1]))
  expect_output(print(sc), "ordinal model on 5 classes")
})

test_that("outcome classes can be set and are checked", {
  sc <- shadow_screen(cont, "y", "observed", "x1", "z_shadow", outcome_classes = 3)
  expect_identical(attr(sc, "outcome_coding")$classes, 3L)
  expect_error(shadow_screen(cont, "y", "observed", "x1", "z_shadow", outcome_classes = 1),
               "outcome_classes")
  # no covariates: the reduced model is the null model
  sc0 <- shadow_screen(cont, "y", "observed", character(), "z_shadow")
  expect_gt(sc0$relevance, 0)
})

test_that("a binary outcome keeps the logistic screen", {
  d <- simulate_shadow(2000, tilt = -1, seed = 2)
  sc <- shadow_screen(d, "y", "observed", c("x1", "x2"), "z_shadow")
  expect_identical(attr(sc, "outcome_coding")$model, "glm")
  expect_no_match(paste(capture.output(print(sc)), collapse = ""), "ordinal")
})

test_that("a heaped outcome is classed between distinct values, never merged silently", {
  d <- cont
  set.seed(4)
  d$y <- ifelse(d$observed == 1, rpois(nrow(d), exp(0.2 + 0.4 * d$z_shadow)), NA)
  sc <- shadow_screen(d, "y", "observed", "x1", "z_shadow")
  yo <- d$y[d$observed == 1]
  # 0 and 1 hold most of the mass: each must stay a class of its own
  expect_gte(attr(sc, "outcome_coding")$classes, 4L)
  expect_output(print(sc), "classes")
  # few distinct values: the values themselves are the levels
  d$y <- ifelse(d$observed == 1, pmin(rpois(nrow(d), 1), 3), NA)
  sc3 <- shadow_screen(d, "y", "observed", "x1", "z_shadow")
  expect_identical(attr(sc3, "outcome_coding")$classes, 4L)
})

# --- categorical outcomes with more than two levels --------------------------

cat3 <- local({
  d <- simulate_shadow(4000, tilt = -1, outcome = "continuous", seed = 22)
  lv <- c("low", "mid", "high")
  d$y3 <- lv[findInterval(d$y, quantile(d$y[d$observed == 1], c(1/3, 2/3))) + 1]
  d$y3[d$observed == 0] <- NA
  d
})

test_that("an ordered factor outcome uses its levels in the ordinal model", {
  d <- cat3
  d$y3 <- factor(d$y3, levels = c("low", "mid", "high"), ordered = TRUE)
  sc <- shadow_screen(d, "y3", "observed", c("x1", "x2"), "z_shadow")
  oc <- attr(sc, "outcome_coding")
  expect_identical(oc$model, "ordinal")
  expect_identical(oc$classes, 3L)
  obs <- d[d$observed == 1, ]
  f0 <- MASS::polr(y3 ~ x1 + x2, data = obs); f1 <- MASS::polr(y3 ~ x1 + x2 + z_shadow, data = obs)
  d0 <- -2 * sum(table(obs$y3) * log(prop.table(table(obs$y3))))
  expect_equal(sc$relevance, (deviance(f0) - deviance(f1)) / d0, tolerance = 1e-5)
  expect_identical(sc$df, 1L)
})

test_that("an unordered factor or character outcome uses a multinomial model", {
  for (d in list(transform(cat3, y3 = factor(y3)), cat3)) {
    sc <- shadow_screen(d, "y3", "observed", c("x1", "x2"), "z_shadow")
    oc <- attr(sc, "outcome_coding")
    expect_identical(oc$model, "multinomial")
    expect_identical(oc$classes, 3L)
    expect_identical(sc$df, 2L)                 # (J - 1) * dim(Z)
    obs <- d[d$observed == 1, ]; obs$y3 <- factor(obs$y3)
    m0 <- nnet::multinom(y3 ~ x1 + x2, data = obs, trace = FALSE)
    m1 <- nnet::multinom(y3 ~ x1 + x2 + z_shadow, data = obs, trace = FALSE)
    d0 <- -2 * sum(table(obs$y3) * log(prop.table(table(obs$y3))))
    expect_equal(sc$relevance, (deviance(m0) - deviance(m1)) / d0, tolerance = 1e-4)
    expect_output(print(sc), "multinomial model on 3 classes")
  }
})

test_that("a factor with two observed levels stays binary", {
  d <- simulate_shadow(2000, tilt = -1, seed = 2)
  d$yf <- factor(ifelse(d$y == 1, "bad", "good"), levels = c("good", "bad", "unused"))
  d$yf[d$observed == 0] <- NA
  sc <- shadow_screen(d, "yf", "observed", c("x1", "x2"), "z_shadow")
  ref <- shadow_screen(d, "y", "observed", c("x1", "x2"), "z_shadow")
  expect_identical(attr(sc, "outcome_coding")$model, "glm")
  expect_equal(sc$relevance, ref$relevance)
})
