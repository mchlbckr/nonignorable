test_that("no value is split across classes and masses are near 1/k", {
  set.seed(41)
  x <- c(rep(10000, 300), rep(5000, 200), rep(35000, 80),
         round(runif(1420, 1000, 35000), -2))
  f <- screen_classes(x, k = 10)
  expect_true(all(tapply(as.integer(f), x, function(v) length(unique(v))) == 1))
  expect_lte(nlevels(f), 10)
  share <- as.vector(table(f)) / length(x)
  expect_true(all(share >= 0.05))
  # quantile classes would put a boundary on a heap
  expect_gt(attr(f, "heaping"), 0.1)
  # classes are ordered by value
  lo <- tapply(x, f, min); hi <- tapply(x, f, max)
  expect_true(all(hi[-length(hi)] < lo[-1]))
})

test_that("values outside the support get their own classes", {
  x <- c(600:661, 662:800, 801:820)
  f <- screen_classes(x, k = 4, support = c(662, 800))
  expect_identical(levels(f)[c(1, nlevels(f))], c("below", "above"))
  expect_true(all(f[x < 662] == "below"))
  expect_true(all(f[x > 800] == "above"))
  expect_equal(nlevels(f), 4 + 2)
})

test_that("few distinct values, NA and input checks", {
  f <- screen_classes(c(1, 2, 2, 3, NA), k = 10)
  expect_identical(levels(f), c("1", "2", "3"))
  expect_true(is.na(f[5]))
  expect_equal(attr(screen_classes(1:100, k = 10), "heaping"), 0)
  expect_error(screen_classes("a"), "numeric")
  expect_error(screen_classes(1:10, k = 0), ">= 1")
  expect_error(screen_classes(1:10, support = c(5, 1)), "increasing")
})

test_that("a hard cut-off in selection triggers the warning; classes catch it", {
  d <- simulate_shadow(6000, tilt = -1, seed = 42)
  d$observed[d$z_score < -0.5] <- 0L          # never accepted below a cut-off
  d$y[d$observed == 0] <- NA
  expect_warning(lin <- shadow_screen(d, "y", "observed", c("x1", "x2"), "z_score"),
                 "never seen among the observed rows")
  cls <- shadow_screen(d, "y", "observed", c("x1", "x2"), "z_score",
                       numeric = "classes")
  expect_gt(cls$exclusion, 1.5 * lin$exclusion)
  cd <- attr(cls, "coding")
  expect_setequal(cd$variable, c("x1", "x2", "z_score"))
  expect_true("below" %in% levels(screen_classes(d$z_score,
                support = range(d$z_score[d$observed == 1]))))
  expect_output(print(cls), "numeric as classes")
})

test_that("linear coding is unchanged and quiet without a cut-off", {
  d <- simulate_shadow(3000, seed = 43)
  expect_silent(sc <- shadow_screen(d, "y", "observed", c("x1", "x2"), "z_shadow"))
  expect_identical(attr(sc, "coding")$coding, rep("linear", 3))
})
