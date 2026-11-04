test_that("plot methods run and return their input", {
  pdf(NULL)
  on.exit(dev.off())
  d <- simulate_shadow(2000, tilt = -1, seed = 31)
  rej <- d[d$observed == 0, ]
  for (thr in list(NULL, 0.3, 0.99)) {
    s <- shadow_scan(rej$eta_obs, rej$obs_rate, threshold = thr)
    expect_identical(plot(s), s)
  }
  s <- shadow_scan(rej$eta_obs, rej$obs_rate, aggregate = "union", threshold = 0.4)
  expect_identical(plot(s, scale = "log"), s)
  sc <- shadow_screen(d, "y", "observed", c("x1", "x2"),
                      c("z_shadow", "z_policy", "z_score"))
  expect_identical(plot(sc), sc)
})

test_that("screen plot handles both scales, zeros and many candidates", {
  pdf(NULL)
  on.exit(dev.off())
  sc <- data.frame(candidate = letters[1:14],
                   relevance = c(0, 10^-(1:13 / 3)), exclusion = 10^-(14:1 / 4),
                   stringsAsFactors = FALSE)
  sc$ratio <- sc$relevance / sc$exclusion
  class(sc) <- c("shadow_screen", "data.frame")
  expect_identical(plot(sc), sc)
  expect_identical(plot(sc, scale = "linear"), sc)
  expect_identical(plot(sc, label = c("a", "n", "zz")), sc)
  sc$relevance <- 0; sc$exclusion <- 0
  expect_error(plot(sc), "All values are zero")
})

test_that("ratio labels stay readable across magnitudes", {
  expect_identical(format_ratio(c(0.0042, 0.95, 16.52, 530.4, Inf)),
                   c("4.2e-03", "0.95", "16.52", "530", "Inf"))
})
