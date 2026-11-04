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
