test_that("stripped scans print and plot", {
  d <- simulate_shadow(2000, tilt = -1, seed = 5)
  rej <- d[d$observed == 0, ]
  s1 <- shadow_scan(rej$eta_obs, rej$obs_rate, threshold = 0.35)
  s2 <- shadow_scan_ci(pmin(pmax(d$eta_obs, 0.01), 0.99), d$obs_rate, d$observed, d$y,
                       threshold = 0.35)
  for (s in list(s1, s2)) {
    st <- strip_scan(s)
    expect_null(st$eta_obs); expect_null(st$r)
    expect_equal(st$curve, s$curve)
    expect_lt(as.numeric(object.size(st)), 20000)
    expect_output(print(st))
    pdf(NULL); expect_no_error(plot(st)); dev.off()
  }
})
