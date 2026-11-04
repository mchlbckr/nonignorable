test_that("simulator is deterministic and leaves the RNG alone", {
  set.seed(99); before <- runif(1); set.seed(99)
  a <- simulate_shadow(500, seed = 7)
  after <- runif(1)
  expect_identical(before, after)
  b <- simulate_shadow(500, seed = 7)
  expect_identical(a, b)
})

test_that("simulator truth is internally consistent", {
  d <- simulate_shadow(2000, tilt = -1, seed = 8)
  expect_equal(attr(d, "tilt"), -1)
  # logit eta_obs = logit eta + tilt, exactly
  expect_equal(qlogis(d$eta_obs), qlogis(d$eta) - 1, tolerance = 1e-10)
  expect_true(all(is.na(d$y) == (d$observed == 0)))
  expect_true(all(d$obs_rate > 0 & d$obs_rate <= 1))
})
