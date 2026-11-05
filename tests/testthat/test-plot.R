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

test_that("the screen plot labels unreadable ratios as such", {
  pdf(NULL)
  on.exit(dev.off())
  d <- simulate_shadow(2000, tilt = -1, seed = 2)
  d$noise <- with_seed(105, rnorm(nrow(d)))
  sc <- shadow_screen(d, "y", "observed", c("x1", "x2"), c("z_shadow", "noise"))
  txt <- character()
  local_mocked_bindings(place_labels = function(x, y, labels, ...) txt <<- labels)
  plot(sc)
  expect_match(txt[2], "noise \\(ratio not readable\\)")
  expect_match(txt[1], "z_shadow \\(ratio 17")
})

# The ACTG 175 screen of the applications vignette, as a hand-built object.
actg_like <- function() {
  sc <- data.frame(candidate = c("cd40", "cd420", "gender", "offtrt"),
                   relevance = c(0.138, 0.238, 5.69e-4, 0.0337),
                   exclusion = c(1.26e-3, 4.79e-3, 5.9e-9, 0.144),
                   stringsAsFactors = FALSE)
  sc$ratio <- sc$relevance / sc$exclusion
  sc$ratio_readable <- c(TRUE, TRUE, FALSE, TRUE)
  class(sc) <- c("shadow_screen", "data.frame")
  sc
}

seg_cross <- function(p1, p2, q1, q2) {          # proper intersection of two segments
  o <- function(a, b, c) sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
  o(p1, p2, q1) * o(p1, p2, q2) < 0 && o(q1, q2, p1) * o(q1, q2, p2) < 0
}
box_overlap <- function(a, b)
  min(a[2], b[2]) > max(a[1], b[1]) && min(a[4], b[4]) > max(a[3], b[3])

test_that("screen labels avoid the annotation and leaders do not cross the line", {
  png(tempfile(fileext = ".png"), width = 7, height = 4.5, units = "in", res = 96)
  on.exit(dev.off())                     # the device knitr draws vignettes on
  got <- NULL
  orig <- shadowscan:::place_labels
  local_mocked_bindings(place_labels = function(...) { got <<- orig(...); got })
  plot(actg_like())
  expect_false(is.null(got$obstacles))
  ann <- got$obstacles[[1]]
  for (b in got$boxes) expect_false(box_overlap(b, ann))
  for (l in got$leaders)
    expect_false(seg_cross(l[1:2], l[3:4], got$line[1:2], got$line[3:4]))
})

test_that("decade labels are thinned when they would collide", {
  # ten decades on 3 inches, labels 0.35 inches wide: every second one
  keep <- shadowscan:::thin_decades(-9:0, inch_per_decade = 0.3, label_width = 0.35)
  expect_equal(keep, seq(0, -9, by = -2)[order(seq(0, -9, by = -2))])
  expect_true(0 %in% keep)
  # ample room: all of them
  expect_equal(shadowscan:::thin_decades(-3:0, 1, 0.35), -3:0)
})

test_that("screen labels do not overlap and leaders do not run through them", {
  png(tempfile(fileext = ".png"), width = 7, height = 4.5, units = "in", res = 96)
  on.exit(dev.off())                     # the device knitr draws vignettes on
  got <- NULL
  orig <- shadowscan:::place_labels
  local_mocked_bindings(place_labels = function(...) { got <<- orig(...); got })
  plot(actg_like())
  through <- function(l, b) {
    t <- seq(0, 1, length.out = 50)
    x <- l[1] + t * (l[3] - l[1]); y <- l[2] + t * (l[4] - l[2])
    any(x > b[1] & x < b[2] & y > b[3] & y < b[4])
  }
  for (j in seq_along(got$boxes)) for (k in seq_along(got$boxes))
    if (j < k) expect_false(box_overlap(got$boxes[[j]], got$boxes[[k]]))
  for (l in got$leaders) {
    hits <- vapply(c(got$boxes, got$obstacles), function(b) through(l, b), TRUE)
    # a leader may touch only its own label, at its end point
    expect_lte(sum(hits), 0)
  }
})
