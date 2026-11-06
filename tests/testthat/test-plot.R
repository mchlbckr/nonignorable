test_that("plot methods run and return their input", {
  pdf(NULL)
  on.exit(dev.off())
  d <- simulate_shadow(2000, tilt = -1, seed = 31)
  rej <- d[d$observed == 0, ]
  for (thr in list(NULL, 0.3, 0.99)) {
    s <- sharp_scan(rej$eta_obs, rej$obs_rate, threshold = thr)
    expect_identical(plot(s), s)
  }
  s <- sharp_scan(rej$eta_obs, rej$obs_rate, aggregate = "union", threshold = 0.4)
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
  expect_match(txt[2], "noise \\(violation at noise level")
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
  orig <- nonignorable:::place_labels
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
  keep <- nonignorable:::thin_decades(-9:0, inch_per_decade = 0.3, label_width = 0.35)
  expect_equal(keep, seq(0, -9, by = -2)[order(seq(0, -9, by = -2))])
  expect_true(0 %in% keep)
  # ample room: all of them
  expect_equal(nonignorable:::thin_decades(-3:0, 1, 0.35), -3:0)
})

test_that("screen labels do not overlap and leaders do not run through them", {
  png(tempfile(fileext = ".png"), width = 7, height = 4.5, units = "in", res = 96)
  on.exit(dev.off())                     # the device knitr draws vignettes on
  got <- NULL
  orig <- nonignorable:::place_labels
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

test_that("the screen plot shows the verdict: rejected candidates grey and labelled", {
  pdf(NULL)
  on.exit(dev.off())
  d <- simulate_shadow(2000, tilt = -1, seed = 2)
  sc <- shadow_screen(d, "y", "observed", c("x1", "x2"), c("z_shadow", "z_policy"))
  expect_identical(sc$verdict, c("pass", "reject"))
  txt <- character()
  local_mocked_bindings(place_labels = function(x, y, labels, ...) txt <<- labels)
  plot(sc)
  expect_no_match(txt[1], "rejected")
  expect_match(txt[2], "rejected")
  st <- nonignorable:::screen_point_style(c("pass", "reject", "pass", "reject"),
                                        c(TRUE, TRUE, FALSE, FALSE), col = "blue")
  expect_equal(st$border[1:2], c("white", "white"))      # filled: border is surface
  expect_equal(st$bg[1], "blue"); expect_equal(st$bg[2], "grey60")
  expect_equal(st$border[3:4], c("blue", "grey60"))      # hollow: border carries colour
})

# --- scan plots: legend placement, cap note, band name ----------------------

capture_legend <- function(expr) {
  got <- NULL
  orig <- nonignorable:::draw_scan_legend
  local_mocked_bindings(draw_scan_legend = function(...) { got <<- list(args = list(...), rect = orig(...)); got$rect })
  force(expr)
  got$usr <- par("usr")
  got
}

legend_hits_ink <- function(x, rect, scale = "ratio") {
  cv <- x$curve
  xs <- if (scale == "ratio") exp(cv$lambda) else cv$lambda
  cols <- intersect(c("lower", "upper", "naive_lower", "naive_upper", "ci_lower", "ci_upper"), names(cv))
  m <- as.matrix(cv[cols])
  lo <- apply(m, 1, min); hi <- apply(m, 1, max)
  b <- c(rect$left, rect$left + rect$w, rect$top - rect$h, rect$top)
  xo <- seq(b[1], b[2], length.out = 200)
  xo <- xo[xo >= min(xs) & xo <= max(xs)]
  l <- approx(xs, lo, xo)$y; h <- approx(xs, hi, xo)$y
  any(h >= b[3] & l <= b[4])
}

test_that("the scan legend avoids every curve", {
  png(tempfile(fileext = ".png"), width = 7, height = 4.5, units = "in", res = 96)
  on.exit(dev.off())
  d <- simulate_shadow(3000, tilt = -1, seed = 31)
  rej <- d[d$observed == 0, ]
  for (thr in list(NULL, 0.3)) {
    s <- sharp_scan(rej$eta_obs, rej$obs_rate, threshold = thr)
    got <- capture_legend(plot(s))
    expect_false(legend_hits_ink(s, got$rect$rect))
  }
  # a band that fills the top left: the legend has to go elsewhere
  s <- sharp_scan(rep(0.5, 50), rep(0.2, 50))
  got <- capture_legend(plot(s))
  expect_false(legend_hits_ink(s, got$rect$rect))
})

test_that("the legend spot prefers a free corner and grows the top when none is free", {
  usr <- c(0, 10, 0, 1)
  xs <- 0:10
  # ink high on the left, low on the right: top right and bottom left are free
  sp <- nonignorable:::scan_legend_spot(3, 0.3, usr, xs, lo = seq(0.6, 0.0, length.out = 11),
                                      hi = seq(0.99, 0.4, length.out = 11))
  expect_identical(sp$corner, "topright")
  # ink everywhere from 0.05 to 0.95: no corner free, the top has to grow
  sp <- nonignorable:::scan_legend_spot(3, 0.3, usr, xs, lo = rep(0.05, 11), hi = rep(0.95, 11))
  expect_identical(sp$corner, "topleft")
  expect_gt(sp$grow_to, 1)
  # an obstacle (the threshold line) rules a corner out
  sp <- nonignorable:::scan_legend_spot(3, 0.3, usr, xs, lo = rep(0.4, 11), hi = rep(0.5, 11),
                                      obstacles = list(c(0, 10, 0.85, 0.85)))
  expect_false(sp$corner %in% c("topleft", "topright"))
})

test_that("plots name the outer band, never a naive one", {
  pdf(NULL)
  on.exit(dev.off())
  d <- simulate_shadow(2000, tilt = -1, seed = 31)
  rej <- d[d$observed == 0, ]
  s <- sharp_scan(rej$eta_obs, rej$obs_rate, threshold = 0.99)
  got <- capture_legend(plot(s))
  labs <- vapply(as.list(got$args$legend), function(e) paste(deparse(e), collapse = ""), "")
  expect_true(any(grepl("Outer band", labs)))
  expect_false(any(grepl("Naive", labs)))
})

test_that("a bound held by the cap is noted in print and legend", {
  pdf(NULL)
  on.exit(dev.off())
  # r = 0.92 everywhere: the lower bound r * eta binds from
  # lambda = logit(eta) - logit(r * eta), here about 0.086
  s <- sharp_scan(rep(0.04, 100), rep(0.92, 100))
  expect_equal(s$held[["lower"]], qlogis(0.04) - qlogis(0.92 * 0.04), tolerance = 1e-8)
  expect_equal(s$held[["upper"]], qlogis(0.92 * 0.04 + 0.08) - qlogis(0.04), tolerance = 1e-8)
  expect_equal(s$capped[["lower"]], 1)
  expect_output(print(s), "held by the cap")
  got <- capture_legend(plot(s))
  labs <- vapply(as.list(got$args$legend), function(e) paste(deparse(e), collapse = ""), "")
  expect_true(any(grepl("held by", labs)))
  # nothing held: no note
  d <- simulate_shadow(2000, tilt = -1, seed = 31)
  rej <- d[d$observed == 0, ]
  s2 <- sharp_scan(rej$eta_obs, rej$obs_rate)
  expect_true(all(is.na(s2$held)))
  expect_no_match(paste(capture.output(print(s2)), collapse = "\n"), "held by")
})

test_that("flip labels get a free strip below curves that run low", {
  usr <- c(0, 10, 0, 1); xs <- 0:10
  box <- list(c(2, 4, 0, 0.05))
  # ink at 0.5 under the label: already clear
  expect_null(nonignorable:::flip_label_floor(box, usr, xs, rep(0.5, 11)))
  # ink at 0.02 under the label: the limit moves down until the label clears it
  lo <- nonignorable:::flip_label_floor(box, usr, xs, rep(0.02, 11))
  expect_lt(lo, 0)
  f <- 0.05 / 1
  expect_equal(0.02, lo + (f + 0.04) * (1 - lo), tolerance = 1e-12)
  # no label: nothing to do
  expect_null(nonignorable:::flip_label_floor(list(), usr, xs, rep(0, 11)))
})

test_that("a scan of unobserved rows is labelled as the missing outcomes", {
  d <- simulate_shadow(800, tilt = -1, seed = 2)
  rej <- d[d$observed == 0, ]
  s <- sharp_scan(rej$eta_obs, rej$obs_rate, threshold = 0.35, unobserved = TRUE,
                  lambda = seq(0, 2, length.out = 21))
  labs <- character()
  local_mocked_bindings(title = function(..., ylab = NULL) labs <<- c(labs, ylab),
                        .package = "graphics")
  pdf(NULL); on.exit(dev.off())
  plot(s)
  expect_true("Rate of the missing outcomes" %in% labs)
})
