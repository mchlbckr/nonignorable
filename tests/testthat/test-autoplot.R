skip_if_not_installed("ggplot2")

# the data of the layer drawing geom `cls`, as given to the layer
layer_of <- function(p, cls) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, cls), TRUE))
  p$layers[i]
}
curve_layer <- function(p) {
  l <- layer_of(p, "GeomLine")
  expect_length(l, 1L)
  l[[1]]$data
}
texts_of <- function(p) unlist(lapply(layer_of(p, "GeomText"), function(l)
  as.character(if (is.null(l$aes_params$label)) l$data$label else l$aes_params$label)))
build_ok <- function(p) {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(b <- ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplot_gtable(b))
  invisible(b)
}

scan_rej <- function(...) {
  d <- simulate_shadow(2000, tilt = -1, seed = 31)
  rej <- d[d$observed == 0, ]
  sharp_scan(rej$eta_obs, rej$obs_rate, ...)
}

test_that("autoplot of a scan returns a ggplot whose curves are x$curve", {
  s <- scan_rej(threshold = 0.3)
  p <- ggplot2::autoplot(s)
  expect_s3_class(p, "ggplot")
  build_ok(p)
  cv <- curve_layer(p)
  for (sr in c("sharp", "naive")) for (e in c("lower", "upper")) {
    rows <- cv[cv$series == sr & cv$edge == e, ]
    col <- paste0(if (sr == "naive") "naive_" else "", e)
    expect_equal(rows$y, s$curve[[col]])
    expect_equal(rows$lambda, s$curve$lambda)
    expect_equal(rows$x, exp(s$curve$lambda))
  }
  expect_false("ci" %in% cv$series)
  band <- layer_of(p, "GeomRibbon")[[1]]$data
  expect_equal(band$lower, s$curve$lower)
  expect_equal(band$upper, s$curve$upper)
  # the MAR point, the threshold and the flip
  pts <- layer_of(p, "GeomPoint")
  mar <- pts[[length(pts)]]$data
  expect_equal(mar$y[mar$key == "point"], s$point)
  expect_equal(layer_of(p, "GeomHline")[[1]]$data$yintercept, 0.3)
  drop <- layer_of(p, "GeomSegment")[[1]]$data
  expect_equal(drop$x, exp(s$flip[["sharp"]]))
  txt <- texts_of(p)
  expect_true(any(grepl("flips at", txt)))
  expect_true("decision threshold" %in% txt)
})

test_that("the log scale puts lambda on the horizontal axis", {
  s <- scan_rej(threshold = 0.3)
  cv <- curve_layer(ggplot2::autoplot(s, scale = "log"))
  expect_equal(cv$x[cv$series == "sharp" & cv$edge == "lower"], s$curve$lambda)
  expect_error(ggplot2::autoplot(s, scale = "nope"))
})

test_that("a scan without threshold has no threshold layers", {
  p <- ggplot2::autoplot(scan_rej())
  build_ok(p)
  expect_length(layer_of(p, "GeomHline"), 0L)
  expect_length(layer_of(p, "GeomText"), 0L)
})

test_that("only the outer band flipping is labelled as such", {
  s <- scan_rej(threshold = 0.99)
  s$flip[["naive"]] <- 1      # as if only the outer band reached it
  s$flip[["sharp"]] <- NA
  p <- ggplot2::autoplot(s)
  build_ok(p)
  txt <- texts_of(p)
  expect_true(any(grepl("outer band flips at 2.7; sharp set never", txt, fixed = TRUE)))
})

test_that("edges held by the cap are marked and named in the legend", {
  d <- simulate_shadow(5000, tilt = -1, seed = 1)
  acc <- d[d$observed == 1, ]
  s <- sharp_scan(acc$eta_obs, acc$obs_rate, aggregate = "union",
                  lambda = seq(0, 6, length.out = 61))
  expect_false(anyNA(s$held))
  p <- ggplot2::autoplot(s)
  build_ok(p)
  pts <- layer_of(p, "GeomPoint")
  pts <- pts[[length(pts)]]$data
  held <- pts[startsWith(pts$key, "held_"), ]
  expect_setequal(held$key, c("held_lower", "held_upper"))
  expect_equal(held$x[held$key == "held_upper"], exp(s$held[["upper"]]))
  sh <- p$scales$get_scales("shape")
  expect_true(any(grepl("edge held by", vapply(sh$labels, deparse1, ""))))
})

test_that("CI, contrast and mean scans plot with their confidence sets", {
  d <- simulate_shadow(1500, tilt = -1, seed = 3)
  f <- ~ x1 + x2 + z_shadow
  inp <- crossfit_inputs(d, "observed",
    fit_eta = function(train, new) predict(glm(update(f, y ~ .), binomial, train),
                                           new, type = "response"),
    fit_r   = function(train, new) predict(glm(update(f, observed ~ .), binomial, train),
                                           new, type = "response"),
    seed = 1)
  s <- sharp_scan_ci(inp$eta_obs, inp$r, d$observed, d$y, target = "unobserved",
                     threshold = 0.35, lambda = seq(0, 2, length.out = 21))
  p <- ggplot2::autoplot(s)
  build_ok(p)
  cv <- curve_layer(p)
  expect_equal(cv$y[cv$series == "ci" & cv$edge == "upper"], s$curve$ci_upper)
  expect_equal(cv$y[cv$series == "ci" & cv$edge == "lower"], s$curve$ci_lower)
  expect_equal(p$labels$y, "Rate of the missing outcomes")

  k <- sharp_scan_contrast(inp$eta_obs, inp$r, d$observed, d$y,
                           target_a = d$x1 > 0, target_b = d$x1 <= 0,
                           lambda = seq(0, 2, length.out = 21), labels = c("A", "B"))
  p <- ggplot2::autoplot(k)
  build_ok(p)
  expect_equal(p$labels$y, "Difference of rates, A - B")
  txt <- texts_of(p)
  expect_true("no difference" %in% txt)

  dc <- simulate_shadow(1000, tilt = -1, outcome = "continuous", seed = 1)
  q <- obs_quantiles(dc, K = 50)
  rj <- dc$observed == 0
  m <- sharp_scan_mean(q[rj, ], dc$obs_rate[rj], threshold = 0.5)
  p <- ggplot2::autoplot(m, ylab = "my mean", main = "title")
  build_ok(p)
  expect_equal(curve_layer(p)$y[1:nrow(m$curve)], m$curve$lower)
  expect_equal(p$labels$y, "my mean")
  expect_equal(p$labels$title, "title")
})

test_that("autoplot of a distribution-function scan has its curves", {
  d <- simulate_shadow(800, tilt = -1, outcome = "continuous", seed = 4)
  xs <- c("x1", "x2", "z_shadow")
  fit <- function(lhs) function(train, newdata)
    predict(glm(reformulate(xs, lhs), binomial, train), newdata, type = "response")
  tg <- seq(-1, 1, by = 0.5)
  cf <- crossfit_cdf(d, "observed", "y", tg, fit_eta = fit("z"),
                     fit_r = fit("observed"), folds = 2, seed = 1)
  k <- sharp_scan_cdf(cf$eta_t, cf$r, d$observed, d$y, t = tg, lambda = 0.5)
  p <- ggplot2::autoplot(k)
  expect_s3_class(p, "ggplot")
  build_ok(p)
  cv <- curve_layer(p)
  expect_setequal(unique(cv$series), c("sharp", "naive", "ci", "point"))
  expect_equal(cv$y[cv$series == "point"], k$curve$point)
  expect_equal(cv$y[cv$series == "ci" & cv$edge == "lower"], k$curve$ci_lower)
  expect_equal(cv$t[cv$series == "sharp" & cv$edge == "upper"], k$curve$t)
})

# The ACTG 175 screen of the applications vignette, as a hand-built object.
screen_like <- function() {
  sc <- data.frame(candidate = c("cd40", "cd420", "gender", "offtrt"),
                   relevance = c(0.138, 0.238, 5.69e-4, 0.0337),
                   exclusion = c(1.26e-3, 4.79e-3, 5.9e-9, 0.144),
                   stringsAsFactors = FALSE)
  sc$ratio <- sc$relevance / sc$exclusion
  sc$ratio_readable <- c(TRUE, TRUE, FALSE, TRUE)
  sc$verdict <- c("pass", "pass", "reject", "reject")
  class(sc) <- c("shadow_screen", "data.frame")
  sc
}

test_that("autoplot of a screen styles points by verdict and readability", {
  sc <- screen_like()
  p <- ggplot2::autoplot(sc)
  expect_s3_class(p, "ggplot")
  build_ok(p)
  pts <- p$data
  expect_equal(pts$exclusion, sc$exclusion)
  expect_equal(pts$relevance, sc$relevance)
  # hollow at noise level, grey when rejected
  expect_equal(pts$fill, c(scan_colors[["sharp"]], scan_colors[["sharp"]], "white", "grey60"))
  expect_equal(pts$border, c("white", "white", "grey60", "white"))
  lab <- layer_of(p, "GeomText")
  labs <- lab[[length(lab)]]$data$label
  expect_match(labs[3], "gender (violation at noise level; rejected)", fixed = TRUE)
  expect_match(labs[1], "cd40 (ratio 110)", fixed = TRUE)
  expect_match(labs[4], "offtrt (ratio 0.23; rejected)", fixed = TRUE)
  # the equality line runs along the diagonal of the square panel
  eq <- layer_of(p, "GeomSegment")
  eq <- eq[[1]]$data
  expect_equal(eq$x, eq$y)
  expect_equal(eq$xend, eq$yend)
  expect_true(p$scales$get_scales("x")$trans$name %in% c("log-10", "log10"))
  expect_equal(p$theme$aspect.ratio, 1)
})

test_that("screen labels stay off the equality line and apart", {
  sc <- screen_like()
  p <- ggplot2::autoplot(sc)
  lab <- layer_of(p, "GeomText")
  lab <- lab[[length(lab)]]$data
  # offtrt sits below the line near the corner: its label drops below with a leader
  i <- lab$candidate == "offtrt"
  expect_lt(lab$y[i], lab$relevance[i])
  expect_false(is.na(lab$x0[i]))
  # without the vertical push two labels of the same column would overlap
  out <- screen_label_layout(c(0.2, 0.2), c(0.6, 0.61), c("aaaa", "bbbb"))
  expect_gte(abs(diff(out$ly)), 0.05 - 1e-9)
  expect_true(all(is.na(out$x0)))
})

test_that("screen autoplot: linear scale, zeros, label choice", {
  sc <- screen_like()
  p <- ggplot2::autoplot(sc, scale = "linear", label = "cd40")
  build_ok(p)
  lab <- layer_of(p, "GeomText")
  expect_equal(lab[[length(lab)]]$data$candidate, "cd40")
  sc$relevance[1] <- 0
  p <- ggplot2::autoplot(sc)
  build_ok(p)
  expect_gt(p$data$relevance[1], 0)        # drawn at the lower limit
  sc$relevance <- 0; sc$exclusion <- 0
  expect_error(ggplot2::autoplot(sc), "All values are zero")
  d <- simulate_shadow(2000, tilt = -1, seed = 31)
  real <- shadow_screen(d, "y", "observed", c("x1", "x2"),
                        c("z_shadow", "z_policy", "z_score"))
  build_ok(ggplot2::autoplot(real, label = character()))
})

test_that("without ggplot2 the methods say what is missing", {
  local_mocked_bindings(has_ggplot2 = function() FALSE)
  expect_error(autoplot.sharp_scan(scan_rej()), "needs the package 'ggplot2'")
  expect_error(autoplot.shadow_screen(screen_like()), "install.packages")
  expect_error(autoplot.sharp_scan_cdf(list()), "use plot\\(\\)")
})

test_that("the methods are registered for ggplot2::autoplot", {
  expect_true(is.function(getS3method("autoplot", "sharp_scan", envir = asNamespace("ggplot2"))))
  expect_true(is.function(getS3method("autoplot", "shadow_screen", envir = asNamespace("ggplot2"))))
  expect_true(is.function(getS3method("autoplot", "sharp_scan_cdf", envir = asNamespace("ggplot2"))))
})
