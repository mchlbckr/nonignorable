# ggplot2 versions of the plot() methods. ggplot2 is in Suggests only: the
# methods are registered lazily for ggplot2::autoplot (S3method(ggplot2::autoplot,
# ...) in NAMESPACE), and every ggplot2 function is called with `ggplot2::`.

has_ggplot2 <- function() requireNamespace("ggplot2", quietly = TRUE)

check_ggplot2 <- function(fun) {
  if (!has_ggplot2())
    stop("`", fun, "()` needs the package 'ggplot2'. Install it with ",
         "install.packages(\"ggplot2\"), or use plot() for the base graphics ",
         "version.", call. = FALSE)
}

# A light theme close to the base plots: axes at left and bottom, no grid,
# the legend below the panel.
autoplot_theme <- function(text_col) {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(colour = text_col),
      axis.text = ggplot2::element_text(colour = text_col),
      axis.line = ggplot2::element_line(colour = text_col),
      axis.ticks = ggplot2::element_line(colour = text_col),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.justification = "left",
      legend.key.width = grid_unit(2, "lines"),
      legend.background = ggplot2::element_blank())
}

grid_unit <- function(x, units) ggplot2::unit(x, units)

# The scan plots merge all their scales into one legend; it takes `ncol`
# columns, set on every scale so that the legends stay merged.
merged_guides <- function(ncol = 2L) {
  g <- ggplot2::guide_legend(ncol = ncol, byrow = TRUE)
  ggplot2::guides(fill = g, colour = g, linetype = g, linewidth = g, shape = g)
}

#' Sensitivity curve of a scan, with ggplot2
#'
#' The ggplot2 version of [plot.sharp_scan()], with the same information: the
#' sharp identified set as a filled band with solid edges, the outer band
#' (without the cap \eqn{\pi \le 1}) dashed, the confidence set dotted where
#' the scan has one, the point that a model on the observed rows asserts
#' (\eqn{\lambda = 0}), and, if the scan has a `threshold`, the decision
#' threshold with a drop line and a label where the set first reaches it.
#' Where the cap \eqn{\pi \le 1} holds an edge of the set flat from some
#' \eqn{\lambda} on (see `held` in [sharp_scan()]), a diamond on that edge
#' marks where, and the legend gives the value.
#'
#' The legend sits below the panel. The result is an ordinary ggplot object:
#' add layers, scales, `labs()` or a `theme()` to it as usual. The methods are
#' registered for [ggplot2::autoplot()] when ggplot2 is loaded; ggplot2 is not
#' required to install or use the rest of the package.
#'
#' The curve layers carry their data in long form: `layer_data()` or
#' `p$layers[[i]]$data` gives `lambda`, the horizontal position `x`, the value
#' `y`, the `series` (`"sharp"`, `"naive"`, `"ci"`) and the `edge`.
#'
#' @param object A `"sharp_scan"` object: from [sharp_scan()],
#'   [sharp_scan_mean()], [sharp_scan_ci()], [sharp_scan_mean_ci()],
#'   [sharp_scan_contrast()] or [sharp_scan_mean_contrast()].
#' @param scale `"ratio"` (default) puts \eqn{e^\lambda} on the horizontal
#'   axis, `"log"` puts \eqn{\lambda}, as in [plot.sharp_scan()].
#' @param colors Named colours `sharp` and `naive`, as in [plot.sharp_scan()].
#' @param text_col Colour for text and reference lines.
#' @param xlab,ylab,main Axis labels and title; `NULL` gives the same defaults
#'   as [plot.sharp_scan()]. `+ ggplot2::labs()` works too.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [plot.sharp_scan()] for base graphics.
#' @exportS3Method ggplot2::autoplot
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   d <- simulate_shadow(5000, tilt = -1, seed = 1)
#'   rej <- d[d$observed == 0, ]
#'   s <- sharp_scan(rej$eta_obs, rej$obs_rate, threshold = 0.35)
#'   autoplot(s)
#'   autoplot(s, scale = "log") + theme_minimal()
#' }
autoplot.sharp_scan <- function(object, scale = c("ratio", "log"),
                                colors = scan_colors, text_col = "grey20",
                                xlab = NULL, ylab = NULL, main = NULL, ...) {
  check_ggplot2("autoplot.sharp_scan")
  .data <- ggplot2::.data
  x <- object
  scale <- match.arg(scale)
  cv <- x$curve
  is_mean <- inherits(x, "sharp_scan_mean")
  is_contrast <- inherits(x, "sharp_scan_contrast")
  is_ci <- inherits(x, "sharp_scan_ci") || is_contrast
  if (is.null(ylab)) ylab <- if (isTRUE(x$unobserved) && !is_contrast)
      paste(if (is_mean) "Mean" else "Rate", "of the missing outcomes")
    else if (is_mean) "Portfolio mean"
    else if (is_contrast) sprintf("Difference of %s, %s - %s",
                                  if (isTRUE(x$of_means)) "means" else "rates",
                                  x$labels[1], x$labels[2])
    else "Portfolio rate"
  if (is.null(xlab)) xlab <- if (scale == "ratio")
    expression("Selection strength " * e^lambda * " (ratio of response propensities)")
  else expression("Selection strength " * lambda * " (log of the ratio)")
  tx <- function(l) if (scale == "ratio") exp(l) else l
  sym <- if (scale == "ratio") quote(e^lambda) else quote(lambda)
  has_ci <- !is.null(cv$ci_lower)
  has_point <- min(cv$lambda) == 0
  naive_col <- grDevices::adjustcolor(colors[["naive"]], 0.75)
  sharp_fill <- grDevices::adjustcolor(colors[["sharp"]], 0.18)

  # the curves in long form: one row per lambda, series and edge
  series <- c("sharp", "naive", if (has_ci) "ci")
  col_of <- c(sharp = "", naive = "naive_", ci = "ci_")
  lines <- do.call(rbind, lapply(series, function(s) do.call(rbind, lapply(
    c("lower", "upper"), function(e) data.frame(
      lambda = cv$lambda, x = tx(cv$lambda), y = cv[[paste0(col_of[[s]], e)]],
      series = s, edge = e, stringsAsFactors = FALSE)))))
  lines$group <- paste(lines$series, lines$edge)
  band <- data.frame(lambda = cv$lambda, x = tx(cv$lambda),
                     lower = cv$lower, upper = cv$upper)

  # the points: lambda = 0, and where the cap starts to hold an edge
  held <- x$held
  held <- if (is.null(held)) numeric() else held[!is.na(held)]
  pts <- data.frame(x = numeric(), y = numeric(), key = character(),
                    stringsAsFactors = FALSE)
  if (has_point)
    pts <- rbind(pts, data.frame(x = tx(0), y = x$point, key = "point"))
  for (side in names(held)) {
    yh <- stats::approx(cv$lambda, cv[[side]], held[[side]], rule = 2)$y
    pts <- rbind(pts, data.frame(x = tx(held[[side]]), y = yh,
                                 key = paste0("held_", side)))
  }

  # one legend for everything: each key has a value in every scale, NA or
  # "blank" where a layer has nothing to show for it
  keys <- c("sharp", "naive", if (has_ci) "ci", if (has_point) "point",
            paste0("held_", names(held)))
  lab_of <- list(
    sharp = if (is_ci) quote("Sharp identified set, estimated")
            else if (!is_mean && identical(x$aggregate, "common")) quote("Sharp identified set (common tilt)")
            else quote("Sharp identified set"),
    naive = if (is_ci) quote("Outer band, plug-in (without " * pi <= 1 * ")")
            else quote("Outer band (without " * pi <= 1 * ")"),
    ci = bquote(.(sprintf("%g%% confidence set", 100 * x$level))),
    point = quote("Model on observed rows " * (lambda == 0)))
  for (side in names(held))
    lab_of[[paste0("held_", side)]] <- bquote(
      .(paste(if (side == "lower") "Lower" else "Upper", "edge held by")) ~
        pi <= 1 ~ "from" ~ .(sym) %~~% .(sprintf("%.2f", tx(held[[side]]))))
  labs <- as.expression(unname(lab_of[keys]))
  val <- function(...) { v <- c(...); v <- v[keys]; names(v) <- keys; v }
  scale_args <- list(name = NULL, limits = keys, breaks = keys, labels = labs,
                     na.value = NA, drop = FALSE)
  scale_args_lt <- scale_args; scale_args_lt$na.value <- "blank"

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$x, ymin = .data$lower,
                                      ymax = .data$upper, fill = "sharp"),
                         data = band, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group,
                                    colour = .data$series, linetype = .data$series,
                                    linewidth = .data$series),
                       data = lines) +
    do.call(ggplot2::scale_fill_manual,
            c(list(values = val(sharp = sharp_fill)), scale_args)) +
    do.call(ggplot2::scale_colour_manual,
            c(list(values = val(sharp = colors[["sharp"]], naive = naive_col,
                                ci = colors[["sharp"]])), scale_args)) +
    do.call(ggplot2::scale_linetype_manual,
            c(list(values = val(sharp = "solid", naive = "dashed", ci = "dotted",
                                point = "blank", held_lower = "blank",
                                held_upper = "blank")), scale_args_lt)) +
    do.call(ggplot2::scale_discrete_manual,
            c(list(aesthetics = "linewidth",
                   values = val(sharp = 0.9, naive = 0.65, ci = 0.6, point = 0,
                                held_lower = 0, held_upper = 0)), scale_args))

  # threshold, drop lines at the flips, and their labels at the bottom
  rows <- 0
  if (!is.null(x$threshold)) {
    rng <- range(lines$x)
    hj <- function(v) if (tx(v) > rng[1] + 0.7 * diff(rng)) 1.05 else -0.05
    flip_lab <- function(prefix, v)
      bquote(.(prefix) ~ .(sym) %~~% .(sprintf("%.1f", tx(v))))
    fs <- x$flip[["sharp"]]; fn <- x$flip[["naive"]]
    fc <- if (is_ci && !is.null(x$flip[["ci"]])) x$flip[["ci"]] else NA_real_
    drops <- data.frame(v = numeric(), col = character(), size = numeric(),
                        stringsAsFactors = FALSE)
    texts <- list()
    if (!is.na(fs)) {
      drops <- rbind(drops, data.frame(v = fs, col = colors[["sharp"]], size = 3.2))
      texts <- c(texts, list(list(v = fs, lab = flip_lab("flips at", fs), row = 1)))
    } else if (!is.na(fn)) {
      drops <- rbind(drops, data.frame(v = fn, col = colors[["naive"]], size = 2.8))
      texts <- c(texts, list(list(
        v = fn, row = 1, lab = bquote(.(sprintf(
          "outer band flips at %.1f; sharp set never", tx(fn)))))))
    }
    if (!is.na(fc)) {
      # where the confidence set reaches the threshold: from there on the
      # data no longer rule out the other decision
      drops <- rbind(drops, data.frame(v = fc, col = colors[["sharp"]], size = 2.4))
      texts <- c(texts, list(list(v = fc, lab = flip_lab("confidence set at", fc),
                                  row = 2)))
    }
    rows <- if (length(texts)) max(vapply(texts, `[[`, 0, "row")) else 0
    drops$x <- tx(drops$v)
    drops$threshold <- rep(x$threshold, nrow(drops))
    p <- p +
      ggplot2::geom_hline(yintercept = x$threshold, linetype = "dashed",
                          colour = text_col, linewidth = 0.4) +
      ggplot2::annotate("text", x = Inf, y = x$threshold, hjust = 1.02, vjust = 1.6,
                        size = 3, colour = text_col,
                        label = if (is_contrast && x$threshold == 0) "no difference"
                                else "decision threshold")
    if (nrow(drops))
      p <- p +
        ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$x,
                                           y = -Inf, yend = .data$threshold),
                              data = drops, linetype = "dotted", colour = text_col,
                              linewidth = 0.4) +
        ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$threshold),
                            data = drops, shape = 21, fill = "white", stroke = 1,
                            colour = drops$col, size = drops$size)
    for (tt in texts)
      p <- p + ggplot2::annotate("text", x = tx(tt$v), y = -Inf,
                                 label = deparse1(tt$lab), parse = TRUE,
                                 hjust = hj(tt$v), vjust = -0.6 - 1.6 * (tt$row - 1),
                                 size = 3, colour = text_col)
  }

  p <- p +
    ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y, shape = .data$key),
                        data = pts, fill = colors[["sharp"]], colour = "white",
                        size = 3.2, stroke = 1.2) +
    do.call(ggplot2::scale_shape_manual,
            c(list(values = val(point = 21, held_lower = 23, held_upper = 23)),
              scale_args)) +
    # room at the bottom for the flip labels
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(
      mult = c(0.04 + 0.08 * rows, 0.04))) +
    ggplot2::labs(x = xlab, y = ylab, title = main) +
    merged_guides() +
    autoplot_theme(text_col)
  p
}

#' Distribution-function scan, with ggplot2
#'
#' The ggplot2 version of the `plot()` method of [sharp_scan_cdf()]: the sharp
#' set for \eqn{P(Y \le t)} at the scan's \eqn{\lambda} as a filled band with
#' solid edges, the outer band dashed, the confidence set dotted and the
#' estimate under MAR (\eqn{\lambda = 0}) as a line, all pointwise in \eqn{t}.
#' The curve layers carry `t`, `y`, `series` and `edge` in long form.
#'
#' @param object A `"sharp_scan_cdf"` object.
#' @inheritParams autoplot.sharp_scan
#' @param xlab,ylab,main Axis labels and title; `ylab = NULL` gives the same
#'   default as `plot()`.
#' @return A ggplot object.
#' @seealso [sharp_scan_cdf()] for `plot()`.
#' @exportS3Method ggplot2::autoplot
autoplot.sharp_scan_cdf <- function(object, colors = scan_colors, text_col = "grey20",
                                    xlab = "t", ylab = NULL, main = NULL, ...) {
  check_ggplot2("autoplot.sharp_scan_cdf")
  .data <- ggplot2::.data
  x <- object
  cv <- x$curve
  if (is.null(ylab))
    ylab <- if (startsWith(x$target, "unobserved"))
      quote(P(Y <= t ~ "|" ~ R == 0)) else quote(P(Y <= t))
  naive_col <- grDevices::adjustcolor(colors[["naive"]], 0.75)
  sharp_fill <- grDevices::adjustcolor(colors[["sharp"]], 0.18)
  cols <- list(sharp = c("lower", "upper"), naive = c("naive_lower", "naive_upper"),
               ci = c("ci_lower", "ci_upper"), point = "point")
  lines <- do.call(rbind, lapply(names(cols), function(s) do.call(rbind, lapply(
    cols[[s]], function(e) data.frame(
      t = cv$t, y = cv[[e]], series = s,
      edge = if (s == "point") "point" else sub(".*_", "", e),
      stringsAsFactors = FALSE)))))
  lines$group <- paste(lines$series, lines$edge)
  band <- data.frame(t = cv$t, lower = cv$lower, upper = cv$upper)
  keys <- names(cols)
  labs <- as.expression(list(
    bquote("Sharp set at" ~ e^lambda == .(format(exp(x$lambda), digits = 3)) *
             ", pointwise in" ~ t),
    quote("Outer band (without " * pi <= 1 * ")"),
    bquote(.(sprintf("%g%% confidence set, pointwise", 100 * x$level))),
    quote("Model on observed rows " * (lambda == 0))))
  scale_args <- list(name = NULL, limits = keys, breaks = keys, labels = labs,
                     na.value = NA)
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$t, ymin = .data$lower,
                                      ymax = .data$upper, fill = "sharp"),
                         data = band, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(x = .data$t, y = .data$y, group = .data$group,
                                    colour = .data$series, linetype = .data$series,
                                    linewidth = .data$series),
                       data = lines) +
    do.call(ggplot2::scale_fill_manual,
            c(list(values = c(sharp = sharp_fill, naive = NA, ci = NA, point = NA)),
              scale_args)) +
    do.call(ggplot2::scale_colour_manual,
            c(list(values = c(sharp = colors[["sharp"]], naive = naive_col,
                              ci = colors[["sharp"]], point = text_col)), scale_args)) +
    do.call(ggplot2::scale_linetype_manual,
            c(list(values = c(sharp = "solid", naive = "dashed", ci = "dotted",
                              point = "solid")), scale_args)) +
    do.call(ggplot2::scale_discrete_manual,
            c(list(aesthetics = "linewidth",
                   values = c(sharp = 0.9, naive = 0.65, ci = 0.6, point = 0.65)),
              scale_args)) +
    ggplot2::expand_limits(y = c(0, 1)) +
    ggplot2::labs(x = xlab, y = ylab, title = main) +
    merged_guides() +
    autoplot_theme(text_col)
}

#' Relevance against exclusion violation, with ggplot2
#'
#' The ggplot2 version of [plot.shadow_screen()]: exclusion violation across,
#' relevance up, both incremental McFadden \eqn{R^2}, on logarithmic axes by
#' default, with the equality line dashed and the panel square. Rejected
#' candidates are grey and labelled "rejected"; hollow points mark candidates
#' whose exclusion violation is at noise level (`ratio_readable` is `FALSE`),
#' whose label gives no ratio. Values of zero are drawn at the lower axis
#' limit. See [plot.shadow_screen()] for how to read the plot.
#'
#' Labels go beside their point, away from the equality line, and are spread
#' vertically where they would overlap. The rules assume a plot of roughly
#' five by five inches; for crowded screens use `label = character()` and add
#' your own labels, for example with `ggrepel::geom_text_repel()`, from the
#' point layer's data (`candidate`, `exclusion`, `relevance`, `label`).
#'
#' @param object A `"shadow_screen"` object.
#' @param scale `"log"` (default) or `"linear"` axes.
#' @param label Candidates to label. Default: all if there are at most
#'   fifteen, else the most relevant one and the one with the highest ratio.
#' @param col Point colour.
#' @param text_col Colour for text and reference lines.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [plot.shadow_screen()] for base graphics.
#' @exportS3Method ggplot2::autoplot
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   d <- simulate_shadow(5000, tilt = -1, seed = 1)
#'   sc <- shadow_screen(d, "y", "observed", c("x1", "x2"),
#'                       c("z_shadow", "z_policy", "z_score"))
#'   ggplot2::autoplot(sc)
#' }
autoplot.shadow_screen <- function(object, scale = c("log", "linear"), label = NULL,
                                   col = scan_colors[["sharp"]], text_col = "grey20",
                                   xlab = expression("Exclusion violation, incremental McFadden " * R^2),
                                   ylab = expression("Relevance, incremental McFadden " * R^2),
                                   main = NULL, ...) {
  check_ggplot2("autoplot.shadow_screen")
  .data <- ggplot2::.data
  x <- object
  scale <- match.arg(scale)
  if (is.null(label)) {
    label <- if (nrow(x) <= 15L) x$candidate else
      unique(c(x$candidate[which.max(x$relevance)],
               x$candidate[which.max(x$ratio)]))
  }
  ex <- x$exclusion; re <- x$relevance
  if (scale == "log") {
    pos <- c(ex, re)[c(ex, re) > 0]
    if (!length(pos)) stop("All values are zero; use scale = \"linear\".", call. = FALSE)
    lo <- 10^floor(log10(min(pos)) - 0.3)
    hi <- 10^ceiling(log10(max(pos)) + 0.05)
    ex <- pmax(ex, lo); re <- pmax(re, lo)
    lim <- c(lo, hi)
    frac <- function(v) (log10(v) - log10(lo)) / (log10(hi) - log10(lo))
    unfrac <- function(f) 10^(log10(lo) + f * (log10(hi) - log10(lo)))
  } else {
    lim <- c(0, max(re, ex) * 1.08)
    frac <- function(v) v / lim[2]
    unfrac <- function(f) f * lim[2]
  }
  readable <- if (is.null(x$ratio_readable)) rep(TRUE, nrow(x)) else x$ratio_readable
  verdict <- if (is.null(x$verdict)) rep("pass", nrow(x)) else x$verdict
  hue <- ifelse(verdict == "reject", "grey60", col)
  pts <- data.frame(candidate = x$candidate, exclusion = ex, relevance = re,
                    ratio = x$ratio, readable = readable, verdict = verdict,
                    fill = ifelse(readable, hue, "white"),
                    border = ifelse(readable, "white", hue),
                    stringsAsFactors = FALSE)
  pts$label <- sprintf("%s (%s%s)", x$candidate,
                       ifelse(readable, sprintf("ratio %s", format_ratio(x$ratio)),
                              "violation at noise level"),
                       ifelse(verdict == "reject", "; rejected", ""))
  lab <- pts[pts$candidate %in% label, , drop = FALSE]
  lab <- cbind(lab, screen_label_layout(frac(lab$exclusion), frac(lab$relevance),
                                        lab$label))
  lab$x <- unfrac(lab$lx); lab$y <- unfrac(lab$ly)
  leaders <- lab[!is.na(lab$x0), , drop = FALSE]
  leaders[c("x0", "y0", "x1", "y1")] <- lapply(leaders[c("x0", "y0", "x1", "y1")], unfrac)

  p <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$exclusion, y = .data$relevance)) +
    ggplot2::annotate("segment", x = lim[1], y = lim[1], xend = lim[2], yend = lim[2],
                      linetype = "dashed", colour = text_col, linewidth = 0.4) +
    ggplot2::annotate("text", x = unfrac(0.03), y = unfrac(0.97),
                      label = "more relevance than violation", hjust = 0, vjust = 1,
                      fontface = "italic", size = 2.9, colour = text_col) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x0, y = .data$y0,
                                       xend = .data$x1, yend = .data$y1),
                          data = leaders, colour = text_col, linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$fill, colour = .data$border),
                        shape = 21, size = 3.6, stroke = 1.2) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                       data = lab, hjust = 0, vjust = 0.5, size = 2.9,
                       colour = text_col) +
    ggplot2::labs(x = xlab, y = ylab, title = main) +
    autoplot_theme(text_col) +
    ggplot2::theme(aspect.ratio = 1)
  if (scale == "log") {
    dec <- seq(log10(lim[1]), log10(lim[2]))
    # ticks at every decade, labels at most about six
    keep <- dec %in% thin_decades(dec, 1, length(dec) / 7.2)
    labs <- parse(text = ifelse(keep, paste0("10^", dec), "''"))
    p <- p +
      ggplot2::scale_x_log10(limits = lim, breaks = 10^dec, labels = labs,
                             expand = ggplot2::expansion(mult = 0.02)) +
      ggplot2::scale_y_log10(limits = lim, breaks = 10^dec, labels = labs,
                             expand = ggplot2::expansion(mult = 0.02))
  } else {
    p <- p +
      ggplot2::scale_x_continuous(limits = lim, expand = ggplot2::expansion(mult = 0.02)) +
      ggplot2::scale_y_continuous(limits = lim, expand = ggplot2::expansion(mult = 0.02))
  }
  p + ggplot2::coord_cartesian(clip = "off")
}

# Where the screen labels go, in fractions of the axis range (fx, fy are the
# points): beside the point and away from the equality line, i.e. to the left
# above it and to the right below it; else to the other side if that is clear
# of the line; else above (or below) the point, as far as it takes to clear
# the line, with a leader. Labels whose boxes would overlap are then pushed
# apart vertically. `char_w` and `line_h` are the width of a character and the
# height of a line as fractions of the axis, for a panel of about four inches
# and 8 pt text. Returns the left end and centre height of each label (`lx`,
# `ly`) and the leader from the point to the label (`x0`, `y0`, `x1`, `y1`,
# `NA` where there is none).
screen_label_layout <- function(fx, fy, labels, char_w = 0.0115, line_h = 0.05,
                                gap = 0.025, right_edge = 1.12) {
  n <- length(fx)
  w <- nchar(labels) * char_w
  lx <- ly <- rep(NA_real_, n)
  for (k in seq_len(n)) {
    above <- fy[k] > fx[k]
    left <- c(fx[k] - gap - w[k], fx[k] - gap)
    right <- c(fx[k] + gap, fx[k] + gap + w[k])
    fits_left <- left[1] >= -0.02
    fits_right <- right[2] <= right_edge
    # beside the point, a label stays clear of the line if the line passes
    # its height outside the label's horizontal span
    clear <- function(span) span[1] > fy[k] + line_h / 2 || span[2] < fy[k] - line_h / 2
    if (above && fits_left || fits_left && clear(left)) {
      lx[k] <- left[1]; ly[k] <- fy[k]
    } else if (!above && fits_right || fits_right && clear(right)) {
      lx[k] <- right[1]; ly[k] <- fy[k]
    } else if (above) {
      lx[k] <- max(fx[k] - w[k] / 2, -0.02)
      ly[k] <- min(max(fy[k] + gap + line_h / 2, lx[k] + w[k] + line_h), 1.1)
    } else {
      lx[k] <- min(fx[k] + w[k] / 2, right_edge) - w[k]
      ly[k] <- max(min(fy[k] - gap - line_h / 2, lx[k] - line_h), -0.05)
    }
  }
  for (iter in seq_len(20)) {
    moved <- FALSE
    o <- order(ly)
    for (a in seq_along(o)[-1]) {
      i <- o[a - 1]; j <- o[a]
      if (min(lx[i] + w[i], lx[j] + w[j]) > max(lx[i], lx[j]) && ly[j] - ly[i] < line_h) {
        push <- (line_h - (ly[j] - ly[i])) / 2
        ly[i] <- ly[i] - push; ly[j] <- ly[j] + push
        moved <- TRUE
      }
    }
    if (!moved) break
  }
  # a leader where the label has moved off its point
  tx <- pmin(pmax(fx, lx), lx + w)
  ty <- pmin(pmax(fy, ly - line_h / 2), ly + line_h / 2)
  far <- sqrt((tx - fx)^2 + (ty - fy)^2) > 1.5 * gap
  data.frame(lx = lx, ly = ly,
             x0 = ifelse(far, fx, NA), y0 = ifelse(far, fy, NA),
             x1 = ifelse(far, tx, NA), y1 = ifelse(far, ty, NA))
}
