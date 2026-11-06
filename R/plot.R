scan_colors <- c(sharp = "#2a78d6", naive = "#eb6834")

#' Sensitivity curve of a scan
#'
#' Plots the portfolio rate (or mean, for [sharp_scan_mean()]) against
#' selection strength: the sharp identified set as a filled band with solid
#' edges, the outer band, the bounds without the cap \eqn{\pi \le 1} (for a
#' rate \eqn{\mathrm{expit}(\ell \pm \lambda)}), dashed and lighter, the point that a
#' model on the observed rows asserts (\eqn{\lambda = 0}) marked, and, if the
#' scan has a `threshold`, the decision threshold with a drop line where the
#' sharp set first reaches it. That crossing is the reason to run the scan: it
#' shows at which strength of selection the decision flips.
#'
#' Where \eqn{\hat r} is small, as for rejected applications, the upper edge of
#' the sharp set coincides with the outer one; see [sharp_scan()].
#'
#' The legend goes into the first corner, from top left, that no curve, band
#' or annotation touches; if none is free, the vertical axis is extended
#' upwards to make room. Where the cap \eqn{\pi \le 1} binds in every target
#' row from some \eqn{\lambda} on, that edge of the set stays flat from there,
#' and the legend says from where (see `held` in [sharp_scan()]).
#'
#' @param x A `"sharp_scan"` or `"sharp_scan_mean"` object.
#' @param scale `"ratio"` (default) puts \eqn{e^\lambda}, the largest
#'   admissible ratio of response propensities, on the horizontal axis,
#'   starting at 1; `"log"` puts \eqn{\lambda}.
#' @param colors Named colours `sharp` and `naive`. The defaults are validated
#'   for colour-vision deficiency; on a dark background use
#'   `c(sharp = "#3987e5", naive = "#d95926")`.
#' @param text_col Colour for all text and reference lines. Text is never drawn
#'   in a series colour.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Passed to [graphics::plot.default()].
#' @return `x`, invisibly.
#' @export
#' @examples
#' d <- simulate_shadow(5000, tilt = -1, seed = 1)
#' rej <- d[d$observed == 0, ]
#' plot(sharp_scan(rej$eta_obs, rej$obs_rate, threshold = 0.35))
plot.sharp_scan <- function(x, scale = c("ratio", "log"),
                             colors = scan_colors, text_col = "grey20",
                             xlab = NULL, ylab = NULL, main = NULL,
                             ...) {
  scale <- match.arg(scale)
  cv <- x$curve
  is_mean <- inherits(x, "sharp_scan_mean")
  is_contrast <- inherits(x, "sharp_scan_contrast")
  is_ci <- inherits(x, "sharp_scan_ci") || inherits(x, "sharp_scan_mean_ci") || is_contrast
  if (is.null(ylab)) ylab <- if (isTRUE(x$unobserved) && !is_contrast)
      paste(if (is_mean) "Mean" else "Rate", "of the missing outcomes")
    else if (is_mean) "Portfolio mean"
    else if (is_contrast) sprintf("Difference of %s, %s - %s",
                                  if (isTRUE(x$of_means)) "means" else "rates",
                                  x$labels[1], x$labels[2])
    else "Portfolio rate"
  tx <- function(l) if (scale == "ratio") exp(l) else l
  xs <- tx(cv$lambda)
  if (is.null(xlab)) xlab <- if (scale == "ratio")
    expression("Selection strength " * e^lambda * " (ratio of response propensities)")
  else expression("Selection strength " * lambda * " (log of the ratio)")
  ylim <- range(cv$naive_lower, cv$naive_upper, cv$ci_lower, cv$ci_upper,
                x$threshold)
  ylim <- ylim + c(-0.04, 0.04) * diff(ylim)
  naive_col <- grDevices::adjustcolor(colors[["naive"]], 0.75)
  has_ci <- !is.null(cv$ci_lower)
  has_point <- min(cv$lambda) == 0
  sym <- if (scale == "ratio") "e^lambda" else "lambda"
  sym_e <- parse(text = sym)[[1]]

  # the legend, assembled first: its size decides where it goes
  held <- x$held
  held_items <- if (!is.null(held)) lapply(names(held)[!is.na(held)], function(side)
    list(lab = bquote(.(paste(if (side == "lower") "Lower" else "Upper",
                              "edge held by")) ~ pi <= 1 ~ "from" ~ .(sym_e) %~~%
                        .(sprintf("%.2f", tx(held[[side]])))),
         fill = NA, border = NA, lty = NA, lwd = NA, pch = NA, pt.bg = NA, col = NA))
  items <- c(list(
    list(lab = if (is_ci) quote("Sharp identified set, estimated")
               else if (!is_mean && identical(x$aggregate, "common")) quote("Sharp identified set (common tilt)")
               else quote("Sharp identified set"),
         fill = grDevices::adjustcolor(colors[["sharp"]], 0.18), border = colors[["sharp"]],
         lty = NA, lwd = NA, pch = NA, pt.bg = NA, col = NA),
    list(lab = if (is_ci) quote("Outer band, plug-in (without " * pi <= 1 * ")")
               else quote("Outer band (without " * pi <= 1 * ")"),
         fill = NA, border = NA, lty = 2, lwd = 1.5, pch = NA, pt.bg = NA, col = naive_col),
    if (has_ci) list(lab = bquote(.(sprintf("%g%% confidence set", 100 * x$level))),
                     fill = NA, border = NA, lty = 3, lwd = 1.3, pch = NA, pt.bg = NA,
                     col = colors[["sharp"]]),
    if (has_point) list(lab = quote("Model on observed rows " * (lambda == 0)),
                        fill = NA, border = NA, lty = NA, lwd = 2, pch = 21,
                        pt.bg = colors[["sharp"]], col = surface_col())
  ), held_items)
  items <- Filter(Negate(is.null), items)
  gc_ <- function(nm) unlist(lapply(items, function(it) if (is.na(it[[nm]])) NA else it[[nm]]))
  leg_args <- list(bty = "n", text.col = text_col, cex = 0.85,
                   legend = as.expression(lapply(items, `[[`, "lab")),
                   fill = gc_("fill"), border = gc_("border"),
                   lty = gc_("lty"), lwd = gc_("lwd"), pch = gc_("pch"),
                   pt.cex = 1.2, pt.bg = gc_("pt.bg"), col = gc_("col"),
                   merge = TRUE, x.intersp = 0.8)

  graphics::plot(xs, cv$upper, type = "n", ylim = ylim, xlab = "", ylab = "",
                 axes = FALSE, ...)

  # Annotations of the threshold, as functions of the user coordinates: they
  # are obstacles for the legend, and drawn after the coordinates are final.
  flip_adj <- function(v, pad, usr) {
    if (tx(v) > usr[1] + 0.7 * (usr[2] - usr[1])) c(1 + pad, -0.6) else c(-pad, -0.6)
  }
  annotations <- function(usr) {
    if (is.null(x$threshold)) return(list())
    lab <- function(v) bquote(.(as.name("flips at")) ~ .(sym_e) %~~% .(sprintf("%.1f", tx(v))))
    fs <- x$flip[["sharp"]]; fn <- x$flip[["naive"]]
    fc <- if (is_ci && !is.null(x$flip[["ci"]])) x$flip[["ci"]] else NA_real_
    out <- list(list(type = "hline"),
                list(type = "text", x = usr[2], y = x$threshold,
                     label = if (is_contrast && x$threshold == 0) "no difference" else "decision threshold",
                     adj = c(1, 1.6), cex = 0.8))
    if (!is.na(fs)) {
      out <- c(out, list(list(type = "drop", v = fs, cex = 1.3, col = colors[["sharp"]]),
                         list(type = "text", x = tx(fs), y = usr[3], label = lab(fs),
                              adj = flip_adj(fs, 0.08, usr), cex = 0.85)))
    } else if (!is.na(fn)) {
      out <- c(out, list(list(type = "drop", v = fn, cex = 1.1, col = colors[["naive"]]),
                         list(type = "text", x = tx(fn), y = usr[3],
                              label = sprintf("outer band flips at %.1f; sharp set never", tx(fn)),
                              adj = flip_adj(fn, 0.05, usr), cex = 0.8)))
    }
    if (!is.na(fc)) {
      # where the confidence set reaches the threshold: from there on the data
      # no longer rule out the other decision
      out <- c(out, list(list(type = "drop", v = fc, cex = 0.9, lwd = 1.5, col = colors[["sharp"]]),
                         list(type = "text", x = tx(fc), y = usr[3],
                              label = bquote("confidence set at" ~ .(sym_e) %~~% .(sprintf("%.1f", tx(fc)))),
                              adj = c(flip_adj(fc, 0.08, usr)[1], -2.2), cex = 0.8)))
    }
    out
  }
  obstacle_boxes <- function(ann, usr) lapply(ann, function(a) switch(a$type,
    hline = c(usr[1], usr[2], x$threshold, x$threshold),
    drop = c(tx(a$v), tx(a$v), usr[3], x$threshold),
    text = {
      w <- graphics::strwidth(a$label, cex = a$cex); h <- graphics::strheight(a$label, cex = a$cex)
      x0 <- a$x - a$adj[1] * w; y0 <- a$y - a$adj[2] * h
      c(x0, x0 + w, y0, y0 + h)
    }))

  ink <- as.matrix(cv[intersect(c("lower", "upper", "naive_lower", "naive_upper",
                                  "ci_lower", "ci_upper"), names(cv))])
  lo <- apply(ink, 1, min); hi <- apply(ink, 1, max)
  if (!is.null(x$threshold)) { lo <- pmin(lo, x$threshold); hi <- pmax(hi, x$threshold) }
  usr <- graphics::par("usr")
  # The flip labels sit at the bottom of the plot. Where a curve runs that low
  # beneath them (a lower bound near 0, say), lower the y limit until the
  # labels have a free strip below the ink.
  new_lo <- flip_label_floor(obstacle_boxes(Filter(function(a) a$type == "text" &&
                                                     a$y == usr[3], annotations(usr)), usr),
                             usr, xs, lo)
  if (!is.null(new_lo)) {
    graphics::plot.window(xlim = usr[1:2], ylim = c(new_lo, usr[4]), xaxs = "i", yaxs = "i")
    usr <- graphics::par("usr")
  }
  size <- do.call(graphics::legend, c(list(x = "topleft", plot = FALSE), leg_args))$rect
  spot <- scan_legend_spot(size$w, size$h, usr, xs, lo, hi,
                           obstacles = obstacle_boxes(annotations(usr), usr))
  if (!is.null(spot$grow_to)) {
    graphics::plot.window(xlim = usr[1:2], ylim = c(usr[3], spot$grow_to),
                          xaxs = "i", yaxs = "i")
    usr <- graphics::par("usr")
  }

  graphics::axis(1, col = text_col, col.axis = text_col)
  graphics::axis(2, las = 1, col = text_col, col.axis = text_col)
  graphics::box(bty = "l", col = text_col)
  graphics::title(xlab = xlab, ylab = ylab, main = main, col.lab = text_col,
                  col.main = text_col)
  graphics::polygon(c(xs, rev(xs)), c(cv$lower, rev(cv$upper)),
                    col = grDevices::adjustcolor(colors[["sharp"]], 0.18),
                    border = NA)
  graphics::lines(xs, cv$naive_upper, lty = 2, lwd = 1.5, col = naive_col)
  graphics::lines(xs, cv$naive_lower, lty = 2, lwd = 1.5, col = naive_col)
  graphics::lines(xs, cv$upper, lwd = 2, col = colors[["sharp"]])
  graphics::lines(xs, cv$lower, lwd = 2, col = colors[["sharp"]])
  if (has_ci) {
    graphics::lines(xs, cv$ci_upper, lwd = 1.3, lty = 3, col = colors[["sharp"]])
    graphics::lines(xs, cv$ci_lower, lwd = 1.3, lty = 3, col = colors[["sharp"]])
  }
  # The lambda = 0 point is explained in the legend, not by a label next to it:
  # when the decision flips near lambda = 0, the flip marker, its drop line and
  # the threshold all sit beside that point and any label there collides.
  if (has_point)
    graphics::points(tx(0), x$point, pch = 21, cex = 1.4, lwd = 2,
                     bg = colors[["sharp"]], col = surface_col())
  for (a in annotations(usr)) switch(a$type,
    hline = graphics::abline(h = x$threshold, lty = 2, col = text_col),
    drop = {
      graphics::segments(tx(a$v), usr[3], tx(a$v), x$threshold, lty = 3, col = text_col)
      graphics::points(tx(a$v), x$threshold, pch = 21, cex = a$cex,
                       lwd = if (is.null(a$lwd)) 2 else a$lwd,
                       bg = surface_col(), col = a$col)
    },
    text = graphics::text(a$x, a$y, a$label, adj = a$adj, col = text_col, cex = a$cex))
  inset <- c(0.02, 0.02)
  do.call(draw_scan_legend, c(list(x = spot$corner, inset = inset), leg_args))
  invisible(x)
}

draw_scan_legend <- function(...) graphics::legend(...)

# The lower y limit that gives the flip labels (boxes c(x0, x1, y0, y1), drawn
# from the bottom edge) a free strip beneath the lowest ink `lo` within their
# horizontal span, or NULL if they are already clear. Label heights are fixed
# in inches, so they keep their share f of the plot height when the limit moves.
flip_label_floor <- function(boxes, usr, xs, lo, pad = 0.04) {
  if (!length(boxes)) return(NULL)
  b <- do.call(rbind, boxes)
  x0 <- max(min(b[, 1]), min(xs)); x1 <- min(max(b[, 2]), max(xs))
  if (x0 > x1) return(NULL)
  xo <- seq(x0, x1, length.out = 200)
  ink <- min(stats::approx(xs, lo, xo)$y)
  H <- usr[4] - usr[3]
  f <- (max(b[, 4]) - usr[3]) / H
  if (ink >= max(b[, 4]) + pad * H) return(NULL)
  (ink - (f + pad) * usr[4]) / (1 - f - pad)
}

# Where the scan legend goes: the first corner whose box no curve, band or
# annotation touches. `lo`, `hi` are the lowest and highest ink at each `xs`;
# `obstacles` are boxes c(x0, x1, y0, y1). If no corner is free, the legend
# goes top left and `grow_to` is the upper y limit that lifts it clear of the
# ink beneath it.
scan_legend_spot <- function(lw, lh, usr, xs, lo, hi, obstacles = list(), pad = 0.02) {
  dx <- pad * diff(usr[1:2]); dy <- pad * diff(usr[3:4])
  boxes <- list(
    topleft     = c(usr[1] + dx, usr[1] + dx + lw, usr[4] - dy - lh, usr[4] - dy),
    topright    = c(usr[2] - dx - lw, usr[2] - dx, usr[4] - dy - lh, usr[4] - dy),
    bottomleft  = c(usr[1] + dx, usr[1] + dx + lw, usr[3] + dy, usr[3] + dy + lh),
    bottomright = c(usr[2] - dx - lw, usr[2] - dx, usr[3] + dy, usr[3] + dy + lh))
  span <- function(b) {
    xo <- seq(b[1], b[2], length.out = 200)
    xo <- xo[xo >= min(xs) & xo <= max(xs)]
    list(lo = stats::approx(xs, lo, xo)$y, hi = stats::approx(xs, hi, xo)$y)
  }
  overlap <- function(a, b)
    min(a[2], b[2]) >= max(a[1], b[1]) && min(a[4], b[4]) >= max(a[3], b[3])
  free <- function(b) {
    s <- span(b)
    !any(s$hi >= b[3] & s$lo <= b[4]) &&
      !any(vapply(obstacles, overlap, TRUE, b = b))
  }
  for (nm in names(boxes)) if (free(boxes[[nm]])) return(list(corner = nm, box = boxes[[nm]]))
  # lift the top: the legend's share of the height stays fixed, the ink under
  # it has to end below its lower edge
  b <- boxes$topleft
  top <- max(span(b)$hi, vapply(obstacles, function(o)
    if (o[2] >= b[1] && o[1] <= b[2]) o[4] else -Inf, 0))
  f <- (lh + 2 * dy) / diff(usr[3:4])
  grow_to <- (top - f * usr[3]) / (1 - f)
  list(corner = "topleft", box = NULL, grow_to = grow_to + 0.01 * (grow_to - usr[3]))
}

#' Relevance against exclusion violation
#'
#' Scatter of the screen: exclusion violation on the horizontal axis,
#' relevance on the vertical, both incremental McFadden \eqn{R^2} on the same
#' scale, with the equality line dashed. Above the line a candidate predicts
#' the outcome better than it predicts selection. The line is not the border
#' between usable and unusable: under strong selection on the outcome a
#' validly excluded candidate shows exclusion signal through the outcome and
#' sits close to it (see [shadow_screen()]). Points are labelled directly
#' rather than through a legend; labels are placed so that they do not overlap
#' each other or the points.
#'
#' Both measures typically span several orders of magnitude, so the default is
#' logarithmic axes, on which the equality line is still the diagonal of the
#' square plotting region. Values of zero are drawn at the lower axis limit.
#' Distance from the line reads as a factor; distance along it as size. Read
#' the relevance itself too: a candidate far above the line with negligible
#' relevance is no shadow variable (see `relevance_chisq_df`). Rejected
#' candidates are drawn grey and labelled "rejected", so that a weak candidate
#' above the line is not mistaken for a pass. Hollow points mark candidates
#' whose exclusion violation is at noise level (`ratio_readable` is `FALSE`):
#' no violation is measured, their horizontal position is not estimated, and
#' their label gives no ratio.
#'
#' @param x A `"shadow_screen"` object.
#' @param scale `"log"` (default) or `"linear"` axes.
#' @param label Candidates to label. Default: all if there are at most fifteen,
#'   else the most relevant one and the one with the highest ratio.
#' @param col Point colour.
#' @param text_col Colour for text and reference lines.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Passed to [graphics::plot.default()].
#' @return `x`, invisibly.
#' @export
#' @examples
#' d <- simulate_shadow(5000, tilt = -1, seed = 1)
#' plot(shadow_screen(d, "y", "observed", c("x1", "x2"),
#'                    c("z_shadow", "z_policy", "z_score")))
plot.shadow_screen <- function(x, scale = c("log", "linear"), label = NULL,
                               col = scan_colors[["sharp"]],
                               text_col = "grey20",
                               xlab = expression("Exclusion violation, incremental McFadden " * R^2),
                               ylab = expression("Relevance, incremental McFadden " * R^2),
                               main = NULL, ...) {
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
  } else {
    lim <- c(0, max(re, ex) * 1.08)
  }
  op <- graphics::par(pty = "s")
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(ex, re, type = "n", xlim = lim, ylim = lim,
                 log = if (scale == "log") "xy" else "",
                 xlab = xlab, ylab = ylab, axes = FALSE,
                 main = main, bty = "l", col.lab = text_col,
                 col.main = text_col, fg = text_col, ...)
  if (scale == "log") {
    dec <- seq(log10(lim[1]), log10(lim[2]))
    labs_all <- parse(text = paste0("10^", dec))
    # ticks at every decade; labels only as densely as they fit
    span <- diff(graphics::grconvertX(lim, "user", "inches"))
    lw <- max(graphics::strwidth(labs_all, units = "inches", cex = graphics::par("cex.axis")))
    keep <- dec %in% thin_decades(dec, span / max(1, length(dec) - 1), lw)
    for (s in 1:2) {
      graphics::axis(s, at = 10^dec, labels = FALSE, col = text_col)
      graphics::axis(s, at = 10^dec[keep], labels = labs_all[keep], las = 1,
                     col = text_col, col.axis = text_col, lwd = 0)
    }
  } else {
    for (s in 1:2) graphics::axis(s, las = 1, col = text_col, col.axis = text_col)
  }
  graphics::box(bty = "l", col = text_col)
  usr <- graphics::par("usr")
  # equality line: the diagonal in the coordinates actually drawn
  graphics::segments(lim[1], lim[1], lim[2], lim[2], lty = 2, col = text_col)
  ux <- function(f) if (scale == "log") 10^(usr[1] + f * (usr[2] - usr[1])) else usr[1] + f * (usr[2] - usr[1])
  uy <- function(f) if (scale == "log") 10^(usr[3] + f * (usr[4] - usr[3])) else usr[3] + f * (usr[4] - usr[3])
  ann <- "more relevance than violation"
  graphics::text(ux(0.04), uy(0.96), ann,
                 adj = c(0, 1), col = text_col, cex = 0.8, font = 3)
  ax0 <- graphics::grconvertX(ux(0.04), "user", "inches")
  ay1 <- graphics::grconvertY(uy(0.96), "user", "inches")
  ann_box <- c(ax0, ax0 + graphics::strwidth(ann, units = "inches", cex = 0.8, font = 3),
               ay1 - graphics::strheight(ann, units = "inches", cex = 0.8, font = 3) * 1.3, ay1)
  # a ratio over an exclusion violation below its noise floor is not an
  # estimate: draw such points hollow and do not print their ratio
  readable <- if (is.null(x$ratio_readable)) rep(TRUE, nrow(x)) else x$ratio_readable
  verdict <- if (is.null(x$verdict)) rep("pass", nrow(x)) else x$verdict
  st <- screen_point_style(verdict, readable, col)
  graphics::points(ex, re, pch = 21, cex = 1.6, lwd = 2, bg = st$bg, col = st$border)
  i <- match(label, x$candidate)
  i <- i[!is.na(i)]
  if (length(i)) {
    labs <- ifelse(readable[i],
                   sprintf("ratio %s", format_ratio(x$ratio[i])),
                   "violation at noise level")
    labs <- sprintf("%s (%s%s)", x$candidate[i], labs,
                    ifelse(verdict[i] == "reject", "; rejected", ""))
    place_labels(ex[i], re[i], labs, px = ex, py = re, line = lim,
                 col = text_col, cex = 0.8, obstacles = list(ann_box))
  }
  invisible(x)
}

format_ratio <- function(r) {
  vapply(r, function(v) {
    if (!is.finite(v)) "Inf"
    else if (v >= 100) formatC(round(v), format = "d", big.mark = ",")
    else if (v >= 0.01) formatC(v, digits = 2, format = "f")
    else formatC(v, digits = 1, format = "e")
  }, "")
}

# Greedy label placement: for each point, in order of position, try sixteen
# positions around it (eight directions, near and far) and take the first
# whose text box overlaps no label already placed, no point, no obstacle box
# and not the equality line from line[1] to line[2], and whose leader line
# crosses neither that line nor another label. Labels stay inside the axes at
# left and bottom and inside the device at right and top. Works in inches, so
# it is the same on linear and log axes. If no position is free, the one with
# the least overlap is used.
place_labels <- function(x, y, labels, px, py, line, col, cex, obstacles = list()) {
  ix <- graphics::grconvertX(x, "user", "inches")
  iy <- graphics::grconvertY(y, "user", "inches")
  pix <- graphics::grconvertX(px, "user", "inches")
  piy <- graphics::grconvertY(py, "user", "inches")
  usr <- graphics::par("usr")
  bx <- graphics::grconvertX(c(usr[1], usr[2]), "user", "inches")
  if (graphics::par("xlog")) bx <- graphics::grconvertX(10^usr[1:2], "user", "inches")
  by <- graphics::grconvertY(c(usr[3], usr[4]), "user", "inches")
  if (graphics::par("ylog")) by <- graphics::grconvertY(10^usr[3:4], "user", "inches")
  # the axes are at the left and bottom only: to the right and top a label
  # may run into the empty margin, up to a small gap before the device edge
  bx[2] <- graphics::grconvertX(1, "ndc", "inches") - 0.05
  by[2] <- graphics::grconvertY(1, "ndc", "inches") - 0.05
  w <- graphics::strwidth(labels, units = "inches", cex = cex)
  h <- graphics::strheight("Xg", units = "inches", cex = cex) * 1.25
  r <- 0.09                                    # point radius plus gap, inches
  # candidate positions: eight directions at the point, then the same eight
  # further out; a label not touching its point gets a short leader line
  dirs <- list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1),
               c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
  cand <- c(lapply(dirs, function(d) list(d = d, s = 1)),
            lapply(dirs, function(d) list(d = d, s = 3)))
  box_of <- function(k, cd) {
    d <- cd$d; g <- r * cd$s
    cx <- ix[k] + d[1] * (g + if (d[1] == 0) 0 else w[k] / 2)
    cy <- iy[k] + d[2] * (g + h / 2)
    c(cx - w[k] / 2, cx + w[k] / 2, cy - h / 2, cy + h / 2)
  }
  overlap <- function(a, b)
    max(0, min(a[2], b[2]) - max(a[1], b[1])) * max(0, min(a[4], b[4]) - max(a[3], b[3]))
  pts <- lapply(seq_along(pix), function(j)
    c(pix[j] - r * 0.7, pix[j] + r * 0.7, piy[j] - r * 0.7, piy[j] + r * 0.7))
  # the equality line as a chain of small boxes
  t <- seq(0, 1, length.out = 120)
  lv <- if (graphics::par("xlog")) 10^(log10(line[1]) + t * diff(log10(line))) else line[1] + t * diff(line)
  lx <- graphics::grconvertX(lv, "user", "inches")
  ly <- graphics::grconvertY(lv, "user", "inches")
  e <- 0.015
  pts <- c(pts, lapply(seq_along(lx), function(j) c(lx[j] - e, lx[j] + e, ly[j] - e, ly[j] + e)))
  diag_seg <- c(lx[1], ly[1], lx[length(lx)], ly[length(ly)])
  # leader from the point's edge to the nearest point of the box, or NULL
  leader_of <- function(k, cd, b) {
    if (cd$s == 1 && any(cd$d == 0)) return(NULL)
    tx <- min(max(ix[k], b[1]), b[2]); ty <- min(max(iy[k], b[3]), b[4])
    v <- c(tx - ix[k], ty - iy[k]); len <- sqrt(sum(v^2))
    if (len <= r) return(NULL)
    c(c(ix[k], iy[k]) + v / len * r * 0.6, tx, ty)
  }
  placed <- list(); leaders <- list()
  for (k in order(ix, iy)) {
    cost <- vapply(cand, function(cd) {
      b <- box_of(k, cd)
      out <- max(0, bx[1] - b[1]) + max(0, b[2] - bx[2]) +
             max(0, by[1] - b[3]) + max(0, b[4] - by[2])
      ld <- leader_of(k, cd, b)
      cross <- !is.null(ld) && (segments_cross(ld, diag_seg) ||
        any(vapply(c(placed, obstacles), function(o) runs_through(ld, o), TRUE)))
      cut <- any(vapply(leaders, function(l) runs_through(l, b), TRUE))
      sum(vapply(c(placed, pts, obstacles), function(o) overlap(b, o), 0)) +
        10 * out + 1 * (cross + cut) + 1e-4 * (cd$s - 1)  # near ring on a tie
    }, 0)
    cd <- cand[[which.min(cost)]]             # first minimum = preferred order
    b <- box_of(k, cd)
    placed[[length(placed) + 1L]] <- b
    ld <- leader_of(k, cd, b)
    if (!is.null(ld)) {
      leaders[[length(leaders) + 1L]] <- ld
      graphics::segments(graphics::grconvertX(ld[1], "inches", "user"),
                         graphics::grconvertY(ld[2], "inches", "user"),
                         graphics::grconvertX(ld[3], "inches", "user"),
                         graphics::grconvertY(ld[4], "inches", "user"),
                         col = col, lwd = 0.6, xpd = NA)
    }
    graphics::text(graphics::grconvertX(mean(b[1:2]), "inches", "user"),
                   graphics::grconvertY(mean(b[3:4]), "inches", "user"),
                   labels[k], col = col, cex = cex, xpd = NA)
  }
  # everything in inches, for tests
  invisible(list(boxes = placed, leaders = leaders, obstacles = obstacles,
                 line = diag_seg))
}

# Proper intersection of two segments, each given as c(x0, y0, x1, y1).
segments_cross <- function(a, b) {
  o <- function(p, q, s) sign((q[1] - p[1]) * (s[2] - p[2]) - (q[2] - p[2]) * (s[1] - p[1]))
  o(a[1:2], a[3:4], b[1:2]) * o(a[1:2], a[3:4], b[3:4]) < 0 &&
    o(b[1:2], b[3:4], a[1:2]) * o(b[1:2], b[3:4], a[3:4]) < 0
}

# Whether segment c(x0, y0, x1, y1) passes through the inside of box
# c(left, right, bottom, top); touching an edge does not count.
runs_through <- function(l, b) {
  t <- seq(0, 1, length.out = 50)
  x <- l[1] + t * (l[3] - l[1]); y <- l[2] + t * (l[4] - l[2])
  any(x > b[1] & x < b[2] & y > b[3] & y < b[4])
}

# Which decades of a log axis get a label: every step-th, counted from the
# top decade, with step the smallest that leaves 20% room between labels.
thin_decades <- function(decades, inch_per_decade, label_width) {
  step <- max(1L, ceiling(1.2 * label_width / inch_per_decade))
  decades[(max(decades) - decades) %% step == 0]
}

# Point style of the screen plot: colour carries the verdict (rejected grey),
# fill carries readability (hollow where the violation is at noise level).
screen_point_style <- function(verdict, readable, col, reject_col = "grey60") {
  hue <- ifelse(verdict == "reject", reject_col, col)
  list(bg = ifelse(readable, hue, surface_col()),
       border = ifelse(readable, surface_col(), hue))
}

surface_col <- function() {
  bg <- graphics::par("bg")
  if (identical(bg, "transparent")) "white" else bg
}
