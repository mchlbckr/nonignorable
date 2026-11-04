scan_colors <- c(sharp = "#2a78d6", naive = "#eb6834")

#' Sensitivity curve of a scan
#'
#' Plots the portfolio rate (or mean, for [shadow_scan_mean()]) against
#' selection strength: the sharp identified set as a filled band with solid
#' edges, the naive band \eqn{\mathrm{expit}(\ell \pm \lambda)} (for a mean:
#' the outer band without the cap \eqn{\pi \le 1}) dashed and lighter, the point that a
#' model on the observed rows asserts (\eqn{\lambda = 0}) marked, and, if the
#' scan has a `threshold`, the decision threshold with a drop line where the
#' sharp set first reaches it. That crossing is the reason to run the scan: it
#' shows at which strength of selection the decision flips.
#'
#' Where \eqn{\hat r} is small, as for rejected applications, the upper edge of
#' the sharp set coincides with the naive one; see [shadow_scan()].
#'
#' @param x A `"shadow_scan"` or `"shadow_scan_mean"` object.
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
#' plot(shadow_scan(rej$eta_obs, rej$obs_rate, threshold = 0.35))
plot.shadow_scan <- function(x, scale = c("ratio", "log"),
                             colors = scan_colors, text_col = "grey20",
                             xlab = NULL, ylab = NULL, main = NULL,
                             ...) {
  scale <- match.arg(scale)
  cv <- x$curve
  is_mean <- inherits(x, "shadow_scan_mean")
  if (is.null(ylab)) ylab <- if (is_mean) "Portfolio mean" else "Portfolio rate"
  tx <- function(l) if (scale == "ratio") exp(l) else l
  xs <- tx(cv$lambda)
  if (is.null(xlab)) xlab <- if (scale == "ratio")
    expression("Selection strength " * e^lambda * " (ratio of response propensities)")
  else expression("Selection strength " * lambda * " (log of the ratio)")
  ylim <- range(cv$naive_lower, cv$naive_upper, cv$ci_lower, cv$ci_upper,
                x$threshold)
  ylim <- ylim + c(-0.04, 0.04) * diff(ylim)

  graphics::plot(xs, cv$upper, type = "n", ylim = ylim, xlab = xlab,
                 ylab = ylab, main = main, las = 1, bty = "l",
                 col.axis = text_col, col.lab = text_col, col.main = text_col,
                 fg = text_col, ...)
  graphics::polygon(c(xs, rev(xs)), c(cv$lower, rev(cv$upper)),
                    col = grDevices::adjustcolor(colors[["sharp"]], 0.18),
                    border = NA)
  naive_col <- grDevices::adjustcolor(colors[["naive"]], 0.75)
  graphics::lines(xs, cv$naive_upper, lty = 2, lwd = 1.5, col = naive_col)
  graphics::lines(xs, cv$naive_lower, lty = 2, lwd = 1.5, col = naive_col)
  graphics::lines(xs, cv$upper, lwd = 2, col = colors[["sharp"]])
  graphics::lines(xs, cv$lower, lwd = 2, col = colors[["sharp"]])
  has_ci <- !is.null(cv$ci_lower)
  if (has_ci) {
    graphics::lines(xs, cv$ci_upper, lwd = 1.3, lty = 3, col = colors[["sharp"]])
    graphics::lines(xs, cv$ci_lower, lwd = 1.3, lty = 3, col = colors[["sharp"]])
  }

  # The lambda = 0 point is explained in the legend, not by a label next to it:
  # when the decision flips near lambda = 0, the flip marker, its drop line and
  # the threshold all sit beside that point and any label there collides.
  has_point <- min(cv$lambda) == 0
  if (has_point)
    graphics::points(tx(0), x$point, pch = 21, cex = 1.4, lwd = 2,
                     bg = colors[["sharp"]], col = surface_col())
  if (!is.null(x$threshold)) {
    usr <- graphics::par("usr")
    graphics::abline(h = x$threshold, lty = 2, col = text_col)
    graphics::text(usr[2], x$threshold, "decision threshold", adj = c(1, 1.6),
                   col = text_col, cex = 0.8)
    fs <- x$flip[["sharp"]]; fn <- x$flip[["naive"]]
    sym <- if (scale == "ratio") "e^lambda" else "lambda"
    lab <- function(v) bquote(.(as.name("flips at")) ~ .(parse(text = sym)[[1]]) %~~% .(sprintf("%.1f", tx(v))))
    # label to the right of the drop line, or to its left near the right edge
    flip_adj <- function(v, pad) {
      if (tx(v) > usr[1] + 0.7 * (usr[2] - usr[1])) c(1 + pad, -0.6) else c(-pad, -0.6)
    }
    if (!is.na(fs)) {
      graphics::segments(tx(fs), usr[3], tx(fs), x$threshold, lty = 3,
                         col = text_col)
      graphics::points(tx(fs), x$threshold, pch = 21, cex = 1.3, lwd = 2,
                       bg = surface_col(), col = colors[["sharp"]])
      graphics::text(tx(fs), usr[3], lab(fs), adj = flip_adj(fs, 0.08),
                     col = text_col, cex = 0.85)
    } else if (!is.na(fn)) {
      graphics::segments(tx(fn), usr[3], tx(fn), x$threshold, lty = 3,
                         col = text_col)
      graphics::points(tx(fn), x$threshold, pch = 21, cex = 1.1, lwd = 2,
                       bg = surface_col(), col = colors[["naive"]])
      graphics::text(tx(fn), usr[3],
                     sprintf("%s band flips at %.1f; sharp set never",
                             if (is_mean) "outer" else "naive", tx(fn)),
                     adj = flip_adj(fn, 0.05), col = text_col, cex = 0.8)
    }
    fc <- if (inherits(x, "shadow_scan_ci")) x$flip[["ci"]] else NA_real_
    if (!is.na(fc)) {
      # where the confidence set reaches the threshold: from there on the data
      # no longer rule out the other decision
      graphics::segments(tx(fc), usr[3], tx(fc), x$threshold, lty = 3,
                         col = text_col)
      graphics::points(tx(fc), x$threshold, pch = 21, cex = 0.9, lwd = 1.5,
                       bg = surface_col(), col = colors[["sharp"]])
      graphics::text(tx(fc), usr[3],
                     bquote("confidence set at" ~ .(parse(text = sym)[[1]]) %~~%
                              .(sprintf("%.1f", tx(fc)))),
                     adj = c(flip_adj(fc, 0.08)[1], -2.2), col = text_col, cex = 0.8)
    }
  }
  items <- list(
    list(lab = if (inherits(x, "shadow_scan_ci")) quote("Sharp identified set, estimated")
               else if (!is_mean && x$aggregate == "common") quote("Sharp identified set (common tilt)")
               else quote("Sharp identified set"),
         fill = grDevices::adjustcolor(colors[["sharp"]], 0.18), border = colors[["sharp"]],
         lty = NA, lwd = NA, pch = NA, pt.bg = NA, col = NA),
    list(lab = if (is_mean) quote("Outer band (without " * pi <= 1 * ")")
               else quote("Naive band " * expit(logit ~ eta[obs] %+-% lambda)),
         fill = NA, border = NA, lty = 2, lwd = 1.5, pch = NA, pt.bg = NA, col = naive_col),
    if (has_ci) list(lab = bquote(.(sprintf("%g%% confidence set", 100 * x$level))),
                     fill = NA, border = NA, lty = 3, lwd = 1.3, pch = NA, pt.bg = NA,
                     col = colors[["sharp"]]),
    if (has_point) list(lab = quote("Model on observed rows " * (lambda == 0)),
                        fill = NA, border = NA, lty = NA, lwd = 2, pch = 21,
                        pt.bg = colors[["sharp"]], col = surface_col())
  )
  items <- Filter(Negate(is.null), items)
  gc_ <- function(nm) unlist(lapply(items, function(it) if (is.na(it[[nm]])) NA else it[[nm]]))
  graphics::legend("topleft", bty = "n", text.col = text_col, cex = 0.85,
                   legend = as.expression(lapply(items, `[[`, "lab")),
                   fill = gc_("fill"), border = gc_("border"),
                   lty = gc_("lty"), lwd = gc_("lwd"), pch = gc_("pch"),
                   pt.cex = 1.2, pt.bg = gc_("pt.bg"), col = gc_("col"),
                   merge = TRUE, x.intersp = 0.8)
  invisible(x)
}

#' Relevance against exclusion violation
#'
#' Scatter of the screen: exclusion violation on the horizontal axis,
#' relevance on the vertical, both incremental McFadden \eqn{R^2} on the same
#' scale, with the equality line dashed. Above the line a candidate predicts
#' the outcome better than it predicts selection. Points are labelled directly
#' rather than through a legend; labels are placed so that they do not overlap
#' each other or the points.
#'
#' Both measures typically span several orders of magnitude, so the default is
#' logarithmic axes, on which the equality line is still the diagonal of the
#' square plotting region. Values of zero are drawn at the lower axis limit.
#' Distance from the line reads as a factor; distance along it as size. Read
#' the relevance itself too: a candidate far above the line with negligible
#' relevance is no shadow variable (see `relevance_chisq_df`).
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
    at <- 10^seq(log10(lim[1]), log10(lim[2]))
    labs_at <- parse(text = paste0("10^", log10(at)))
    for (s in 1:2) graphics::axis(s, at = at, labels = labs_at, las = 1,
                                  col = text_col, col.axis = text_col)
  } else {
    for (s in 1:2) graphics::axis(s, las = 1, col = text_col, col.axis = text_col)
  }
  graphics::box(bty = "l", col = text_col)
  usr <- graphics::par("usr")
  # equality line: the diagonal in the coordinates actually drawn
  graphics::segments(lim[1], lim[1], lim[2], lim[2], lty = 2, col = text_col)
  ux <- function(f) if (scale == "log") 10^(usr[1] + f * (usr[2] - usr[1])) else usr[1] + f * (usr[2] - usr[1])
  uy <- function(f) if (scale == "log") 10^(usr[3] + f * (usr[4] - usr[3])) else usr[3] + f * (usr[4] - usr[3])
  graphics::text(ux(0.04), uy(0.96), "more relevance than violation",
                 adj = c(0, 1), col = text_col, cex = 0.8, font = 3)
  graphics::points(ex, re, pch = 21, cex = 1.6, lwd = 2,
                   bg = col, col = surface_col())
  i <- match(label, x$candidate)
  i <- i[!is.na(i)]
  if (length(i)) {
    labs <- sprintf("%s (ratio %s)", x$candidate[i], format_ratio(x$ratio[i]))
    place_labels(ex[i], re[i], labs, px = ex, py = re, line = lim,
                 col = text_col, cex = 0.8)
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

# Greedy label placement: for each point, in order of position, try eight
# positions around it and take the first whose text box overlaps neither a
# label already placed, nor any point, nor the equality line from line[1] to
# line[2], nor the plot border. Works in inches, so it is the same on linear
# and log axes. If no position is free, the one with the least overlap is used.
place_labels <- function(x, y, labels, px, py, line, col, cex) {
  ix <- graphics::grconvertX(x, "user", "inches")
  iy <- graphics::grconvertY(y, "user", "inches")
  pix <- graphics::grconvertX(px, "user", "inches")
  piy <- graphics::grconvertY(py, "user", "inches")
  usr <- graphics::par("usr")
  bx <- graphics::grconvertX(c(usr[1], usr[2]), "user", "inches")
  if (graphics::par("xlog")) bx <- graphics::grconvertX(10^usr[1:2], "user", "inches")
  by <- graphics::grconvertY(c(usr[3], usr[4]), "user", "inches")
  if (graphics::par("ylog")) by <- graphics::grconvertY(10^usr[3:4], "user", "inches")
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
  placed <- list()
  for (k in order(ix, iy)) {
    cost <- vapply(cand, function(cd) {
      b <- box_of(k, cd)
      out <- max(0, bx[1] - b[1]) + max(0, b[2] - bx[2]) +
             max(0, by[1] - b[3]) + max(0, b[4] - by[2])
      sum(vapply(c(placed, pts), function(o) overlap(b, o), 0)) +
        10 * out + 1e-4 * (cd$s - 1)          # prefer the near ring on a tie
    }, 0)
    cd <- cand[[which.min(cost)]]             # first minimum = preferred order
    b <- box_of(k, cd)
    placed[[length(placed) + 1L]] <- b
    if (cd$s > 1 || all(cd$d != 0)) {
      # leader from the point's edge to the nearest point of the box
      tx <- min(max(ix[k], b[1]), b[2]); ty <- min(max(iy[k], b[3]), b[4])
      v <- c(tx - ix[k], ty - iy[k]); len <- sqrt(sum(v^2))
      if (len > r) {
        p0 <- c(ix[k], iy[k]) + v / len * r * 0.6
        graphics::segments(graphics::grconvertX(p0[1], "inches", "user"),
                           graphics::grconvertY(p0[2], "inches", "user"),
                           graphics::grconvertX(tx, "inches", "user"),
                           graphics::grconvertY(ty, "inches", "user"),
                           col = col, lwd = 0.6, xpd = NA)
      }
    }
    graphics::text(graphics::grconvertX(mean(b[1:2]), "inches", "user"),
                   graphics::grconvertY(mean(b[3:4]), "inches", "user"),
                   labels[k], col = col, cex = cex, xpd = NA)
  }
  invisible(NULL)
}

surface_col <- function() {
  bg <- graphics::par("bg")
  if (identical(bg, "transparent")) "white" else bg
}
