scan_colors <- c(sharp = "#2a78d6", naive = "#eb6834")

#' Sensitivity curve of a scan
#'
#' Plots the portfolio rate against selection strength: the sharp identified
#' set as a filled band with solid edges, the naive band
#' \eqn{\mathrm{expit}(\ell \pm \lambda)} dashed and lighter, the point that a
#' model on the observed rows asserts (\eqn{\lambda = 0}) marked, and, if the
#' scan has a `threshold`, the decision threshold with a drop line where the
#' sharp set first reaches it. That crossing is the reason to run the scan: it
#' shows at which strength of selection the decision flips.
#'
#' Where \eqn{\hat r} is small, as for rejected applications, the upper edge of
#' the sharp set coincides with the naive one; see [shadow_scan()].
#'
#' @param x A `"shadow_scan"` object.
#' @param scale `"odds"` (default) puts \eqn{e^\lambda}, the largest admissible
#'   ratio of observation propensities, on the horizontal axis, starting at 1;
#'   `"logodds"` puts \eqn{\lambda}.
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
plot.shadow_scan <- function(x, scale = c("odds", "logodds"),
                             colors = scan_colors, text_col = "grey20",
                             xlab = NULL, ylab = "Portfolio rate", main = NULL,
                             ...) {
  scale <- match.arg(scale)
  cv <- x$curve
  tx <- function(l) if (scale == "odds") exp(l) else l
  xs <- tx(cv$lambda)
  if (is.null(xlab)) xlab <- if (scale == "odds")
    expression("Selection strength " * e^lambda * " (ratio of observation propensities)")
  else expression("Selection strength " * lambda * " (log scale)")
  ylim <- range(cv$naive_lower, cv$naive_upper, x$threshold)
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

  if (min(cv$lambda) == 0) {
    graphics::points(tx(0), x$point, pch = 21, cex = 1.4, lwd = 2,
                     bg = colors[["sharp"]], col = surface_col())
    graphics::text(tx(0), x$point, "model on observed rows",
                   adj = c(-0.08, 2.2), col = text_col, cex = 0.8)
  }
  if (!is.null(x$threshold)) {
    usr <- graphics::par("usr")
    graphics::abline(h = x$threshold, lty = 2, col = text_col)
    graphics::text(usr[2], x$threshold, "decision threshold", adj = c(1, 1.6),
                   col = text_col, cex = 0.8)
    fs <- x$flip[["sharp"]]; fn <- x$flip[["naive"]]
    sym <- if (scale == "odds") "e^lambda" else "lambda"
    lab <- function(v) bquote(.(as.name("flips at")) ~ .(parse(text = sym)[[1]]) %~~% .(sprintf("%.1f", tx(v))))
    if (!is.na(fs)) {
      graphics::segments(tx(fs), usr[3], tx(fs), x$threshold, lty = 3,
                         col = text_col)
      graphics::points(tx(fs), x$threshold, pch = 21, cex = 1.4, lwd = 2,
                       bg = colors[["sharp"]], col = surface_col())
      graphics::text(tx(fs), usr[3], lab(fs), adj = c(-0.08, -0.6),
                     col = text_col, cex = 0.85)
    } else if (!is.na(fn)) {
      graphics::segments(tx(fn), usr[3], tx(fn), x$threshold, lty = 3,
                         col = text_col)
      graphics::points(tx(fn), x$threshold, pch = 21, cex = 1.1, lwd = 2,
                       bg = colors[["naive"]], col = surface_col())
      graphics::text(tx(fn), usr[3],
                     sprintf("naive band flips at %.1f; sharp set never",
                             tx(fn)),
                     adj = c(-0.05, -0.6), col = text_col, cex = 0.8)
    }
  }
  graphics::legend("topleft", bty = "n", text.col = text_col, cex = 0.85,
                   legend = c(
                     bquote("Sharp identified set (" * .(x$aggregate) * " tilt)"),
                     quote("Naive band " * expit(logit ~ eta[obs] %+-% lambda))
                   ),
                   fill = c(grDevices::adjustcolor(colors[["sharp"]], 0.18), NA),
                   border = c(colors[["sharp"]], NA),
                   lty = c(NA, 2), lwd = c(NA, 1.5),
                   col = c(NA, naive_col), merge = TRUE, x.intersp = 0.8)
  invisible(x)
}

#' Relevance against exclusion violation
#'
#' Scatter of the screen: exclusion violation on the horizontal axis,
#' relevance on the vertical, both incremental McFadden \eqn{R^2} on the same
#' scale, with the equality line dashed from the origin. Above the line a
#' candidate predicts the outcome better than it predicts selection. Points are
#' labelled directly rather than through a legend.
#'
#' @param x A `"shadow_screen"` object.
#' @param label Candidates to label. Default: all if there are at most eight,
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
plot.shadow_screen <- function(x, label = NULL, col = scan_colors[["sharp"]],
                               text_col = "grey20",
                               xlab = expression("Exclusion violation, incremental McFadden " * R^2),
                               ylab = expression("Relevance, incremental McFadden " * R^2),
                               main = NULL, ...) {
  if (is.null(label)) {
    label <- if (nrow(x) <= 8L) x$candidate else
      unique(c(x$candidate[which.max(x$relevance)],
               x$candidate[which.max(x$ratio)]))
  }
  lim <- c(0, max(x$relevance, x$exclusion) * 1.08)
  op <- graphics::par(pty = "s")
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(x$exclusion, x$relevance, type = "n", xlim = lim, ylim = lim,
                 xlab = xlab, ylab = ylab,
                 main = main, las = 1, bty = "l", col.axis = text_col,
                 col.lab = text_col, col.main = text_col, fg = text_col, ...)
  top <- lim[2]
  graphics::segments(0, 0, top, top, lty = 2, col = text_col)
  graphics::text(top * 0.04, top * 0.96, "more relevance than violation",
                 adj = c(0, 1), col = text_col, cex = 0.8, font = 3)
  graphics::points(x$exclusion, x$relevance, pch = 21, cex = 1.6, lwd = 2,
                   bg = col, col = surface_col())
  i <- match(label, x$candidate)
  i <- i[!is.na(i)]
  if (length(i)) {
    labs <- sprintf("%s (ratio %s)", x$candidate[i],
                    formatC(x$ratio[i], digits = 2, format = "f"))
    right <- x$exclusion[i] > 0.55 * lim[2]
    graphics::text(x$exclusion[i], x$relevance[i], labs,
                   pos = ifelse(right, 2, 4), col = text_col, cex = 0.8,
                   offset = 0.7, xpd = NA)
  }
  invisible(x)
}

surface_col <- function() {
  bg <- graphics::par("bg")
  if (identical(bg, "transparent")) "white" else bg
}
