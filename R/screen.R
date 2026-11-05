#' Screen candidate shadow variables
#'
#' `shadow_screen()` asks whether a variable \eqn{Z} is fit to serve as a
#' **shadow variable** for a binary outcome \eqn{Y} that is missing not at
#' random: relevant for \eqn{Y} given \eqn{X}, and excluded from the response
#' model given \eqn{X} and \eqn{Y}. It is a different tool from
#' [shadow_scan()]: the screen presupposes a shadow-variable model that
#' point-identifies the selection parameter; the scan needs no shadow variable
#' and returns an identified set. Do not read the two as one procedure.
#'
#' For each candidate four nested `glm`s are fitted: \eqn{Y} on \eqn{X} and on
#' \eqn{X + Z} **among the observed rows**, and \eqn{R} on \eqn{X} and on
#' \eqn{X + Z} on the **pooled** sample. Two axes are reported, both as
#' incremental McFadden \eqn{R^2}, so that they share one scale:
#' \deqn{\Delta R^2 = \frac{D(X) - D(X + Z)}{D_0},}
#' the drop in deviance from adding \eqn{Z}, divided by the null deviance of
#' the same response on the same rows.
#' * **relevance**: \eqn{\Delta R^2} of \eqn{Z} in the outcome model;
#' * **exclusion**: \eqn{\Delta R^2} of \eqn{Z} in the response model, i.e.
#'   the size of the exclusion violation;
#' * **ratio**: relevance / exclusion. Above one, \eqn{Z} predicts the outcome
#'   better than it predicts selection.
#'
#' As a second column for each axis the likelihood-ratio \eqn{\chi^2} per
#' degree of freedom is reported, the form in which relevance is read against
#' the conventional first-stage standard of 10. It is the measured
#' \eqn{\Delta R^2} divided by its noise floor: without any effect,
#' \eqn{\chi^2} has mean df, so the same coding yields \eqn{\Delta R^2}
#' of about df / \eqn{D_0} from noise alone. `min_relevance = 10` asks for
#' ten times that floor. The floor grows with the degrees of freedom, so a
#' finer coding raises \eqn{\Delta R^2} by itself; \eqn{\chi^2}/df does not.
#' The same holds for the denominator of the ratio: where the exclusion
#' \eqn{\chi^2}/df is below 1, the violation is not told apart from noise and
#' the ratio is a division by noise. Such ratios are flagged in
#' `ratio_readable`, by `print()` and by hollow points in the plot. In large
#' samples every violation is significant, so exclusion should be read as an
#' effect size, not a test.
#'
#' Exclusion is measured *marginally* over \eqn{Y}: the pooled response model
#' cannot condition on an outcome that is missing. Under nonignorable selection
#' a relevant \eqn{Z} therefore shows some exclusion signal through \eqn{Y}
#' even when it is validly excluded; the screen is a screen, not a test of the
#' restriction.
#'
#' **Numeric characteristics.** With `numeric = "linear"` (the default) a
#' numeric column enters linearly. Where the selection decision uses a
#' characteristic through a threshold -- a lender that declines every
#' application below a score cut-off -- a linear term cannot represent the
#' step, and the exclusion violation is understated, often by an order of
#' magnitude. The screen warns when a numeric characteristic takes values
#' never seen among the observed rows on more than 1% of the sample, the
#' visible trace of such a rule.
#'
#' `numeric = "classes"` codes every numeric covariate and candidate with
#' more than two distinct values by [screen_classes()]: about `k` classes of
#' equal mass, cut only between distinct values, plus a class `"below"` or
#' `"above"` for values outside the range seen among the observed rows. The
#' relevance hardly changes; the exclusion axis does. The cost is degrees of
#' freedom: a candidate with ten classes has nine, so its relevance
#' \eqn{\chi^2} per degree of freedom is lower than with one linear term,
#' and a candidate that is relevant only through a smooth trend can fall
#' below `min_relevance`. Read both codings when they disagree. The coding
#' used is stored in `attr(x, "coding")`. Flexible coding pays only when the
#' sample carries the degrees of freedom: the screen warns when a model has
#' fewer than 10 cases of the rarer class per parameter, separately for the
#' outcome model, fitted to the observed rows, and the response model.
#'
#' Two different questions sit behind the coding. Whether a characteristic
#' needs flexible coding is answered by the gap in support: the share of rows
#' outside the range seen among the observed rows (`outside` in
#' `attr(x, "coding")`), which a threshold in the selection leaves and a
#' threshold that only defines the population does not. Whether classes are
#' well defined is answered by the mass on single values (`heaping`, reported
#' for `numeric = "classes"` and by [screen_classes()]). Classes
#' are cut only between distinct values, so heavy atoms give classes of
#' unequal mass, not arbitrary ones. A loan amount in round numbers shows
#' both: no gap in support, so a linear term suffices; where an analysis
#' codes it in classes all the same, for instance to compare two codings
#' that differ in nothing else, the heaping shows whether that result is
#' well defined. For other cut-offs or class counts
#' per characteristic, code the columns yourself with [screen_classes()] and
#' pass factors.
#'
#' @param data A data frame.
#' @param outcome Name of the outcome column. May be `NA` where unobserved.
#' @param response Name of the response indicator column (0/1 or logical;
#'   1 = outcome observed).
#' @param covariates Character vector of covariate names \eqn{X}; may be empty.
#' @param candidates Character vector of candidate names \eqn{Z}. Factors are
#'   allowed and enter with their degrees of freedom.
#' @param family Family of the outcome model, passed to `glm()`. The response
#'   model is always binomial. The McFadden measure presumes a likelihood
#'   model for discrete outcomes; for a continuous outcome another effect-size
#'   measure would be needed.
#' @param min_relevance Minimum relevance \eqn{\chi^2} per degree of freedom
#'   for a candidate to pass. Default 10.
#' @param min_ratio Minimum ratio of relevance to exclusion for a candidate to
#'   pass. Default 1, the equality line.
#' @param numeric `"linear"` (default) or `"classes"`; see the section on
#'   numeric characteristics.
#' @param k Target number of classes for `numeric = "classes"`. Default 10.
#'
#' @return A data frame of class `"shadow_screen"`, one row per candidate,
#'   with `candidate`, `relevance`, `exclusion`, `ratio`, `relevance_chisq_df`,
#'   `exclusion_chisq_df`, `df`, `n_observed`, `n`, `epp_relevance` and
#'   `epp_exclusion` (cases of the rarer class per parameter of the outcome
#'   and of the response model; `NA` for a target that is not binary),
#'   `ratio_readable` (`FALSE` where the exclusion \eqn{\chi^2}/df is below 1),
#'   `verdict` (`"pass"` or
#'   `"reject"`) and `reason`. Attribute `"coding"`: one row per numeric
#'   characteristic with its treatment, the number of classes, the share of
#'   rows outside the range seen among the observed rows, and, for classes,
#'   the heaping share of [screen_classes()].
#' @seealso [plot.shadow_screen()], [screen_classes()]
#' @export
#' @examples
#' d <- simulate_shadow(5000, tilt = -1, seed = 1)
#' sc <- shadow_screen(d, outcome = "y", response = "observed",
#'                     covariates = c("x1", "x2"),
#'                     candidates = c("z_shadow", "z_policy", "z_score"))
#' sc
#' plot(sc)
#' # numeric characteristics as classes, robust to threshold rules
#' shadow_screen(d, "y", "observed", c("x1", "x2"), "z_score",
#'               numeric = "classes")
shadow_screen <- function(data, outcome, response, covariates = character(),
                          candidates, family = stats::binomial(),
                          min_relevance = 10, min_ratio = 1,
                          numeric = c("linear", "classes"), k = 10) {
  numeric <- match.arg(numeric)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (missing(candidates) || length(candidates) == 0L)
    stop("Give at least one candidate.", call. = FALSE)
  needed <- c(outcome, response, covariates, candidates)
  miss <- setdiff(needed, names(data))
  if (length(miss))
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "),
         call. = FALSE)
  if (length(intersect(candidates, covariates)))
    stop("A candidate cannot also be a covariate.", call. = FALSE)
  rv <- data[[response]]
  if (is.logical(rv)) rv <- as.integer(rv)
  if (!all(rv[!is.na(rv)] %in% c(0, 1)))
    stop("`response` must be 0/1 or logical.", call. = FALSE)
  data[[response]] <- rv
  if (any(rv == 1 & is.na(data[[outcome]]), na.rm = TRUE))
    stop("The outcome is missing for some rows marked as observed.",
         call. = FALSE)

  coded <- code_numeric(data, response, c(covariates, candidates), numeric, k)
  data <- coded$data

  rows <- lapply(candidates, function(z) {
    cols <- c(response, covariates, z)
    # keep only the columns needed, so that subsetting a wide frame stays cheap
    pooled <- data[stats::complete.cases(data[cols]), c(outcome, cols),
                   drop = FALSE]
    obs <- pooled[pooled[[response]] == 1, , drop = FALSE]
    rel <- incr_fit(obs, outcome, covariates, z, family)
    exc <- incr_fit(pooled, response, covariates, z, stats::binomial())
    data.frame(
      candidate = z,
      relevance = rel$r2, exclusion = exc$r2,
      ratio = if (exc$r2 > 0) rel$r2 / exc$r2 else Inf,
      relevance_chisq_df = rel$chisq / rel$df,
      exclusion_chisq_df = exc$chisq / exc$df,
      df = rel$df, n_observed = nrow(obs), n = nrow(pooled),
      ratio_readable = exc$chisq / exc$df >= 1,
      epp_relevance = rel$epp, epp_exclusion = exc$epp,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  thin <- function(v) out$candidate[!is.na(v) & v < 10]
  for (side in list(list(v = out$epp_relevance, m = "outcome model (observed rows)"),
                    list(v = out$epp_exclusion, m = "response model (all rows)")))
    if (length(thin(side$v)))
      warning("Fewer than 10 cases of the rarer class per parameter in the ",
              side$m, " for: ", paste(sprintf("%s (%.1f)", thin(side$v),
                                              side$v[!is.na(side$v) & side$v < 10]),
                                      collapse = ", "),
              ". The incremental McFadden R^2 has no penalty for degrees of ",
              "freedom, so extra parameters can manufacture relevance or ",
              "violation. Prefer linear coding or a smaller `k`.", call. = FALSE)
  weak <- out$relevance_chisq_df < min_relevance
  viol <- out$ratio <= min_ratio
  out$verdict <- ifelse(weak | viol, "reject", "pass")
  out$reason <- ifelse(weak & viol, "weak and violates exclusion",
                ifelse(weak, "weak", ifelse(viol, "violates exclusion", "")))
  attr(out, "settings") <- list(min_relevance = min_relevance,
                                min_ratio = min_ratio, numeric = numeric, k = k)
  attr(out, "coding") <- coded$coding
  class(out) <- c("shadow_screen", "data.frame")
  out
}

# Numeric characteristics with more than two distinct values: report the share
# of rows outside the range seen among the observed rows (the trace of a
# threshold rule) and, for numeric = "classes", replace them by screen_classes().
code_numeric <- function(data, response, vars, numeric, k) {
  num <- vars[vapply(vars, function(v) is.numeric(data[[v]]) &&
                       length(unique(stats::na.omit(data[[v]]))) > 2L, logical(1))]
  if (!length(num)) return(list(data = data, coding = NULL))
  obs <- data[[response]] == 1 & !is.na(data[[response]])
  coding <- do.call(rbind, lapply(num, function(v) {
    x <- data[[v]]
    sup <- range(x[obs], na.rm = TRUE)
    outside <- mean(x < sup[1] | x > sup[2], na.rm = TRUE)
    if (numeric == "classes") {
      f <- screen_classes(x, k = k, support = sup)
      data[[v]] <<- f
      data.frame(variable = v, coding = "classes", classes = nlevels(f),
                 outside = outside, heaping = attr(f, "heaping"),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(variable = v, coding = "linear", classes = NA_integer_,
                 outside = outside, heaping = NA_real_, stringsAsFactors = FALSE)
    }
  }))
  if (numeric == "linear") {
    flag <- coding$variable[coding$outside > 0.01]
    if (length(flag))
      warning("Numeric characteristic(s) with values never seen among the ",
              "observed rows on more than 1% of the sample: ",
              paste(sprintf("%s (%.0f%%)", flag,
                            100 * coding$outside[coding$outside > 0.01]),
                    collapse = ", "),
              ". A threshold rule in the selection is likely; a linear term ",
              "understates its exclusion violation. Consider numeric = \"classes\".",
              call. = FALSE)
  }
  list(data = data, coding = coding)
}

# Incremental McFadden R^2 of z. Only deviances are kept: the design matrix is
# built once, each fit is reduced to its deviance and rank at once, and no glm
# object (model frame, QR, fitted values) is held while the next one is fitted.
# With millions of rows that is the difference between one and several copies
# of an n x p matrix.
incr_fit <- function(d, lhs, covariates, z, family) {
  f1 <- stats::reformulate(c(covariates, z), lhs)
  mf <- stats::model.frame(f1, d)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(f1, mf)
  z_term <- match(z, attr(stats::terms(f1), "term.labels"))
  in_z <- attr(X, "assign") == z_term
  rm(mf)
  fit_dev <- function(X) {
    fit <- suppressWarnings(stats::glm.fit(X, y, family = family))
    c(deviance = fit$deviance, null = fit$null.deviance, rank = fit$rank)
  }
  full <- fit_dev(X)
  X <- X[, !in_z, drop = FALSE]
  red <- fit_dev(X)
  rm(X)
  # events per parameter: cases of the rarer class of a binary target per
  # coefficient of the full model; NA for targets that are not binary
  yy <- if (is.factor(y)) as.integer(y) - 1L else y
  epp <- if (family$family == "binomial" && all(yy %in% c(0, 1)))
    min(sum(yy == 1), sum(yy == 0)) / full[["rank"]] else NA_real_
  chisq <- max(red[["deviance"]] - full[["deviance"]], 0)
  df <- full[["rank"]] - red[["rank"]]
  if (df < 1)
    stop("Candidate `", z, "` adds no degrees of freedom (constant or ",
         "collinear with the covariates).", call. = FALSE)
  list(r2 = chisq / red[["null"]], chisq = chisq, df = df, epp = epp)
}

#' @export
print.shadow_screen <- function(x, digits = 3, ...) {
  s <- attr(x, "settings")
  cat("<shadow_screen> relevance vs exclusion violation, incremental McFadden R^2\n")
  cat(sprintf("  pass: relevance chi2/df >= %s and ratio > %s\n",
              format(s$min_relevance), format(s$min_ratio)))
  cd <- attr(x, "coding")
  if (!is.null(cd) && nrow(cd)) {
    if (all(cd$coding == "classes")) {
      cat(sprintf("  numeric as classes: %s\n", paste(sprintf("%s (%d)",
          cd$variable, cd$classes), collapse = ", ")))
    } else {
      cat(sprintf("  numeric as linear terms: %s\n",
                  paste(cd$variable, collapse = ", ")))
    }
  }
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(df, "settings") <- NULL
  show <- df[c("candidate", "relevance", "exclusion", "ratio",
               "relevance_chisq_df", "verdict", "reason")]
  num <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], signif, digits = digits)
  print(show, row.names = FALSE)
  if (!is.null(x$ratio_readable) && any(!x$ratio_readable))
    cat(sprintf(paste0("  ratio not readable (exclusion chi2/df below 1, so the ",
                       "violation is not told apart from noise): %s\n"),
                paste(x$candidate[!x$ratio_readable], collapse = ", ")))
  invisible(x)
}
