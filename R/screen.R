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
#' the conventional first-stage standard of 10. In large samples every
#' violation is significant, so exclusion should be read as an effect size, not
#' a test.
#'
#' Exclusion is measured *marginally* over \eqn{Y}: the pooled response model
#' cannot condition on an outcome that is missing. Under nonignorable selection
#' a relevant \eqn{Z} therefore shows some exclusion signal through \eqn{Y}
#' even when it is validly excluded; the screen is a screen, not a test of the
#' restriction.
#'
#' **Code numeric characteristics flexibly.** Every column enters the models
#' as it is given: a numeric column linearly, a factor with its levels. Where
#' the selection decision uses a characteristic through a threshold -- a
#' lender that declines every application below a score cut-off -- a linear
#' term cannot represent the step, and the exclusion violation of that
#' characteristic is understated, often by an order of magnitude. Pass such
#' characteristics as factors (for instance quantile classes plus the known
#' cut-off) or as spline terms built beforehand, both as candidates and as
#' covariates. The relevance hardly changes; the exclusion axis does. Check
#' that the result does not depend on how ties at class boundaries are
#' assigned.
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
#'
#' @return A data frame of class `"shadow_screen"`, one row per candidate,
#'   with `candidate`, `relevance`, `exclusion`, `ratio`, `relevance_chisq_df`,
#'   `exclusion_chisq_df`, `df`, `n_observed`, `n`, `verdict` (`"pass"` or
#'   `"reject"`) and `reason`.
#' @seealso [plot.shadow_screen()]
#' @export
#' @examples
#' d <- simulate_shadow(5000, tilt = -1, seed = 1)
#' sc <- shadow_screen(d, outcome = "y", response = "observed",
#'                     covariates = c("x1", "x2"),
#'                     candidates = c("z_shadow", "z_policy", "z_score"))
#' sc
#' plot(sc)
shadow_screen <- function(data, outcome, response, covariates = character(),
                          candidates, family = stats::binomial(),
                          min_relevance = 10, min_ratio = 1) {
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
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  weak <- out$relevance_chisq_df < min_relevance
  viol <- out$ratio <= min_ratio
  out$verdict <- ifelse(weak | viol, "reject", "pass")
  out$reason <- ifelse(weak & viol, "weak and violates exclusion",
                ifelse(weak, "weak", ifelse(viol, "violates exclusion", "")))
  attr(out, "settings") <- list(min_relevance = min_relevance,
                                min_ratio = min_ratio)
  class(out) <- c("shadow_screen", "data.frame")
  out
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
  chisq <- max(red[["deviance"]] - full[["deviance"]], 0)
  df <- full[["rank"]] - red[["rank"]]
  if (df < 1)
    stop("Candidate `", z, "` adds no degrees of freedom (constant or ",
         "collinear with the covariates).", call. = FALSE)
  list(r2 = chisq / red[["null"]], chisq = chisq, df = df)
}

#' @export
print.shadow_screen <- function(x, digits = 3, ...) {
  s <- attr(x, "settings")
  cat("<shadow_screen> relevance vs exclusion violation, incremental McFadden R^2\n")
  cat(sprintf("  pass: relevance chi2/df >= %s and ratio > %s\n",
              format(s$min_relevance), format(s$min_ratio)))
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(df, "settings") <- NULL
  show <- df[c("candidate", "relevance", "exclusion", "ratio",
               "relevance_chisq_df", "verdict", "reason")]
  num <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], signif, digits = digits)
  print(show, row.names = FALSE)
  invisible(x)
}
