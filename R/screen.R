#' Screen candidate shadow variables
#'
#' `shadow_screen()` asks whether a variable \eqn{Z} is fit to serve as a
#' **shadow variable** for an outcome \eqn{Y} that is missing not at
#' random, binary, categorical or numeric: relevant for \eqn{Y} given \eqn{X}, and excluded from the response
#' model given \eqn{X} and \eqn{Y}. It is a different tool from
#' [sharp_scan()]: the screen presupposes a shadow-variable model that
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
#' `ratio_readable`, reported by `print()` as a "violation at noise level"
#' and drawn as hollow points in the plot. For the candidate this is good news:
#' no violation is measured. Only the ratio is not an estimate. In large
#' samples every violation is significant, so exclusion should be read as an
#' effect size, not a test.
#'
#' Exclusion is measured *marginally* over \eqn{Y}: the pooled response model
#' cannot condition on an outcome that is missing. Under nonignorable selection
#' a relevant \eqn{Z} therefore shows some exclusion signal through \eqn{Y}
#' even when it is validly excluded; the screen is a screen, not a test of the
#' restriction. The ratio therefore counts only together with the relevance
#' requirement, never on its own: the stronger the selection on the outcome
#' and the association, the closer a validly excluded candidate moves to the
#' equality line (Bücker 2026, section "The screen when the outcome is
#' continuous", measures this on a design where exclusion holds exactly), and
#' for a candidate without relevance the ratio is noise
#' that shifts with the coding. The equality line is not the border between
#' usable and unusable candidates; a ratio near one does not by itself show a
#' violation, a ratio well below one does.
#'
#' **A numeric outcome.** When the observed outcome is numeric with more than
#' two values and `family` is left at its binomial default, relevance comes
#' from a proportional-odds model ([MASS::polr()]) on `outcome_classes`
#' classes of about equal mass of the observed outcome, or on its values if
#' it has no more; exclusion is unchanged, and is fitted once per candidate.
#' Adding \eqn{Z} costs its own degrees of freedom whatever the number of
#' classes, while the null deviance, twice \eqn{n} times the entropy of the
#' class shares, grows with it (up to \eqn{2n\log k} for \eqn{k} classes). So \eqn{\Delta R^2} *falls* as the
#' classes get finer while \eqn{\chi^2}/df rises towards the continuous
#' association: compare candidates under one coding only, and read
#' \eqn{\chi^2}/df, not \eqn{\Delta R^2}, across codings. Five classes is a
#' choice for the stability of the ordinal model (few cut points, the
#' proportional-odds assumption over few cuts), not an optimum. Bücker (2026),
#' section "The screen when the outcome is continuous", shows the behaviour
#' of both measures in the number of classes and checks the screen against a
#' known truth.
#'
#' **A categorical outcome.** A factor or character outcome with more than two
#' observed levels is screened in the same way, with its levels as classes:
#' an ordered factor in the proportional-odds model, an unordered one in a
#' multinomial logit ([nnet::multinom()]), where a candidate costs its
#' degrees of freedom once per level beyond the first. A factor with two
#' observed levels is binary.
#'
#' The screen **measures relevance; it does not certify identification**.
#' Identification with a shadow variable needs more than dependence of
#' \eqn{Z} on \eqn{Y} given \eqn{X}: in a parametric tilt model a rank
#' condition, in a nonparametric one completeness, which can fail in simple
#' models once the outcome or the candidate is continuous; Beppu and Morikawa
#' (2024) give a condition for categorical candidates that can be checked on
#' the observed data. The screen checks neither. And under selection on the
#' outcome a relevant candidate shows exclusion signal through \eqn{Y} even
#' when it is validly excluded, the more so the stronger the selection and
#' the association.
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
#'   model is always binomial. With the binomial default a numeric outcome
#'   with more than two values is screened by an ordinal model on classes,
#'   and a factor outcome with more than two levels by an ordinal or a
#'   multinomial model; see the sections on a numeric and a categorical
#'   outcome. Another family (say `poisson()` for
#'   counts) fits that `glm` instead.
#' @param min_relevance Minimum relevance \eqn{\chi^2} per degree of freedom
#'   for a candidate to pass. Default 10.
#' @param min_ratio Minimum ratio of relevance to exclusion for a candidate to
#'   pass, applied only together with `min_relevance`. Default 1, the
#'   equality line; see the details for why it is no test of exclusion.
#' @param numeric `"linear"` (default) or `"classes"`; see the section on
#'   numeric characteristics.
#' @param k Target number of classes for `numeric = "classes"`. Default 10.
#' @param outcome_classes Number of classes of about equal mass for a
#'   numeric outcome, cut only between distinct values as in
#'   [screen_classes()]. Default 5. An outcome with at most this many distinct
#'   values (a count, a rating) uses them as its levels; heavy atoms can leave
#'   fewer classes than asked, and `print()` says so.
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
#'   the heaping share of [screen_classes()]. Attribute `"outcome_coding"`:
#'   `model` (`"glm"`, `"ordinal"` or `"multinomial"`), the number of outcome `classes`, the
#'   number `requested` and the number of `distinct` observed values.
#' @references
#' Bücker, M. (2026). The price of a weak shadow variable: inference under
#' nonignorable missingness as relevance vanishes. \doi{10.48550/arXiv.2610.14695}
#'
#' Beppu, K. and Morikawa, K. (2024). Verifiable identification condition for
#' nonignorable nonresponse data with categorical instrumental variables.
#' *Statistical Theory and Related Fields*, 8, 40--50.
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
#' # a continuous outcome: relevance from an ordinal model on five classes
#' dc <- simulate_shadow(3000, tilt = -1, outcome = "continuous", seed = 1)
#' dc$y[dc$observed == 0] <- NA
#' shadow_screen(dc, "y", "observed", c("x1", "x2"), c("z_shadow", "z_policy"))
shadow_screen <- function(data, outcome, response, covariates = character(),
                          candidates, family = stats::binomial(),
                          min_relevance = 10, min_ratio = 1,
                          numeric = c("linear", "classes"), k = 10,
                          outcome_classes = 5) {
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

  # a numeric outcome with more than two values, under the default binomial
  # family: relevance from an ordinal model on classes of the outcome
  y_obs <- data[[outcome]][rv == 1 & !is.na(rv)]
  ordinal <- is.numeric(y_obs) && length(unique(y_obs)) > 2L &&
    identical(family$family, "binomial")
  outcome_coding <- list(model = "glm", classes = NA_integer_)
  # a factor or character outcome with more than two observed levels: an
  # ordered factor keeps its order in the ordinal model, an unordered one goes
  # into a multinomial model; glm() would silently fit the first level
  # against all others
  categorical <- (is.factor(y_obs) || is.character(y_obs)) &&
    length(unique(y_obs)) > 2L && identical(family$family, "binomial")
  if (categorical) {
    yv <- data[[outcome]]
    yv[rv != 1 | is.na(rv)] <- NA
    yv <- if (is.ordered(yv)) droplevels(yv) else factor(yv)
    data[[outcome]] <- yv
    outcome_coding <- list(model = if (is.ordered(yv)) "ordinal" else "multinomial",
                           classes = nlevels(yv), requested = NA_integer_,
                           distinct = nlevels(yv))
  }
  if (ordinal) {
    if (!is.numeric(outcome_classes) || length(outcome_classes) != 1L ||
        is.na(outcome_classes) || outcome_classes < 2 ||
        outcome_classes != round(outcome_classes))
      stop("`outcome_classes` must be a single whole number >= 2.", call. = FALSE)
    # few distinct values: they are already ordinal levels; otherwise classes
    # of about equal mass, cut only between distinct values, so that an atom
    # (a count of 0, a rounded score) is never merged with its neighbour
    yv <- data[[outcome]]
    yv[rv != 1 | is.na(rv)] <- NA
    vals <- sort(unique(y_obs))
    yc <- if (length(vals) <= outcome_classes) factor(yv, levels = vals)
          else screen_classes(yv, k = outcome_classes, support = range(y_obs))
    data[[outcome]] <- factor(yc, levels = levels(yc), ordered = TRUE)
    outcome_coding <- list(model = "ordinal", classes = nlevels(data[[outcome]]),
                           requested = as.integer(outcome_classes),
                           distinct = length(vals))
  }

  coded <- code_numeric(data, response, c(covariates, candidates), numeric, k)
  data <- coded$data

  rows <- lapply(candidates, function(z) {
    cols <- c(response, covariates, z)
    # keep only the columns needed, so that subsetting a wide frame stays cheap
    pooled <- data[stats::complete.cases(data[cols]), c(outcome, cols),
                   drop = FALSE]
    obs <- pooled[pooled[[response]] == 1, , drop = FALSE]
    rel <- switch(outcome_coding$model,
                  ordinal = incr_polr(obs, outcome, covariates, z),
                  multinomial = incr_multinom(obs, outcome, covariates, z),
                  incr_fit(obs, outcome, covariates, z, family))
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
  attr(out, "outcome_coding") <- outcome_coding
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

# Incremental McFadden R^2 of z in a proportional-odds model for an ordered
# outcome. The null deviance is that of the class shares, 2 n log k for k
# classes of equal mass; df counts the columns z adds, whatever the number of
# classes, since the k - 1 cut points sit in both models.
incr_polr <- function(d, lhs, covariates, z) {
  d[[lhs]] <- droplevels(d[[lhs]])
  tab <- table(d[[lhs]]); tab <- tab[tab > 0]
  null <- -2 * sum(tab * log(tab / sum(tab)))
  rank_of <- function(vars) if (!length(vars)) 1L else
    qr(stats::model.matrix(stats::reformulate(vars), d))$rank
  dev_of <- function(vars) if (!length(vars)) null else
    MASS::polr(stats::reformulate(vars, lhs), data = d)$deviance
  df <- rank_of(c(covariates, z)) - rank_of(covariates)
  if (df < 1)
    stop("Candidate `", z, "` adds no degrees of freedom (constant or ",
         "collinear with the covariates).", call. = FALSE)
  chisq <- max(dev_of(covariates) - dev_of(c(covariates, z)), 0)
  list(r2 = chisq / null, chisq = chisq, df = df, epp = NA_real_)
}

# The same for an unordered outcome in a multinomial logit: z costs its
# columns once per non-reference level, (J - 1) * dim(z) degrees of freedom.
incr_multinom <- function(d, lhs, covariates, z) {
  d[[lhs]] <- droplevels(d[[lhs]])
  tab <- table(d[[lhs]]); tab <- tab[tab > 0]
  null <- -2 * sum(tab * log(tab / sum(tab)))
  rank_of <- function(vars) if (!length(vars)) 1L else
    qr(stats::model.matrix(stats::reformulate(vars), d))$rank
  dev_of <- function(vars) if (!length(vars)) null else
    nnet::multinom(stats::reformulate(vars, lhs), data = d, trace = FALSE,
                   maxit = 1000)$deviance
  df <- (rank_of(c(covariates, z)) - rank_of(covariates)) * (length(tab) - 1L)
  if (df < 1)
    stop("Candidate `", z, "` adds no degrees of freedom (constant or ",
         "collinear with the covariates).", call. = FALSE)
  chisq <- max(dev_of(covariates) - dev_of(c(covariates, z)), 0)
  list(r2 = chisq / null, chisq = chisq, df = as.integer(df), epp = NA_real_)
}

# "pass", or "reject: " and the reasons in short form, for one column of print
verdict_text <- function(verdict, reason) {
  short <- gsub("violates exclusion", "exclusion", gsub(" and ", ", ", reason))
  ifelse(verdict == "reject" & nzchar(short), paste0("reject: ", short), verdict)
}

#' @export
print.shadow_screen <- function(x, digits = 3, ...) {
  s <- attr(x, "settings")
  cat("<shadow_screen> relevance vs exclusion violation, incremental McFadden R^2\n")
  cat(sprintf("  pass: relevance chi2/df >= %s and ratio > %s\n",
              format(s$min_relevance), format(s$min_ratio)))
  oc <- attr(x, "outcome_coding")
  if (identical(oc$model, "multinomial"))
    cat(strwrap(sprintf(paste0("relevance: multinomial model on %d classes of the observed ",
                               "outcome; a candidate has (classes - 1) degrees of freedom ",
                               "per column, so read chi2/df"), oc$classes),
                indent = 2, exdent = 4), sep = "\n")
  if (identical(oc$model, "ordinal") && is.na(oc$requested))
    cat(strwrap(sprintf("relevance: ordinal model on the %d levels of the observed outcome",
                        oc$classes), indent = 2, exdent = 4), sep = "\n")
  if (identical(oc$model, "ordinal") && !is.na(oc$requested))
    cat(strwrap(sprintf(paste0("relevance: ordinal model on %d classes of the observed ",
                               "outcome%s; R^2 falls as classes get finer, so compare chi2/df ",
                               "across codings"),
                        oc$classes,
                        if (oc$distinct <= oc$requested) " (its distinct values)"
                        else if (oc$classes != oc$requested)
                          sprintf(" (%d requested; atoms are not split)", oc$requested)
                        else ""), indent = 2, exdent = 4), sep = "\n")
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
  # each number formatted on its own, so that one large value does not pad
  # the whole column with zeros; short header, so the table fits 80 columns
  num <- list(relevance = df$relevance, exclusion = df$exclusion, ratio = df$ratio,
              "chi2/df" = df$relevance_chisq_df)
  num <- lapply(num, function(v) vapply(v, function(z) format(z, digits = digits), ""))
  wd <- pmax(nchar(names(num)), vapply(num, function(v) max(nchar(v)), 1L))
  show <- data.frame(df$candidate, Map(formatC, num, width = wd),
                     verdict_text(df$verdict, df$reason),
                     check.names = FALSE, stringsAsFactors = FALSE)
  names(show) <- c("candidate", mapply(formatC, names(num), width = wd), "verdict")
  print(show, row.names = FALSE, right = FALSE)
  if (!is.null(x$ratio_readable) && any(!x$ratio_readable))
    cat(strwrap(sprintf(paste0("violation at noise level (exclusion chi2/df below 1; ",
                               "the ratio divides by noise): %s"),
                        paste(x$candidate[!x$ratio_readable], collapse = ", ")),
                indent = 2, exdent = 4), sep = "\n")
  invisible(x)
}
