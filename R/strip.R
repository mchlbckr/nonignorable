# Drop the row-level vectors of a scan, keeping what print() and plot() need,
# so that results on large data can be stored and shipped small.
strip_scan <- function(x) {
  x[c("eta_obs", "r", "weights", "binding")] <- list(NULL)
  x
}
