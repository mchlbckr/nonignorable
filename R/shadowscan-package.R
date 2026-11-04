#' shadowscan: screens and sharp scans for nonignorable binary outcomes
#'
#' If a binary outcome is missing not at random, every model fitted on the
#' observed rows silently asserts *one* value for how observability depends on
#' the outcome: none. Nothing in the data selects that value. The package
#' offers two tools for two different questions.
#'
#' * [shadow_screen()]: *is this variable fit to serve as a shadow variable?*
#'   Relevance for the outcome against violation of the exclusion restriction,
#'   on one scale. This belongs to a point-identified shadow-variable model.
#' * [shadow_scan()]: *how far can the answer move?* The sharp identified set
#'   under a bound on the logit tilt of the selection mechanism, with no shadow
#'   variable at all.
#'
#' The two rest on different assumptions and are not two halves of one
#' procedure. The help pages state the formulas, so that every number can be
#' checked by hand. [simulate_shadow()] provides data with known truth.
#'
#' @keywords internal
"_PACKAGE"
