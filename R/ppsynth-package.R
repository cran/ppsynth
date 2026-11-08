#' ppsynth: Partial-Pooling Synthetic Control with Multiple Outcomes
#'
#' R implementation of the estimator in Grier (2026). See
#' \code{\link{ppsynth}} for the long-panel interface, which mirrors the
#' Stata command of the same name, and \code{\link{fit_ppsynth}} /
#' \code{\link{placebo_inference}} for the matrix-level engine.
#'
#' Two engines are available and certified against each other and against
#' the Python engine behind the Stata command: \code{engine = "cpp"}
#' (default; compiled inner solver) and \code{engine = "r"} (pure R
#' reference). Both run the same algorithm with the same solver settings.
#'
#' @keywords internal
#' @useDynLib ppsynth, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
