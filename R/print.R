# print.R -- Stata-style display for ppsynth objects (mirrors ppsynth.ado's
# results block, including the 2.6.0 pooling headline, the -massprofile-
# table, the -effects- tables, and the cautions).

fmt_g <- function(x, width, digits) {
  # Stata %w.dg: general format with `digits` significant digits
  formatC(x, width = width, digits = digits, format = "g")
}
fmt_f <- function(x, width, digits) formatC(x, width = width, digits = digits, format = "f")
center <- function(s, width) {
  s <- as.character(s)
  if (nchar(s) >= width) return(s)
  left <- (width - nchar(s)) %/% 2
  paste0(strrep(" ", left), s, strrep(" ", width - nchar(s) - left))
}
lam_label <- function(l) if (is.infinite(l) || l >= 1e5) "inf" else trimws(formatC(l, digits = 3, format = "g"))
hline <- function(n = 64) strrep("-", n)

#' Print a ppsynth fit
#'
#' Reproduces the Stata command's results display.
#'
#' @param x a \code{"ppsynth"} object.
#' @param massprofile logical; also print the aggregation mass profile over
#'   the lambda grid (the Stata \code{massprofile} option).
#' @param effects logical; also print the per-period effects and, with
#'   inference, the per-period placebo p-values (the Stata \code{effects}
#'   option).
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.ppsynth <- function(x, massprofile = FALSE, effects = FALSE, ...) {
  K <- length(x$outcomes)
  doinf <- !is.null(x$p_outcome)
  cat("\n")
  cat(sprintf("Partial-pooling synthetic control (%d outcomes, method = %s)\n", K, x$method))
  cat(hline(), "\n", sep = "")
  cat(sprintf("  pre: %d   post: %d   donors: %d\n", x$T0, x$T1, x$N0))
  cat(hline(), "\n", sep = "")
  if (doinf) {
    cat(center("outcome", 14), formatC("mean post-treatment effect", width = 28),
        formatC("pre-RMSPE", width = 13), formatC("placebo p", width = 12), "\n", sep = "")
    for (k in seq_len(K)) {
      cat(center(x$outcomes[k], 14), fmt_g(x$tau[k], 28, 4), fmt_g(x$pre_rmspe[k], 13, 4),
          fmt_f(x$p_outcome[k], 12, 3), "\n", sep = "")
    }
    cat(hline(), "\n", sep = "")
    cat(sprintf("  joint p (mean stat) = %6.3f    joint p (max stat) = %6.3f\n",
                x$p_joint_mean, x$p_joint_max))
  } else {
    cat(center("outcome", 14), formatC("mean post-treatment effect", width = 28),
        formatC("pre-RMSPE", width = 13), "\n", sep = "")
    for (k in seq_len(K)) {
      cat(center(x$outcomes[k], 14), fmt_g(x$tau[k], 28, 4), fmt_g(x$pre_rmspe[k], 13, 4),
          "\n", sep = "")
    }
    cat(hline(), "\n", sep = "")
  }
  if (x$method == "aggregate" && !is.null(x$agg_mass)) {
    if (isTRUE(x$pooling_degenerate)) {
      cat("  weights already common at lambda = 0; pooling scale undefined.\n")
      cat(hline(), "\n", sep = "")
    } else if (!is.null(x$realized_pooling)) {
      cat(sprintf("  Realized weight pooling: %.1f%%\n", x$realized_pooling))
      if (!is.null(x$pooling_q25)) {
        cat(sprintf("  Middle 50%% of predictive mass: %.0f%% to %.0f%% pooling\n",
                    x$pooling_q25, x$pooling_q75))
      }
      cat(hline(), "\n", sep = "")
    }
    hints <- character(0)
    if (!doinf) hints <- c(hints, "inference = TRUE for placebo p-values")
    if (!massprofile) hints <- c(hints, "massprofile = TRUE for the lambda grid")
    if (length(hints)) cat(sprintf("  (add %s)\n", paste(hints, collapse = ", ")))
    if (massprofile && !isTRUE(x$pooling_degenerate)) {
      M <- length(x$agg_mass)
      labrow <- paste(formatC(vapply(x$lam_grid, lam_label, ""), width = 6), collapse = "")
      poolrow <- paste(formatC(trimws(formatC(x$pooling_grid, format = "f", digits = 0)), width = 6), collapse = "")
      massrow <- paste(formatC(trimws(formatC(100 * x$agg_mass, format = "f", digits = 1)), width = 6), collapse = "")
      cat("\n  Aggregation mass profile:\n")
      cat(formatC("lambda", width = 10), labrow, "\n", sep = "")
      cat(formatC("pool %", width = 10), poolrow, "\n", sep = "")
      cat(formatC("mass %", width = 10), massrow, "\n", sep = "")
      cat("    Note: pooling per lambda is panel-specific; see ?ppsynth.\n")
    }
  }
  if (effects) {
    hdr <- paste0("    period", paste(formatC(abbreviate(x$outcomes, 9), width = 10), collapse = ""))
    cat("\n  per-period effects (treated minus synthetic, original scale)\n")
    cat(hdr, "\n", sep = "")
    ptimes <- as.numeric(rownames(x$effects))
    for (i in seq_len(nrow(x$effects))) {
      cat(formatC(trimws(formatC(ptimes[i], digits = 10, format = "g")), width = 10),
          paste(fmt_g(x$effects[i, ], 10, 3), collapse = ""), "\n", sep = "")
    }
    if (doinf && !is.null(x$pp_pvals)) {
      cat(sprintf("\n  per-period placebo p-values (%s convention)\n", x$convention))
      cat(hdr, "\n", sep = "")
      for (i in seq_len(nrow(x$pp_pvals))) {
        cat(formatC(trimws(formatC(ptimes[i], digits = 10, format = "g")), width = 10),
            paste(fmt_f(x$pp_pvals[i, ], 10, 2), collapse = ""), "\n", sep = "")
      }
    }
  }
  if (!is.null(x$support_warning)) cat("\n  ", x$support_warning, "\n", sep = "")
  if (!is.null(x$cv_warning)) cat("\n  ", x$cv_warning, "\n", sep = "")
  invisible(x)
}
