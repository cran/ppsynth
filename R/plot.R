# plot.R -- base-graphics counterparts of the Stata -fig- graphs:
# mass profile (aggregate), CV curve (cv), treated-vs-synthetic paths, and
# gap paths with placebo donors (inference).

#' Plot a ppsynth fit
#'
#' Base-graphics versions of the Stata command's \code{fig} graphs. Each
#' selected graph is drawn on its own page; wrap the call in \code{pdf()} or
#' \code{png()} to save.
#'
#' @param x a \code{"ppsynth"} object.
#' @param which which graphs to draw: any of \code{"paths"} (treated vs.
#'   synthetic, one panel per outcome), \code{"gaps"} (treated gap against
#'   the placebo donors' gaps; needs \code{inference = TRUE}), \code{"mass"}
#'   (aggregation mass over the lambda grid; aggregate method) and
#'   \code{"cv"} (leave-future-out loss over the grid; cv method). The
#'   default draws every graph that applies, in the Stata order.
#' @param likesynth logical; draw the treatment line at the first treated
#'   period (synth presentation) instead of the last pre-treatment period.
#'   Defaults to the value used in the fit.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
plot.ppsynth <- function(x, which = c("mass", "cv", "paths", "gaps"),
                         likesynth = x$likesynth, ...) {
  which <- match.arg(which, several.ok = TRUE)
  doinf <- !is.null(x$plgaps)
  tm <- x$times
  pre_times <- tm[tm < x$trperiod]
  if (!likesynth && length(pre_times)) {
    xlpos <- max(pre_times); lnote <- "vertical line = last pre-treatment period"
  } else {
    xlpos <- x$trperiod; lnote <- "vertical line = first treated period"
  }
  K <- length(x$outcomes)
  grid_labels <- vapply(x$lam_grid, function(l)
    if (is.infinite(l) || l >= 1e6) "inf" else trimws(formatC(l, digits = 3, format = "g")), "")

  if ("mass" %in% which && x$method == "aggregate" && !is.null(x$agg_mass)) {
    op <- graphics::par(mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op), add = TRUE)
    graphics::barplot(x$agg_mass, names.arg = grid_labels, las = 2,
                      ylab = "EWA mass", xlab = "", col = "grey70", border = NA,
                      main = "Aggregation mass over the lambda grid")
    graphics::mtext("lambda", side = 1, line = 3)
  }
  if ("cv" %in% which && x$method == "cv" && !is.null(x$cv_curve)) {
    op <- graphics::par(mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op), add = TRUE)
    M <- length(x$cv_curve)
    graphics::plot(seq_len(M), x$cv_curve, type = "b", pch = 16, xaxt = "n",
                   ylab = "leave-future-out loss", xlab = "lambda",
                   main = "CV loss over the lambda grid")
    graphics::axis(1, at = seq_len(M), labels = grid_labels, las = 2)
  }
  if ("paths" %in% which) {
    nr <- ceiling(sqrt(K)); nc <- ceiling(K / nr)
    op <- graphics::par(mfrow = c(nr, nc), mar = c(3, 3, 2.5, 1), oma = c(2.5, 0, 2.5, 0))
    on.exit(graphics::par(op), add = TRUE)
    for (k in seq_len(K)) {
      yr <- range(c(x$y_treated[, k], x$y_synth[, k]))
      graphics::plot(tm, x$y_treated[, k], type = "l", lwd = 2.2, col = "navy",
                     ylim = yr, xlab = "", ylab = "", main = x$outcomes[k], cex.main = 0.95)
      graphics::lines(tm, x$y_synth[, k], lty = 2, col = "maroon", lwd = 1.6)
      graphics::abline(v = xlpos, lty = 3, col = "grey50")
    }
    graphics::mtext("Treated (solid) vs. synthetic (dash)", outer = TRUE, side = 3, line = 0.8, font = 2)
    graphics::mtext(lnote, outer = TRUE, side = 1, line = 1, cex = 0.8)
  }
  if ("gaps" %in% which) {
    if (!doinf) {
      message("  (gap/placebo plot requires inference = TRUE; skipped)")
    } else {
      nr <- ceiling(sqrt(K)); nc <- ceiling(K / nr)
      op <- graphics::par(mfrow = c(nr, nc), mar = c(3, 3, 2.5, 1), oma = c(2.5, 0, 2.5, 0))
      on.exit(graphics::par(op), add = TRUE)
      for (k in seq_len(K)) {
        G <- x$plgaps[[k]]
        yr <- range(c(G, x$att_path[, k]))
        graphics::matplot(tm, G, type = "l", lty = 1, lwd = 0.6, col = "grey80",
                          ylim = yr, xlab = "", ylab = "", main = x$outcomes[k], cex.main = 0.95)
        graphics::lines(tm, x$att_path[, k], lwd = 2.2, col = "navy")
        graphics::abline(h = 0, lty = 2, col = "grey50")
        graphics::abline(v = xlpos, lty = 3, col = "grey50")
      }
      graphics::mtext("Gap paths: treated (navy) vs. placebo donors (gray)", outer = TRUE,
                      side = 3, line = 0.8, font = 2)
      graphics::mtext(paste("whole-pipeline placebo gaps.", lnote), outer = TRUE,
                      side = 1, line = 1, cex = 0.8)
    }
  }
  invisible(x)
}
