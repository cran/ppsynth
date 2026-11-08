# ppsynth.R -- long-panel front end, mirroring the Stata command -ppsynth-
# (ppsynth.ado v2.6.1): option set, panel construction, per-period effects
# and p-values, and the returned object.

#' Partial-pooling synthetic control with multiple outcomes
#'
#' R port of the Stata command \code{ppsynth} (and of the Python engine it
#' calls). Estimates synthetic-control weights for two or more outcomes with
#' partial pooling: each outcome gets its own simplex weight vector, shrunk
#' toward a common centroid (the cross-outcome mean weight vector) by a fusion
#' penalty of strength lambda. lambda = 0 reproduces separate per-outcome SCM;
#' lambda = Inf reproduces common (concatenated) weights as in Sun, Ben-Michael
#' and Feller (REStat 2025). Both corners are computed exactly; interior grid
#' values are computed by block-coordinate descent, warm-started along the
#' grid.
#'
#' @param data a data frame holding a balanced long panel.
#' @param outcomes character vector of outcome column names (K >= 1; the
#'   method is designed for K >= 2).
#' @param id name of the numeric panel-unit identifier column.
#' @param time name of the numeric time column.
#' @param trunit value of \code{id} identifying the treated unit.
#' @param trperiod first post-treatment value of \code{time}.
#' @param method \code{"aggregate"} (default) or \code{"cv"}; see Details.
#' @param lambda fix the pooling strength instead of tuning it; pass
#'   \code{Inf} for the common-weights corner. Snapped to the nearest grid
#'   point. Overrides \code{method}.
#' @param holdout holdout periods for \code{method = "cv"}; default
#'   \code{max(2, T0 \%/\% 4)}.
#' @param inference logical; run whole-pipeline placebo inference (reruns the
#'   entire pipeline, tuning included, for every donor treated as a placebo).
#'   Slow in this pure-R engine: see the README for timings.
#' @param convention p-value convention, \code{"phipson_smyth"} (default) or
#'   \code{"strict"}; see Details.
#' @param likesynth logical; synth/synth_runner presentation: implies
#'   \code{convention = "strict"} and draws the treatment line in plots at the
#'   first treated period rather than the last pre-treatment period.
#' @param lam_grid the lambda grid; default \code{\link{LAM_GRID_DEFAULT}}.
#' @param engine \code{"cpp"} (default; compiled inner solver) or \code{"r"}
#'   (pure-R reference solver). Same algorithm and settings; both are
#'   certified against the Python engine behind the Stata command.
#' @param cores number of cores to use (default 1): the leave-one-out folds
#'   of the aggregation step and, with \code{inference = TRUE}, the placebo
#'   units are independent and are spread across cores. Results do not
#'   depend on \code{cores}.
#' @param verbose logical; print the inference progress notes.
#'
#' @details
#' Two tuning methods, selected by \code{method}:
#' \describe{
#' \item{\code{"aggregate"} (default)}{exponentially weighted average over the
#'   lambda-grid "experts", each scored by leave-one-out prediction on the
#'   pre-period. The reported weights are the mass-weighted average of expert
#'   weights, a convex combination and hence itself on the simplex. This is
#'   the paper's primary estimator; it carries a pathwise regret guarantee
#'   relative to every expert on the grid, including both corners. The
#'   default display reports the final estimator's realized weight pooling
#'   and the middle 50\% of predictive mass on the realized-pooling scale;
#'   \code{print(x, massprofile = TRUE)} adds the mass over the grid.}
#' \item{\code{"cv"}}{leave-future-out cross-validation selects a single lambda
#'   from the grid. A companion/diagnostic, not the primary estimator: its
#'   loss is summed across outcomes with no down-weighting, so a single
#'   poorly-fitting outcome can dominate the selection (see the fittability
#'   guardrail).}
#' }
#' \code{lambda} overrides both and fixes lambda directly.
#'
#' Donor-support check: for each outcome, the number of pre-treatment periods
#' in which the treated unit's value falls outside the range spanned by the
#' donors. Nonnegative weights summing to one cannot reach outside that
#' range. Runs on every fit, fires on any violation, report-only; the caution
#' is printed and returned in \code{support_warning}.
#'
#' Fittability guardrail: if any outcome's pre-RMSPE exceeds 0.9 times its own
#' pre-period standard deviation, that outcome is essentially unfittable from
#' this donor pool. Under \code{method = "cv"} a caution naming the outcome(s)
#' is printed and returned in \code{cv_warning}.
#'
#' \code{inference = TRUE} reruns the entire chosen pipeline for each donor
#' treated in turn as a placebo, the actual treated unit returned to the donor
#' pool. Per-outcome p-values compare the treated unit's post/pre RMSPE ratio
#' against the placebo distribution; joint p-values use the mean and max of
#' the ratios studentized over all N0+1 units. \code{convention}:
#' \describe{
#' \item{\code{"phipson_smyth"} (default)}{p = (1 + count of placebo ratios at
#'   or above the treated ratio) / (N0+1). Always valid (Phipson & Smyth
#'   2010); minimum attainable value 1/(N0+1).}
#' \item{\code{"strict"}}{p = (count strictly above) / (N0+1), the synth-style
#'   strict-exceedance count. Can be exactly zero and is anti-conservative;
#'   provided for comparability with legacy synth output.}
#' }
#' Per-period placebo p-values (\code{pp_pvals}) rank, at each post period, the
#' treated unit's absolute gap standardized by its own pre-RMSPE among the
#' placebos' (descriptive; shown by \code{print(x, effects = TRUE)}).
#'
#' @return An object of class \code{"ppsynth"}: a list with the quantities the
#'   Stata command returns in \code{r()}, under the same names where
#'   possible: \code{gamma} (N0 x K donor weights, rows named by donor id),
#'   \code{gamma_common}, \code{tau}, \code{pre_rmspe}, \code{effects}
#'   (T1 x K post-period effects), \code{pp_pvals} (with inference),
#'   \code{lam_grid}, \code{agg_mass}, \code{cv_curve}, \code{p_outcome},
#'   \code{p_joint_max}, \code{p_joint_mean}, \code{y_treated},
#'   \code{y_synth}, \code{att_path} (all T x K), \code{times},
#'   \code{plgaps} (list of K T x N0 placebo gap matrices, with inference),
#'   \code{lambda}, \code{T0}, \code{T1}, \code{N0}, \code{method},
#'   \code{outcomes}, \code{cv_warning}, \code{support_warning},
#'   \code{realized_pooling}, \code{pooling_q25}, \code{pooling_q75},
#'   \code{pooling_grid}, \code{pooling_degenerate}, \code{convention},
#'   plus \code{fit} (the raw \code{\link{fit_ppsynth}} result) and
#'   \code{inference} (the raw \code{\link{placebo_inference}} result or
#'   \code{NULL}).
#'
#' @references Grier, K. B. (2026). "Correlated Weights: Partial Pooling for
#'   Synthetic Control with Multiple Outcomes." Working paper,
#'   \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7149398}.
#'
#'   Sun, L., Ben-Michael, E. and Feller, A. (2025). "Using Multiple Outcomes
#'   to Improve the Synthetic Control Method." Review of Economics and
#'   Statistics.
#'
#'   Phipson, B. and Smyth, G. K. (2010). "Permutation p-values should never
#'   be zero." Statistical Applications in Genetics and Molecular Biology 9(1).
#'
#' @examples
#' # a small simulated panel shipped with the package (unit 1 treated from 2010)
#' data(ppsynth_test_panel)
#' fit <- ppsynth(ppsynth_test_panel, outcomes = c("gdp", "health", "invest"),
#'                id = "id", time = "year", trunit = 1, trperiod = 2010)
#' print(fit, massprofile = TRUE)
#' fit$gamma
#'
#' # the Stata worked example (Section 5A of the paper): German reunification,
#' # five outcomes, 20 donors. The aggregate fit takes a few seconds.
#' \donttest{
#' data(germany_panel)
#' outs <- c("proprights", "money", "inttrade", "regulation", "csh_g")
#' fit <- ppsynth(germany_panel, outcomes = outs, id = "id", time = "year",
#'                trunit = 8, trperiod = 1990)
#' print(fit, massprofile = TRUE)
#' plot(fit, which = "paths")
#' }
#' # Placebo inference reruns the whole pipeline for every donor (a minute or
#' # two on this panel; `cores` spreads the donors across cores):
#' \donttest{
#' inf <- ppsynth(germany_panel, outcomes = outs, id = "id", time = "year",
#'                trunit = 8, trperiod = 1990, inference = TRUE, cores = 2)
#' print(inf, effects = TRUE)
#' plot(inf, which = "gaps")
#' }
#' @seealso \code{\link{print.ppsynth}}, \code{\link{plot.ppsynth}},
#'   \code{\link{fit_ppsynth}}, \code{\link{placebo_inference}}
#' @export
ppsynth <- function(data, outcomes, id, time, trunit, trperiod,
                    method = c("aggregate", "cv"), lambda = NULL,
                    holdout = NULL, inference = FALSE,
                    convention = c("phipson_smyth", "strict"),
                    likesynth = FALSE, lam_grid = NULL,
                    engine = c("cpp", "r"), cores = 1L, verbose = TRUE) {
  method <- match.arg(method)
  engine <- match.arg(engine)
  conv_given <- !missing(convention)
  convention <- match.arg(convention)
  if (likesynth) {
    if (conv_given && convention != "strict")
      stop("likesynth implies convention = \"strict\"; drop one of the two")
    convention <- "strict"
  }
  data <- as.data.frame(data)
  for (v in c(id, time, outcomes))
    if (!v %in% names(data)) stop(sprintf("column '%s' not found in data", v))

  # markout: drop rows with missing id/time/outcomes (Stata's markout)
  keep <- stats::complete.cases(data[, c(id, time, outcomes), drop = FALSE])
  data <- data[keep, , drop = FALSE]

  ids <- as.numeric(data[[id]])
  times <- as.numeric(data[[time]])
  if (!any(ids == trunit)) stop(sprintf("treated unit %s not found in %s", format(trunit), id))

  utime <- sort(unique(times))
  uids <- sort(unique(ids))
  donors <- uids[uids != trunit]
  T0 <- sum(utime < trperiod)
  Tn <- length(utime)
  K <- length(outcomes)

  Y_tr <- vector("list", K); Y_dn <- vector("list", K)
  trow <- match(times, utime)
  dcol <- match(ids, donors)
  is_tr <- ids == trunit
  for (k in seq_len(K)) {
    y <- as.numeric(data[[outcomes[k]]])
    tr <- rep(NA_real_, Tn)
    dn <- matrix(NA_real_, Tn, length(donors))
    tr[trow[is_tr]] <- y[is_tr]
    dn[cbind(trow[!is_tr], dcol[!is_tr])] <- y[!is_tr]
    if (anyNA(tr) || anyNA(dn))
      stop(sprintf("ppsynth: unbalanced panel or missing values in %s", outcomes[k]))
    Y_tr[[k]] <- tr; Y_dn[[k]] <- dn
  }

  inf <- NULL
  if (inference) {
    if (verbose) {
      message(sprintf("  Running placebo inference: %d full pipeline re-estimations.",
                      length(uids)))
      message("  Each re-estimation reruns the whole pipeline (tuning included).")
    }
    t_start <- proc.time()[["elapsed"]]
    inf <- placebo_inference(Y_tr, Y_dn, T0, method = method, lam = lambda,
                             lam_grid = lam_grid, holdout = holdout,
                             convention = convention, outcome_names = outcomes,
                             verbose = verbose, engine = engine, cores = cores)
    el <- round(proc.time()[["elapsed"]] - t_start)
    if (verbose) {
      fmt <- if (el >= 60) sprintf("%dm %ds", el %/% 60, el %% 60) else sprintf("%ds", el)
      message(sprintf("  placebo inference completed in %s.", fmt))
    }
    res <- inf$treated_fit
  } else {
    res <- fit_ppsynth(Y_tr, Y_dn, T0, method = method, lam = lambda,
                       lam_grid = lam_grid, holdout = holdout,
                       outcome_names = outcomes, engine = engine, cores = cores)
  }

  dn_names <- vapply(donors, function(dd)
    if (dd == round(dd)) as.character(as.integer(dd)) else as.character(dd), "")
  gamma <- res$gamma
  dimnames(gamma) <- list(dn_names, outcomes)
  gamma_common <- stats::setNames(res$gamma_common, dn_names)
  att <- do.call(cbind, res$att_path); colnames(att) <- outcomes
  ytr <- do.call(cbind, Y_tr); colnames(ytr) <- outcomes
  rownames(att) <- rownames(ytr) <- utime
  post <- (T0 + 1):Tn

  out <- list(
    gamma = gamma, gamma_common = gamma_common,
    tau = stats::setNames(res$tau, outcomes),
    pre_rmspe = stats::setNames(res$pre_rmspe, outcomes),
    effects = att[post, , drop = FALSE],
    pp_pvals = NULL,
    lam_grid = res$lam_grid, agg_mass = res$agg_mass, cv_curve = res$cv_curve,
    p_outcome = NULL, p_joint_max = NULL, p_joint_mean = NULL,
    y_treated = ytr, y_synth = ytr - att, att_path = att, times = utime,
    plgaps = NULL,
    lambda = res$lam, T0 = T0, T1 = Tn - T0, N0 = length(donors),
    method = method, outcomes = outcomes,
    cv_warning = res$cv_warning, support_warning = res$support_warning,
    realized_pooling = res$realized_pooling,
    pooling_q25 = res$pooling_q25, pooling_q75 = res$pooling_q75,
    pooling_grid = res$pooling_grid, pooling_degenerate = res$pooling_degenerate,
    convention = convention, likesynth = likesynth,
    trunit = trunit, trperiod = trperiod, donors = donors,
    engine = engine, fit = res, inference = inf)

  if (inference) {
    out$p_outcome <- stats::setNames(inf$p_outcome, outcomes)
    out$p_joint_max <- inf$p_joint_max
    out$p_joint_mean <- inf$p_joint_mean
    # placebo gap paths, one T x N0 matrix per outcome (donor order = column order)
    plgaps <- vector("list", K)
    pp <- matrix(NA_real_, Tn - T0, K, dimnames = list(utime[post], outcomes))
    for (k in seq_len(K)) {
      G <- do.call(cbind, lapply(inf$placebo_gaps, function(pg) pg[[k]]))
      dimnames(G) <- list(utime, dn_names)
      plgaps[[k]] <- G
      # per-period p-values: |standardized gap| of the treated ranked among
      # the placebos', each unit standardized by its own pre-RMSPE
      prem <- sqrt(colMeans(G[seq_len(T0), , drop = FALSE]^2))
      prem[prem <= 1e-12] <- NA_real_
      ztr <- abs(att[post, k]) / max(res$pre_rmspe[k], 1e-12)
      zpl <- sweep(abs(G[post, , drop = FALSE]), 2, prem, "/")
      J <- ncol(zpl)
      if (convention == "strict") {
        pp[, k] <- rowSums(zpl > ztr, na.rm = TRUE) / (J + 1)
      } else {
        pp[, k] <- (1 + rowSums(zpl >= ztr, na.rm = TRUE)) / (J + 1)
      }
    }
    names(plgaps) <- outcomes
    out$plgaps <- plgaps
    out$pp_pvals <- pp
  }
  class(out) <- "ppsynth"
  out
}
