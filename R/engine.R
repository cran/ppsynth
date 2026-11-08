# engine.R -- pure-R port of ppsynth_core.py (PPSYNTH v2.6.1, release 2.6.0
# engine). Partial-pooling synthetic control with multiple outcomes.
#
# This file is a line-by-line port of the Python reference engine. Function
# names, solver settings (iteration caps, tolerances, warm-start order), the
# eta rule, the pooling-scale conventions, and the warning texts are kept
# identical so that the R engine can be certified against the Python one.
#
# Estimator family: per-outcome simplex weights gamma_k with fusion penalty
#     sum_k (1/T0)||y_k - X_k g_k||^2 + lam * sum_k ||g_k - gbar||^2,
# gbar = mean_k g_k. lam = 0 is separate SCM; lam = Inf is common
# (concatenated) weights as in Sun, Ben-Michael & Feller (REStat 2025).
# Exact corner solutions are included as grid endpoints.

#' Default lambda grid
#'
#' The ten-point pooling grid used by the aggregation estimator: exact
#' separate fits at 0, exact common weights at \code{Inf}, and eight interior
#' penalties warm-started along the grid.
#' @export
LAM_GRID_DEFAULT <- c(0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, Inf)

#' Fittability threshold
#'
#' An outcome whose pre-treatment RMSPE exceeds this multiple of its own
#' pre-period standard deviation is flagged as poorly fit (see Section 5.3 of
#' the working paper).
#' @export
FITTABILITY_THRESHOLD <- 0.9

# ---------------------------------------------------------------------------
# solvers
# ---------------------------------------------------------------------------

# Euclidean projection onto the probability simplex.
project_simplex <- function(v) {
  n <- length(v)
  u <- sort(v, decreasing = TRUE)
  css <- cumsum(u)
  rho <- max(which(u + (1 - css) / seq_len(n) > 0))
  theta <- (css[rho] - 1) / rho
  pmax(v - theta, 0)
}

# Accelerated projected gradient (FISTA) for
#   min_x  x'Ax - 2b'x + lam ||x - center||^2   s.t. x in simplex.
fista <- function(A, b, lam, center, x0, Lmax, n_iter = 400L, tol = 1e-9) {
  step <- 1 / (2 * Lmax + 2 * lam + 1e-12)
  x <- x0
  y <- x0
  t <- 1
  for (it in seq_len(n_iter)) {
    grad <- 2 * (drop(A %*% y) - b)
    if (lam > 0) grad <- grad + 2 * lam * (y - center)
    xn <- project_simplex(y - step * grad)
    tn <- 0.5 * (1 + sqrt(1 + 4 * t * t))
    y <- xn + ((t - 1) / tn) * (xn - x)
    if (max(abs(xn - x)) < tol) return(xn)
    x <- xn
    t <- tn
  }
  x
}

# Sequential mean over a list of K equal-length vectors (matches numpy's
# axis-0 reduce order: row0 + row1 + ... then divide by K).
list_mean <- function(lst) {
  s <- lst[[1]]
  K <- length(lst)
  if (K > 1) for (k in 2:K) s <- s + lst[[k]]
  s / K
}

# Block-coordinate descent over outcomes for a fixed interior lambda.
# grams: list of K lists (A, b, L); warm: list of K weight vectors.
bcd <- function(grams, lam, warm, n_outer = 300L, tol = 1e-11) {
  K <- length(grams)
  gam <- warm
  gbar <- list_mean(gam)
  for (it in seq_len(n_outer)) {
    move <- 0
    for (k in seq_len(K)) {
      g <- grams[[k]]
      gn <- fista(g$A, g$b, lam, gbar, gam[[k]], g$L, n_iter = 800L, tol = 1e-9)
      move <- max(move, max(abs(gn - gam[[k]])))
      gam[[k]] <- gn
    }
    gbar <- list_mean(gam)
    if (move < tol) break
  }
  gam
}

# Gram matrices (A = X'X/T, b = X'y/T, L = largest eigenvalue of A).
grams_of <- function(ys, Xs) {
  out <- vector("list", length(ys))
  for (k in seq_along(ys)) {
    X <- Xs[[k]]
    y <- ys[[k]]
    Tn <- nrow(X)
    A <- (t(X) %*% X) / Tn
    b <- drop(t(X) %*% y) / Tn
    L <- eigen(A, symmetric = TRUE, only.values = TRUE)$values[1]
    out[[k]] <- list(A = A, b = b, L = L)
  }
  out
}

# Fit every expert on the grid, with EXACT corner solutions at 0 and Inf.
# Returns an M x K x N array. engine = "cpp" runs the compiled solver
# (src/solvers.cpp); engine = "r" runs the pure-R solver below. Gram
# matrices and their largest eigenvalues are computed in R for both engines.
fit_all_lams <- function(ys, Xs, lam_grid, engine = "cpp") {
  K <- length(ys)
  N <- ncol(Xs[[1]])
  M <- length(lam_grid)
  grams <- grams_of(ys, Xs)
  Acat <- list_mean(lapply(grams, `[[`, "A"))
  bcat <- list_mean(lapply(grams, `[[`, "b"))
  Lcat <- eigen(Acat, symmetric = TRUE, only.values = TRUE)$values[1]
  if (engine == "cpp") {
    return(fit_all_lams_cpp(lapply(grams, `[[`, "A"), lapply(grams, `[[`, "b"),
                            vapply(grams, `[[`, 0, "L"), Acat, bcat, Lcat,
                            as.numeric(lam_grid)))
  }
  if (engine != "r") stop(sprintf("unknown engine '%s'", engine))
  unif <- rep(1 / N, N)
  zero <- rep(0, N)
  sep <- lapply(grams, function(g)
    fista(g$A, g$b, 0, zero, unif, g$L, n_iter = 8000L, tol = 1e-12))
  gcat <- fista(Acat, bcat, 0, zero, unif, Lcat, n_iter = 8000L, tol = 1e-12)
  out <- array(0, dim = c(M, K, N))
  warm <- sep
  for (j in seq_len(M)) {
    lam <- lam_grid[j]
    if (lam == 0) {
      for (k in seq_len(K)) out[j, k, ] <- sep[[k]]
    } else if (is.infinite(lam)) {
      for (k in seq_len(K)) out[j, k, ] <- gcat
    } else {
      warm <- bcd(grams, lam, warm)
      for (k in seq_len(K)) out[j, k, ] <- warm[[k]]
    }
  }
  out
}

# Population standard deviation (numpy's default ddof = 0).
pop_sd <- function(x) sqrt(mean((x - mean(x))^2))

# lapply over independent work items on `cores` cores. Fork-based on Unix,
# socket cluster on Windows; cores <= 1 is a plain lapply. Every item is
# computed by the same deterministic code whichever path is taken, so the
# results do not depend on `cores`.
pp_lapply <- function(X, FUN, cores = 1L) {
  cores <- as.integer(cores)
  if (is.na(cores) || cores <= 1L || length(X) <= 1L) return(lapply(X, FUN))
  cores <- min(cores, length(X))
  if (.Platform$OS.type == "unix") {
    return(parallel::mclapply(X, FUN, mc.cores = cores))
  }
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::parLapply(cl, X, FUN)
}

# Weight dispersion V(G) = sum_k ||g_k - gbar||^2 for a K x N matrix.
disp_of <- function(G) {
  G <- as.matrix(G)
  ctr <- colMeans(G)
  sum(sweep(G, 2, ctr)^2)
}

# ---------------------------------------------------------------------------
# main estimator
# ---------------------------------------------------------------------------

#' Fit the partial-pooling synthetic control (matrix interface)
#'
#' Direct port of \code{fit_ppsynth()} in \code{ppsynth_core.py}. Most users
#' will call \code{\link{ppsynth}} on a long panel instead; this function is
#' the engine underneath it and takes the data as matrices.
#'
#' @param Y_treated list of K numeric vectors of length T (the treated unit,
#'   one vector per outcome).
#' @param Y_donors list of K numeric matrices, each T x N0 (the donors, one
#'   matrix per outcome, columns in the same donor order for every outcome).
#' @param T0 number of pre-treatment periods (at least 4).
#' @param method \code{"aggregate"} (default; leave-one-out exponentially
#'   weighted aggregation over the lambda grid) or \code{"cv"}
#'   (leave-future-out selection of a single lambda).
#' @param lam optional fixed lambda; overrides \code{method}. Pass \code{Inf}
#'   for the common-weights corner. Snapped to the nearest grid point.
#' @param lam_grid the lambda grid; default \code{\link{LAM_GRID_DEFAULT}}.
#' @param holdout number of holdout periods for \code{method = "cv"};
#'   default \code{max(2, T0 \%/\% 4)}.
#' @param outcome_names optional outcome labels, used only in the caution
#'   messages.
#' @param engine \code{"cpp"} (default; compiled inner solver) or \code{"r"}
#'   (pure-R reference solver). Same algorithm and settings; see the
#'   certification report shipped in \code{inst/certification}.
#' @param cores number of cores for the leave-one-out folds of the
#'   aggregation step (default 1). Results do not depend on \code{cores}.
#'
#' @details
#' Fittability guardrail: any outcome whose pre-RMSPE exceeds
#' \code{FITTABILITY_THRESHOLD} (0.9) times its own pre-period sd is listed in
#' \code{poor_fit_outcomes}. Aggregation is not corrected for this (its LOO
#' scoring already down-weights a bad expert), but under \code{method = "cv"}
#' the \code{cv_warning} field is populated, since CV's summed
#' leave-future-out loss has no such damping.
#'
#' Support check: for each outcome, each pre-treatment period where the
#' treated value falls outside the donor min-max band is a period no convex
#' combination of donors can match in levels. The per-outcome share of such
#' periods is returned in \code{support_outside_share} (period indices in
#' \code{support_outside_periods}, 1-based), and any violation at all
#' populates \code{support_warning}. Report-only; the estimator is never
#' modified.
#'
#' @return A list with elements \code{gamma} (N0 x K weights),
#'   \code{gamma_common} (N0; cross-outcome mean weights), \code{tau} (K mean
#'   post-treatment effects, original scale), \code{att_path} (list of K
#'   gap paths of length T), \code{method}, \code{lam} (effective lambda),
#'   \code{lam_grid}, \code{agg_mass} (aggregation mass over the grid, or
#'   \code{NULL}), \code{eta}, \code{cv_curve} (or \code{NULL}),
#'   \code{gamma_path} (M x K x N0 array of every expert's weights),
#'   \code{pre_rmspe}, \code{poor_fit_outcomes}, \code{cv_warning},
#'   \code{support_outside_share}, \code{support_outside_periods},
#'   \code{support_warning}, \code{realized_pooling}, \code{pooling_grid},
#'   \code{pooling_q25}, \code{pooling_q75}, \code{pooling_degenerate}.
#' @examples
#' # the shipped test panel as matrices: unit 1 treated from 2010 (T0 = 10)
#' data(ppsynth_test_panel)
#' outs <- c("gdp", "health", "invest")
#' Y_tr <- lapply(outs, function(v) ppsynth_test_panel[[v]][ppsynth_test_panel$id == 1])
#' Y_dn <- lapply(outs, function(v)
#'   sapply(2:9, function(u) ppsynth_test_panel[[v]][ppsynth_test_panel$id == u]))
#' fit <- fit_ppsynth(Y_tr, Y_dn, T0 = 10, outcome_names = outs)
#' fit$agg_mass
#' fit$realized_pooling
#' @seealso \code{\link{ppsynth}}, \code{\link{placebo_inference}}
#' @export
fit_ppsynth <- function(Y_treated, Y_donors, T0, method = "aggregate",
                        lam = NULL, lam_grid = NULL, holdout = NULL,
                        outcome_names = NULL, engine = c("cpp", "r"),
                        cores = 1L) {
  engine <- match.arg(engine)
  K <- length(Y_treated)
  if (is.null(lam_grid)) lam_grid <- LAM_GRID_DEFAULT
  lam_grid <- as.numeric(lam_grid)
  if (is.null(holdout)) holdout <- max(2L, T0 %/% 4L)
  if (T0 < 4) stop("need at least 4 pre-treatment periods")
  Y_donors <- lapply(Y_donors, as.matrix)

  pre <- seq_len(T0)
  mus <- vapply(seq_len(K), function(k) mean(Y_treated[[k]][pre]), 0)
  sds <- vapply(seq_len(K), function(k) pop_sd(Y_treated[[k]][pre]) + 1e-12, 0)
  ys_pre <- lapply(seq_len(K), function(k) (Y_treated[[k]][pre] - mus[k]) / sds[k])
  Xs_pre <- lapply(seq_len(K), function(k) (Y_donors[[k]][pre, , drop = FALSE] - mus[k]) / sds[k])
  N0 <- ncol(Xs_pre[[1]])
  M <- length(lam_grid)

  fits <- fit_all_lams(ys_pre, Xs_pre, lam_grid, engine)   # M x K x N0
  agg_mass <- NULL; eta <- NULL; cv_curve <- NULL

  if (!is.null(lam)) {
    g_fin <- ifelse(is.infinite(lam_grid), 1e18, lam_grid)
    l_fin <- if (is.infinite(lam)) 1e18 else lam
    j <- which.min(abs(g_fin - l_fin))
    gamma_mk <- matrix(fits[j, , ], nrow = K)
    lam_eff <- lam_grid[j]
    method_used <- method
  } else if (method == "aggregate") {
    # leave-one-out losses for every expert (folds are independent)
    loo_fold <- function(t) {
      keep <- pre != t
      f_t <- fit_all_lams(lapply(ys_pre, function(y) y[keep]),
                          lapply(Xs_pre, function(X) X[keep, , drop = FALSE]),
                          lam_grid, engine)
      col <- numeric(M)
      for (j in seq_len(M)) {
        s <- 0
        for (k in seq_len(K)) {
          e <- ys_pre[[k]][t] - sum(Xs_pre[[k]][t, ] * f_t[j, k, ])
          s <- s + e * e
        }
        col[j] <- s
      }
      col
    }
    losses <- do.call(cbind, pp_lapply(seq_len(T0), loo_fold, cores))
    Lbar <- mean(losses)
    eta <- 1 / (2 * Lbar + 1e-12)                      # disclosed rule
    cum <- rowSums(losses)
    agg_mass <- exp(-eta * (cum - min(cum)))
    agg_mass <- agg_mass / sum(agg_mass)
    gamma_mk <- matrix(0, K, N0)
    for (m in seq_len(M)) gamma_mk <- gamma_mk + agg_mass[m] * matrix(fits[m, , ], nrow = K)
    # effective lambda: mass-weighted median over the grid
    finite <- ifelse(is.infinite(lam_grid), 1e6, lam_grid)
    ord <- order(finite)
    cmass <- cumsum(agg_mass[ord])
    lam_eff <- lam_grid[ord][which(cmass >= 0.5)[1]]
    method_used <- method
  } else if (method == "cv") {
    fit_T <- T0 - holdout
    if (fit_T < 3) stop("too few pre-periods for CV; use method='aggregate'")
    f_cv <- fit_all_lams(lapply(ys_pre, function(y) y[seq_len(fit_T)]),
                         lapply(Xs_pre, function(X) X[seq_len(fit_T), , drop = FALSE]),
                         lam_grid, engine)
    hold <- (fit_T + 1):T0
    cv_curve <- vapply(seq_len(M), function(j) {
      s <- 0
      for (k in seq_len(K)) {
        r <- ys_pre[[k]][hold] - drop(Xs_pre[[k]][hold, , drop = FALSE] %*% f_cv[j, k, ])
        s <- s + mean(r^2)
      }
      s
    }, 0)
    j <- which.min(cv_curve)
    gamma_mk <- matrix(fits[j, , ], nrow = K)
    lam_eff <- lam_grid[j]
    method_used <- method
  } else {
    stop(sprintf("unknown method '%s'", method))
  }

  # realized weight pooling (2.6.0): percent reduction in weight dispersion
  # relative to separate fitting, 100*(1 - sqrt(V(G)/V(G_separate))).
  # Quartiles of predictive mass use a first-crossing convention on the grid
  # atoms, no interpolation.
  V0 <- disp_of(matrix(fits[1, , ], nrow = K))
  if (!is.finite(V0) || V0 < 1e-12) {
    pooling_degenerate <- TRUE
    pooling_grid <- NULL; realized_pooling <- NULL
    pooling_q25 <- NULL; pooling_q75 <- NULL
  } else {
    pooling_degenerate <- FALSE
    pooling_grid <- vapply(seq_len(M), function(m)
      100 * (1 - sqrt(max(disp_of(matrix(fits[m, , ], nrow = K)), 0) / V0)), 0)
    realized_pooling <- 100 * (1 - sqrt(max(disp_of(gamma_mk), 0) / V0))
    if (is.null(agg_mass)) {
      pooling_q25 <- NULL; pooling_q75 <- NULL
    } else {
      cm <- cumsum(agg_mass)
      pooling_q25 <- pooling_grid[which(cm >= 0.25)[1]]
      pooling_q75 <- pooling_grid[which(cm >= 0.75)[1]]
    }
  }

  gamma <- t(gamma_mk)                                # N0 x K
  gamma_common <- rowMeans(gamma)

  Tn <- length(Y_treated[[1]])
  post <- (T0 + 1):Tn
  tau <- numeric(K); att_path <- vector("list", K); pre_rmspe <- numeric(K)
  for (k in seq_len(K)) {
    synth <- drop(Y_donors[[k]] %*% gamma[, k])       # weights sum to 1 ->
    gap <- Y_treated[[k]] - synth                      # original scale is exact
    att_path[[k]] <- gap
    tau[k] <- mean(gap[post])
    pre_rmspe[k] <- sqrt(mean(gap[pre]^2))
  }

  # fittability guardrail
  poor_fit <- pre_rmspe > (FITTABILITY_THRESHOLD * sds)
  poor_fit_idx <- which(poor_fit)
  cv_warning <- NULL
  if (method == "cv" && length(poor_fit_idx) > 0) {
    names_ <- if (!is.null(outcome_names)) outcome_names else as.character(seq_len(K) - 1L)
    bad <- paste(names_[poor_fit_idx], collapse = ", ")
    one <- length(poor_fit_idx) == 1
    cv_warning <- sprintf(
      "CAUTION: outcome%s %s %s poorly (pre-RMSPE exceeds %gx %s own pre-period sd). Consider removing %s from the model.",
      if (one) "" else "s", bad, if (one) "fits" else "fit",
      FITTABILITY_THRESHOLD, if (one) "its" else "their", if (one) "it" else "them")
  }

  # support check (per-period convex hull; boundary counts as inside;
  # fires on any violation; report-only)
  support_outside_periods <- vector("list", K)
  support_outside_count <- integer(K)
  for (k in seq_len(K)) {
    D <- Y_donors[[k]][pre, , drop = FALSE]
    lo <- apply(D, 1, min)
    hi <- apply(D, 1, max)
    yt <- Y_treated[[k]][pre]
    out_t <- which(yt < lo | yt > hi)
    support_outside_periods[[k]] <- out_t
    support_outside_count[k] <- length(out_t)
  }
  support_outside_share <- support_outside_count / T0
  support_warning <- NULL
  if (any(support_outside_count > 0)) {
    names_s <- if (!is.null(outcome_names)) outcome_names else as.character(seq_len(K) - 1L)
    bad_s <- which(support_outside_count > 0)
    if (length(bad_s) == 1) {
      k <- bad_s[1]
      support_warning <- sprintf(
        "CAUTION: outcome '%s' lies outside the support of the donor outcomes in %d of %d pre-treatment periods (%.0f%%). Synthetic control weights cannot extrapolate beyond the donor range. Check the data and consider revising the donor pool before proceeding.",
        names_s[k], support_outside_count[k], T0, 100 * support_outside_share[k])
    } else {
      parts <- paste(sprintf("'%s' (%d of %d, %.0f%%)", names_s[bad_s],
                             support_outside_count[bad_s], T0,
                             100 * support_outside_share[bad_s]), collapse = ", ")
      support_warning <- sprintf(
        "CAUTION: outcomes %s lie outside the support of the donor outcomes in the pre-treatment period. Synthetic control weights cannot extrapolate beyond the donor range. Check the data and consider revising the donor pool before proceeding.",
        parts)
    }
  }

  list(gamma = gamma, gamma_common = gamma_common, tau = tau,
       att_path = att_path, method = method_used, lam = lam_eff,
       lam_grid = lam_grid, agg_mass = agg_mass, eta = eta,
       cv_curve = cv_curve, gamma_path = fits, pre_rmspe = pre_rmspe,
       poor_fit_outcomes = poor_fit_idx, cv_warning = cv_warning,
       support_outside_share = support_outside_share,
       support_outside_periods = support_outside_periods,
       support_warning = support_warning,
       realized_pooling = realized_pooling,
       pooling_grid = pooling_grid, pooling_q25 = pooling_q25,
       pooling_q75 = pooling_q75,
       pooling_degenerate = pooling_degenerate,
       T0 = T0, N0 = N0, K = K, engine = engine)
}

# ---------------------------------------------------------------------------
# placebo (in-space permutation) inference
# ---------------------------------------------------------------------------

rmspe_stats <- function(fit, T0, eps = 1e-12) {
  K <- length(fit$att_path)
  Tn <- length(fit$att_path[[1]])
  ratios <- numeric(K)
  for (k in seq_len(K)) {
    gap <- fit$att_path[[k]]
    pre_ <- sqrt(mean(gap[seq_len(T0)]^2)) + eps
    post_ <- sqrt(mean(gap[(T0 + 1):Tn]^2))
    ratios[k] <- post_ / pre_
  }
  ratios
}

#' Whole-pipeline placebo inference (matrix interface)
#'
#' Reruns the ENTIRE chosen pipeline (leave-one-out aggregation or CV
#' included) for each donor treated in turn as a placebo, with the actual
#' treated unit returned to the donor pool. Direct port of
#' \code{placebo_inference()} in \code{ppsynth_core.py}.
#'
#' @inheritParams fit_ppsynth
#' @param convention \code{"phipson_smyth"} (default;
#'   p = (1 + #\{>=\})/(N0+1), always valid, minimum 1/(N0+1)) or
#'   \code{"strict"} (p = #\{>\}/(N0+1); synth-style strict exceedance, can be
#'   zero, anti-conservative).
#' @param verbose print a progress line after each placebo fit (ignored when
#'   \code{cores > 1}).
#' @param cores number of cores: the placebo units are fitted in parallel
#'   (each with a single-core leave-one-out step). Results do not depend on
#'   \code{cores}.
#'
#' @details Per-outcome p-values rank the treated unit's post/pre RMSPE ratio
#' against the placebo ratios. Joint p-values use the max and mean of the
#' per-outcome ratios studentized over ALL N0+1 units (treated and donors),
#' which preserves exchangeability under the sharp null.
#'
#' @return A list: \code{p_outcome} (K), \code{p_joint_max},
#'   \code{p_joint_mean}, \code{stats_treated} (K), \code{stats_placebo}
#'   (N0 x K), \code{treated_fit} (the \code{\link{fit_ppsynth}} result for
#'   the treated unit), \code{placebo_gaps} (list over donors of K-lists of
#'   gap paths), \code{convention}.
#' @examples
#' data(ppsynth_test_panel)
#' outs <- c("gdp", "health", "invest")
#' Y_tr <- lapply(outs, function(v) ppsynth_test_panel[[v]][ppsynth_test_panel$id == 1])
#' Y_dn <- lapply(outs, function(v)
#'   sapply(2:9, function(u) ppsynth_test_panel[[v]][ppsynth_test_panel$id == u]))
#' inf <- placebo_inference(Y_tr, Y_dn, T0 = 10, outcome_names = outs)
#' inf$p_outcome
#' inf$p_joint_max
#' @seealso \code{\link{ppsynth}}, \code{\link{fit_ppsynth}}
#' @export
placebo_inference <- function(Y_treated, Y_donors, T0, method = "aggregate",
                              lam = NULL, lam_grid = NULL, holdout = NULL,
                              convention = "phipson_smyth",
                              outcome_names = NULL, verbose = FALSE,
                              engine = c("cpp", "r"), cores = 1L) {
  engine <- match.arg(engine)
  K <- length(Y_treated)
  Y_donors <- lapply(Y_donors, as.matrix)
  N0 <- ncol(Y_donors[[1]])
  cores <- as.integer(cores)

  fit_unit <- function(y_tr, y_dn, fold_cores)
    fit_ppsynth(y_tr, y_dn, T0, method = method, lam = lam,
                lam_grid = lam_grid, holdout = holdout,
                outcome_names = outcome_names, engine = engine,
                cores = fold_cores)

  # the treated fit may use every core on its folds; the placebo units are
  # then spread across cores, one core each
  treated_fit <- fit_unit(Y_treated, Y_donors, cores)
  stats_tr <- rmspe_stats(treated_fit, T0)

  placebo_one <- function(j) {
    y_tr_j <- lapply(seq_len(K), function(k) Y_donors[[k]][, j])
    y_dn_j <- lapply(seq_len(K), function(k)
      cbind(Y_donors[[k]][, -j, drop = FALSE], Y_treated[[k]]))
    fit_j <- fit_unit(y_tr_j, y_dn_j, 1L)
    if (verbose && cores <= 1L) message(sprintf("placebo %d/%d done", j, N0))
    list(stats = rmspe_stats(fit_j, T0), gaps = fit_j$att_path)
  }
  pl <- pp_lapply(seq_len(N0), placebo_one, cores)
  stats_pl <- do.call(rbind, lapply(pl, `[[`, "stats"))
  placebo_gaps <- lapply(pl, `[[`, "gaps")

  if (convention == "strict") {
    p_outcome <- vapply(seq_len(K), function(k)
      sum(stats_pl[, k] > stats_tr[k]) / (N0 + 1), 0)
  } else {
    p_outcome <- vapply(seq_len(K), function(k)
      (1 + sum(stats_pl[, k] >= stats_tr[k])) / (N0 + 1), 0)
  }

  allstats <- rbind(stats_pl, stats_tr)
  mu <- colMeans(allstats)
  sd_ <- apply(allstats, 2, pop_sd) + 1e-12
  z_tr <- (stats_tr - mu) / sd_
  z_pl <- sweep(sweep(stats_pl, 2, mu), 2, sd_, "/")
  zmax_pl <- apply(z_pl, 1, max)
  zmean_pl <- rowMeans(z_pl)
  if (convention == "strict") {
    p_joint_max <- sum(zmax_pl > max(z_tr)) / (N0 + 1)
    p_joint_mean <- sum(zmean_pl > mean(z_tr)) / (N0 + 1)
  } else {
    p_joint_max <- (1 + sum(zmax_pl >= max(z_tr))) / (N0 + 1)
    p_joint_mean <- (1 + sum(zmean_pl >= mean(z_tr))) / (N0 + 1)
  }

  list(p_outcome = p_outcome, p_joint_max = p_joint_max,
       p_joint_mean = p_joint_mean, stats_treated = stats_tr,
       stats_placebo = stats_pl, treated_fit = treated_fit,
       placebo_gaps = placebo_gaps, convention = convention)
}
