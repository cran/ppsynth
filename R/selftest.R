# selftest.R -- R counterparts of the three Python self-tests that
# ppsynth_setup runs (_ppsynth_selftest, _ppsynth_selftest_guardrail,
# _ppsynth_selftest_support). These check invariants, not exact floats; the
# exact-float certification against the Python engine lives in the package
# tests (tests/testthat/test-certification.R) and in tools/gate.R.

# Fixed pseudo-random normals for the self-test designs: a 31-bit linear
# congruential sequence (Park-Miller multiplier) started at `stream`, pushed
# through Box-Muller. Pure arithmetic on doubles, so it touches neither
# set.seed() nor .Random.seed and gives the same numbers on every platform.
det_normal <- function(n, stream) {
  m <- 2147483647
  x <- as.numeric(stream)
  for (i in 1:20) x <- (16807 * x) %% m      # burn-in past the small start value
  u <- numeric(2 * n)
  for (i in seq_len(2 * n)) {
    x <- (16807 * x) %% m
    u[i] <- x / m
  }
  u1 <- u[seq(1, 2 * n, by = 2)]
  u2 <- u[seq(2, 2 * n, by = 2)]
  sqrt(-2 * log(u1)) * cos(2 * pi * u2)
}

#' Run the ppsynth self-tests
#'
#' Deterministic toy fits that check the estimator's invariants (simplex
#' weights, sane aggregation mass, finite effects, in-sample loss monotone
#' in lambda), the fittability guardrail (aggregation is not flagged;
#' CV carries the caution) and the donor-support check (fires on any
#' violation, counts periods correctly, boundary counts as inside).
#'
#' The toy designs are built from a fixed internal number sequence, so the
#' function neither sets nor disturbs R's random-number generator.
#'
#' @param engine \code{"cpp"} (default) or \code{"r"}; the engine to test.
#' @return \code{TRUE} invisibly; stops with a message on the first failure.
#' @examples
#' ppsynth_selftest()
#' @export
ppsynth_selftest <- function(engine = c("cpp", "r")) {
  engine <- match.arg(engine)
  # 1. toy fit, invariants
  T0 <- 10; T1 <- 4; J <- 5; K <- 2
  X <- lapply(seq_len(K), function(k) apply(matrix(det_normal((T0 + T1) * J, 10L + k), T0 + T1, J), 2, cumsum))
  w <- c(0.4, 0.3, 0.2, 0.1, 0.0)
  Y <- lapply(seq_along(X), function(k) drop(X[[k]] %*% w) + 0.1 * det_normal(T0 + T1, 20L + k))
  fit <- fit_ppsynth(Y, X, T0 = T0, method = "aggregate", engine = engine)
  g <- fit$gamma
  if (abs(sum(g) - K) >= 1e-6) stop("weights do not sum to one per outcome")
  if (min(g) < -1e-9) stop("negative weight")
  m <- fit$agg_mass
  if (abs(sum(m) - 1) >= 1e-8 || min(m) < 0) stop("bad mass")
  if (!all(is.finite(fit$tau))) stop("non-finite effect")
  # monotonicity regression test for the corner-convergence bug class
  ys <- lapply(seq_len(K), function(k) { y <- Y[[k]][1:T0]; (y - mean(y)) / (pop_sd(y) + 1e-12) })
  Xs <- lapply(seq_len(K), function(k) { y <- Y[[k]][1:T0]; (X[[k]][1:T0, ] - mean(y)) / (pop_sd(y) + 1e-12) })
  fits <- fit_all_lams(ys, Xs, c(0, 1e-4, 0.1), engine)
  prev <- -Inf
  for (j in 1:3) {
    G <- matrix(fits[j, , ], nrow = K)
    ins <- sum(vapply(seq_len(K), function(k) sum((ys[[k]] - drop(Xs[[k]] %*% G[k, ]))^2), 0))
    if (ins < prev - 1e-6) stop("in-sample SSE not monotone in lambda")
    prev <- ins
  }

  # 2. guardrail: one outcome pure noise, one clean
  T0 <- 14; T1 <- 6; J <- 5
  Xclean <- apply(matrix(det_normal((T0 + T1) * J, 31L), T0 + T1, J), 2, cumsum)
  w <- c(0.5, 0.3, 0.2, 0, 0)
  Yclean <- drop(Xclean %*% w) + 0.05 * det_normal(T0 + T1, 32L)
  Xnoise <- matrix(det_normal((T0 + T1) * J, 33L), T0 + T1, J)
  Ynoise <- 5 * det_normal(T0 + T1, 34L)
  Y <- list(Yclean, Ynoise); X <- list(Xclean, Xnoise)
  fit_agg <- fit_ppsynth(Y, X, T0 = T0, method = "aggregate", engine = engine)
  if (!identical(as.integer(fit_agg$poor_fit_outcomes), 2L)) stop("noise outcome should flag")
  if (!is.null(fit_agg$cv_warning)) stop("aggregation must not carry a CV warning")
  fit_cv <- fit_ppsynth(Y, X, T0 = T0, method = "cv", outcome_names = c("clean", "noisy"), engine = engine)
  if (!identical(as.integer(fit_cv$poor_fit_outcomes), 2L)) stop("noise outcome should flag under cv")
  if (is.null(fit_cv$cv_warning)) stop("CV must be annotated when hijacked")
  if (!grepl("noisy", fit_cv$cv_warning)) stop("CV warning must name the outcome")

  # 3. support check
  T0 <- 20; T1 <- 5; J <- 6
  X <- apply(matrix(det_normal((T0 + T1) * J, 41L), T0 + T1, J), 2, cumsum)
  Yin <- drop(X %*% rep(1 / J, J))
  Yout <- Yin
  push <- c(3, 8, 12, 17)                      # 1-based
  Yout[push] <- apply(X[push, , drop = FALSE], 1, max) + 1
  fit <- fit_ppsynth(list(Yin, Yout), list(X, X), T0 = T0, outcome_names = c("inside", "poking"), engine = engine)
  if (fit$support_outside_share[1] != 0) stop("interior combo flagged")
  if (round(fit$support_outside_share[2] * T0) != 4) stop("wrong violation count")
  if (!identical(as.integer(fit$support_outside_periods[[2]]), as.integer(push))) stop("wrong periods")
  wmsg <- fit$support_warning
  if (is.null(wmsg) || !grepl("poking", wmsg) || grepl("inside", strsplit(wmsg, "outcome '")[[1]][1]))
    stop("warning must name only the offender")
  if (!grepl("4 of 20", wmsg) || !grepl("20%", wmsg)) stop("warning must report count and share")
  fit_b <- fit_ppsynth(list(X[, 1], Yin), list(X, X), T0 = T0, engine = engine)
  if (!is.null(fit_b$support_warning)) stop("band edge must count as inside")
  Yone <- Yin; Yone[6] <- max(X[6, ]) + 1
  fit_1 <- fit_ppsynth(list(Yone), list(X), T0 = T0, outcome_names = "once", engine = engine)
  if (is.null(fit_1$support_warning) || !grepl("1 of 20", fit_1$support_warning))
    stop("must fire on any violation")

  message(sprintf("ppsynth self-tests passed (invariants, guardrail, support check; engine = %s).", engine))
  invisible(TRUE)
}
