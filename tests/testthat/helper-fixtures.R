# helper-fixtures.R -- read the Python reference fixtures shipped in
# inst/certification/fixtures (written by tools/make_fixtures.py from
# ppsynth_core.py, PPSYNTH v2.6.1).

fixture_path <- function(name) {
  p <- system.file("certification", "fixtures", paste0(name, ".json"), package = "ppsynth")
  if (p == "") p <- file.path("..", "..", "inst", "certification", "fixtures", paste0(name, ".json"))
  p
}

fx_num <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) x <- unlist(x)
  x[x == "inf"] <- "Inf"
  as.numeric(x)
}

fx_depth <- function(x) if (!is.list(x)) 0L else 1L + fx_depth(x[[1]])
fx_arr <- function(x) {
  d <- fx_depth(x)
  if (d <= 1) return(fx_num(x))
  if (d == 2) return(do.call(rbind, lapply(x, fx_num)))
  M <- length(x); K <- length(x[[1]]); N <- length(x[[1]][[1]])
  a <- array(0, c(M, K, N))
  for (m in seq_len(M)) for (k in seq_len(K)) a[m, k, ] <- fx_num(x[[m]][[k]])
  a
}

read_fixture <- function(name) {
  skip_if_not_installed("jsonlite")
  p <- fixture_path(name)
  skip_if(!file.exists(p), paste("fixture missing:", name))
  fx <- jsonlite::fromJSON(p, simplifyVector = FALSE)
  fx$data <- list(Y_treated = lapply(fx$Y_treated, fx_num),
                  Y_donors = lapply(fx$Y_donors, fx_arr),
                  T0 = fx$T0, names = unlist(fx$outcome_names))
  fx
}

# Certification tolerance for cross-language floating-point agreement.
# The two engines run the same algorithm with the same settings, but their
# BLAS kernels and reduction orders differ at the last bit; the measured
# maximum deviation over all fixtures is recorded in
# inst/certification/results/ and this constant sits above it.
CERT_TOL <- 1e-8

expect_fit_matches <- function(rf, pf, tol = CERT_TOL) {
  expect_equal(rf$gamma, fx_arr(pf$gamma), tolerance = tol, ignore_attr = TRUE)
  expect_equal(rf$gamma_common, fx_num(pf$gamma_common), tolerance = tol, ignore_attr = TRUE)
  expect_equal(rf$tau, fx_num(pf$tau), tolerance = tol, ignore_attr = TRUE)
  expect_equal(unlist(rf$att_path), unlist(lapply(pf$att_path, fx_num)), tolerance = tol, ignore_attr = TRUE)
  expect_equal(rf$lam, fx_num(pf$lam), ignore_attr = TRUE)
  expect_equal(rf$pre_rmspe, fx_num(pf$pre_rmspe), tolerance = tol, ignore_attr = TRUE)
  expect_equal(rf$gamma_path, fx_arr(pf$gamma_path), tolerance = tol, ignore_attr = TRUE)
  expect_equal(as.numeric(rf$poor_fit_outcomes) - 1, fx_num(pf$poor_fit_outcomes), ignore_attr = TRUE)
  expect_equal(rf$cv_warning, pf$cv_warning)
  expect_equal(rf$support_outside_share, fx_num(pf$support_outside_share), ignore_attr = TRUE)
  expect_equal(unlist(lapply(rf$support_outside_periods, function(p) p - 1)),
               unlist(lapply(pf$support_outside_periods, fx_num)), ignore_attr = TRUE)
  expect_equal(rf$support_warning, pf$support_warning)
  expect_equal(rf$pooling_degenerate, isTRUE(pf$pooling_degenerate))
  if (!is.null(pf$agg_mass)) {
    expect_equal(rf$agg_mass, fx_num(pf$agg_mass), tolerance = tol, ignore_attr = TRUE)
    expect_equal(rf$eta, fx_num(pf$eta), tolerance = tol, ignore_attr = TRUE)
    expect_equal(rf$pooling_q25, fx_num(pf$pooling_q25), tolerance = tol, ignore_attr = TRUE)
    expect_equal(rf$pooling_q75, fx_num(pf$pooling_q75), tolerance = tol, ignore_attr = TRUE)
  } else {
    expect_null(rf$agg_mass)
  }
  if (!is.null(pf$cv_curve)) {
    expect_equal(rf$cv_curve, fx_num(pf$cv_curve), tolerance = tol, ignore_attr = TRUE)
  } else {
    expect_null(rf$cv_curve)
  }
  if (!isTRUE(pf$pooling_degenerate)) {
    expect_equal(rf$realized_pooling, fx_num(pf$realized_pooling), tolerance = tol, ignore_attr = TRUE)
    expect_equal(rf$pooling_grid, fx_num(pf$pooling_grid), tolerance = tol, ignore_attr = TRUE)
  }
}
