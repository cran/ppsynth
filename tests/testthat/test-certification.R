# Certification against the Python reference engine (ppsynth_core.py /
# ppsynth_fast.py, PPSYNTH v2.6.1). The compiled engine runs every fixture
# on every check (seconds). The pure-R engine and the Germany placebo
# inference run only with PPSYNTH_FULL_GATE=1 (minutes to an hour);
# tools/gate.R runs the same comparisons and records the deviations in
# inst/certification/results/.

full_gate <- identical(Sys.getenv("PPSYNTH_FULL_GATE"), "1")
engines <- if (full_gate) c("cpp", "r") else "cpp"

fit_spec <- function(name) {
  if (name %in% c("aggregate", "aggregate_fast", "main", "boundary", "single")) return(list(method = "aggregate"))
  if (name == "cv") return(list(method = "cv"))
  if (name == "cv_holdout3") return(list(method = "cv", holdout = 3))
  if (startsWith(name, "lam_")) { v <- sub("lam_", "", name); return(list(lam = if (v == "inf") Inf else as.numeric(v))) }
  stop(name)
}
inf_spec <- function(name) {
  if (name == "lam_1_ps") return(list(lam = 1, convention = "phipson_smyth"))
  list(method = sub("_(ps|strict)$", "", name),
       convention = if (grepl("_strict$", name)) "strict" else "phipson_smyth")
}
expect_inf_matches <- function(ri, pi, tol = CERT_TOL) {
  expect_equal(ri$p_outcome, fx_num(pi$p_outcome), ignore_attr = TRUE)
  expect_equal(ri$p_joint_max, fx_num(pi$p_joint_max))
  expect_equal(ri$p_joint_mean, fx_num(pi$p_joint_mean))
  expect_equal(ri$stats_treated, fx_num(pi$stats_treated), tolerance = tol, ignore_attr = TRUE)
  expect_equal(ri$stats_placebo, fx_arr(pi$stats_placebo), tolerance = tol, ignore_attr = TRUE)
  expect_equal(unlist(ri$placebo_gaps), unlist(lapply(pi$placebo_gaps, function(pg) lapply(pg, fx_num))),
               tolerance = tol, ignore_attr = TRUE)
  expect_fit_matches(ri$treated_fit, pi$treated_fit, tol)
}

for (eng in engines) {

  test_that(paste("toy self-test design matches the Python engine, engine =", eng), {
    fx <- read_fixture("selftest_toy")
    d <- fx$data
    rf <- fit_ppsynth(d$Y_treated, d$Y_donors, T0 = d$T0, method = "aggregate", engine = eng)
    expect_fit_matches(rf, fx$fits$aggregate)
    K <- length(d$Y_treated); T0 <- d$T0
    ys <- lapply(1:K, function(k) { y <- d$Y_treated[[k]][1:T0]; (y - mean(y)) / (ppsynth:::pop_sd(y) + 1e-12) })
    Xs <- lapply(1:K, function(k) { y <- d$Y_treated[[k]][1:T0]; (d$Y_donors[[k]][1:T0, ] - mean(y)) / (ppsynth:::pop_sd(y) + 1e-12) })
    mf <- ppsynth:::fit_all_lams(ys, Xs, fx_num(fx$monotone_grid), eng)
    expect_equal(mf, fx_arr(fx$monotone_fits), tolerance = CERT_TOL, ignore_attr = TRUE)
  })

  test_that(paste("guardrail design matches the Python engine, engine =", eng), {
    fx <- read_fixture("selftest_guardrail")
    d <- fx$data
    rf <- fit_ppsynth(d$Y_treated, d$Y_donors, T0 = d$T0, method = "aggregate", engine = eng)
    expect_fit_matches(rf, fx$fits$aggregate)
    rf <- fit_ppsynth(d$Y_treated, d$Y_donors, T0 = d$T0, method = "cv",
                      outcome_names = c("clean", "noisy"), engine = eng)
    expect_fit_matches(rf, fx$fits$cv)
    expect_match(rf$cv_warning, "noisy")
  })

  test_that(paste("support-check designs match the Python engine, engine =", eng), {
    fx <- read_fixture("selftest_support")
    d <- fx$data
    rf <- fit_ppsynth(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = c("inside", "poking"), engine = eng)
    expect_fit_matches(rf, fx$fits$main)
    rf <- fit_ppsynth(list(fx_num(fx$X0), d$Y_treated[[1]]), list(d$Y_donors[[1]], d$Y_donors[[1]]),
                      T0 = d$T0, engine = eng)
    expect_fit_matches(rf, fx$fits$boundary)
    rf <- fit_ppsynth(list(fx_num(fx$Yone)), list(d$Y_donors[[1]]), T0 = d$T0, outcome_names = "once", engine = eng)
    expect_fit_matches(rf, fx$fits$single)
  })

  test_that(paste("Germany panel fits match the Python engine, engine =", eng), {
    if (eng == "r") skip_if_not(full_gate, "pure-R Germany fits take minutes; PPSYNTH_FULL_GATE=1")
    fx <- read_fixture("germany_fits")
    d <- fx$data
    for (nm in names(fx$fits)) {
      rf <- do.call(fit_ppsynth, c(list(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names,
                                        engine = eng), fit_spec(nm)))
      expect_fit_matches(rf, fx$fits[[nm]])
    }
  })

  test_that(paste("test panel fits and placebo inference match the Python engine, engine =", eng), {
    if (eng == "r") skip_if_not(full_gate, "pure-R placebo inference takes an hour; PPSYNTH_FULL_GATE=1")
    fx <- read_fixture("testpanel")
    d <- fx$data
    for (nm in names(fx$fits)) {
      rf <- do.call(fit_ppsynth, c(list(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names,
                                        engine = eng), fit_spec(nm)))
      expect_fit_matches(rf, fx$fits[[nm]])
    }
    for (nm in names(fx$inference)) {
      ri <- do.call(placebo_inference, c(list(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names,
                                              engine = eng), inf_spec(nm)))
      expect_inf_matches(ri, fx$inference[[nm]])
    }
  })
}

test_that("Germany placebo inference matches the numba engine (PPSYNTH_FULL_GATE=1)", {
  skip_if_not(full_gate, "set PPSYNTH_FULL_GATE=1 to run the Germany placebo gate")
  fx <- read_fixture("germany_inference")
  d <- fx$data
  rf <- fit_ppsynth(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names)
  expect_fit_matches(rf, fx$fits$aggregate_fast)
  for (nm in names(fx$inference)) {
    ri <- do.call(placebo_inference, c(list(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names,
                                            cores = 2L), inf_spec(nm)))
    expect_inf_matches(ri, fx$inference[[nm]])
  }
})

test_that("results do not depend on cores", {
  fx <- read_fixture("testpanel")
  d <- fx$data
  f1 <- fit_ppsynth(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names, cores = 1L)
  f2 <- fit_ppsynth(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names, cores = 2L)
  expect_identical(f1$gamma, f2$gamma)
  expect_identical(f1$agg_mass, f2$agg_mass)
  expect_identical(f1$gamma_path, f2$gamma_path)
  i1 <- placebo_inference(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names, lam = 1, cores = 1L)
  i2 <- placebo_inference(d$Y_treated, d$Y_donors, T0 = d$T0, outcome_names = d$names, lam = 1, cores = 2L)
  expect_identical(i1$stats_placebo, i2$stats_placebo)
  expect_identical(i1$p_outcome, i2$p_outcome)
})
