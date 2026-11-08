# Front-end tests: panel construction, returned object, print and plot.

test_that("self-tests pass", {
  expect_true(suppressMessages(ppsynth_selftest("cpp")))
  skip_if_not(identical(Sys.getenv("PPSYNTH_FULL_GATE"), "1"), "pure-R self-tests take minutes; PPSYNTH_FULL_GATE=1")
  expect_true(suppressMessages(ppsynth_selftest("r")))
})

test_that("ppsynth() builds the panel like the Stata glue and returns the r() set", {
  data(ppsynth_test_panel, package = "ppsynth")
  fit <- ppsynth(ppsynth_test_panel, outcomes = c("gdp", "health", "invest"),
                 id = "id", time = "year", trunit = 1, trperiod = 2010, lambda = Inf)
  expect_s3_class(fit, "ppsynth")
  expect_equal(fit$T0, 10); expect_equal(fit$T1, 6); expect_equal(fit$N0, 8)
  expect_equal(dim(fit$gamma), c(8, 3))
  expect_equal(rownames(fit$gamma), as.character(2:9))
  expect_equal(colnames(fit$gamma), c("gdp", "health", "invest"))
  # common-weights corner: every outcome carries the same weights
  expect_equal(fit$gamma[, 1], fit$gamma[, 2], ignore_attr = TRUE)
  expect_equal(colSums(fit$gamma), rep(1, 3), ignore_attr = TRUE)
  expect_true(all(fit$gamma >= -1e-9))
  expect_equal(fit$lambda, Inf)
  expect_equal(dim(fit$effects), c(6, 3))
  expect_equal(fit$y_synth + fit$att_path, fit$y_treated)
  expect_null(fit$p_outcome)
  # the matrix engine gives the same answer on the same panel
  ytr <- lapply(c("gdp", "health", "invest"), function(v) ppsynth_test_panel[[v]][ppsynth_test_panel$id == 1])
  ydn <- lapply(c("gdp", "health", "invest"), function(v)
    sapply(2:9, function(u) ppsynth_test_panel[[v]][ppsynth_test_panel$id == u]))
  eng <- fit_ppsynth(ytr, ydn, T0 = 10, lam = Inf)
  expect_equal(unname(fit$gamma), eng$gamma)
  expect_equal(unname(fit$tau), eng$tau)
})

test_that("print and plot run", {
  data(ppsynth_test_panel, package = "ppsynth")
  fit <- ppsynth(ppsynth_test_panel, outcomes = c("gdp", "health", "invest"),
                 id = "id", time = "year", trunit = 1, trperiod = 2010, lambda = 0.3)
  out <- capture.output(print(fit, effects = TRUE))
  expect_true(any(grepl("Partial-pooling synthetic control \\(3 outcomes, method = aggregate\\)", out)))
  expect_true(any(grepl("pre: 10   post: 6   donors: 8", out)))
  expect_true(any(grepl("per-period effects", out)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(fit, which = "paths"))
  expect_message(plot(fit, which = "gaps"), "requires inference")
})

test_that("errors mirror the Stata command", {
  data(ppsynth_test_panel, package = "ppsynth")
  expect_error(ppsynth(ppsynth_test_panel, outcomes = "gdp", id = "id", time = "year",
                       trunit = 99, trperiod = 2010), "treated unit 99 not found")
  bad <- ppsynth_test_panel[-5, ]
  expect_error(ppsynth(bad, outcomes = "gdp", id = "id", time = "year",
                       trunit = 1, trperiod = 2010, lambda = 0), "unbalanced panel")
  expect_error(ppsynth(ppsynth_test_panel, outcomes = "gdp", id = "id", time = "year",
                       trunit = 1, trperiod = 2010, likesynth = TRUE,
                       convention = "phipson_smyth"), "likesynth implies")
})
