## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 5)

## -----------------------------------------------------------------------------
library(ppsynth)
data(germany_panel)
outs <- c("proprights", "money", "inttrade", "regulation", "csh_g")

fit <- ppsynth(germany_panel, outcomes = outs,
               id = "id", time = "year", trunit = 8, trperiod = 1990)
fit

## -----------------------------------------------------------------------------
print(fit, massprofile = TRUE)

## -----------------------------------------------------------------------------
round(fit$gamma, 2)
fit$tau

## -----------------------------------------------------------------------------
plot(fit, which = "paths")   # treated (solid) vs synthetic (dash), one panel per outcome
plot(fit, which = "mass")    # aggregation mass over the lambda grid

## ----eval = FALSE-------------------------------------------------------------
#  inf <- ppsynth(germany_panel, outcomes = outs,
#                 id = "id", time = "year", trunit = 8, trperiod = 1990,
#                 inference = TRUE, cores = 2)
#  print(inf, effects = TRUE)
#  plot(inf, which = "gaps")

## -----------------------------------------------------------------------------
Y_tr <- lapply(outs, function(v) germany_panel[[v]][germany_panel$id == 8])
Y_dn <- lapply(outs, function(v)
  sapply(setdiff(sort(unique(germany_panel$id)), 8),
         function(u) germany_panel[[v]][germany_panel$id == u]))
eng <- fit_ppsynth(Y_tr, Y_dn, T0 = 20)
round(eng$agg_mass, 3)

