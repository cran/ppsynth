#' German reunification panel
#'
#' The worked example of Section 5A of the working paper: five outcomes for
#' 21 OECD countries, 1970-2004. Germany (\code{id = 8}) is treated at
#' reunification (first treated year 1990), giving 20 pre-treatment years and
#' 15 post years. The outcomes are four subcomponents of the historical index
#' of economic liberty (property rights, sound money, freedom to trade,
#' regulation) and the government consumption share. Identical to
#' \code{germany_panel.dta} shipped with the Stata package.
#'
#' @format A data frame with 735 rows and 8 columns: \code{id}, \code{country},
#'   \code{year}, \code{proprights}, \code{money}, \code{inttrade},
#'   \code{regulation}, \code{csh_g}.
#' @source \url{https://github.com/angusgrier/PPSYNTH}
"germany_panel"

#' Small simulated test panel
#'
#' A 9-unit, 16-period simulated panel with three outcomes
#' (\code{gdp}, \code{health}, \code{invest}), shipped with the Stata package
#' as \code{ppsynth_test_panel.dta}. Used by the certification tests with
#' unit 1 treated from 2010.
#'
#' @format A data frame with 144 rows and 5 columns: \code{id}, \code{year},
#'   \code{gdp}, \code{health}, \code{invest}.
#' @source \url{https://github.com/angusgrier/PPSYNTH}
"ppsynth_test_panel"
