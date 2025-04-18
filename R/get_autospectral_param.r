# get_autospill_param.r


#' Get Autospectral Parameters
#'
#' This function retrieves autospectral parameters for a specified cytometer.
#'     It sets up directory parameters for figures and tables if required.
#'
#' @title Get Autospectral Parameters
#' @description Retrieves autospectral parameters for a specified cytometer.
#' @importsFrom base get0
#' @param cytometer The type of cytometer (default is "aurora").
#' @param figures Logical indicating whether to set up directory parameters for
#'     figures and tables (default is TRUE).
#' @return A list of autospectral parameters.
#' @examples
#' \dontrun{
#' get.autospectral.param(cytometer = "aurora", figures = TRUE)
#' }


get.autospectral.param <- function( cytometer = "aurora", figures = TRUE )
{
  autosp.param <- get.autospectral.param.minimal()

  if( figures ){

    # directory parameters

    autosp.param$figure.scatter.dir.base <- "figure_scatter"

    autosp.param$figure.gate.dir <- "figure_gate"
    autosp.param$figure.af.dir <- "figure_autofluorescence"
    autosp.param$figure.peacoqc.dir <- "figure_peacoQC"
    autosp.param$figure.clean.control.dir <- "figure_clean_controls"
    autosp.param$figure.spectral.ribbon.dir <- "figure_spectral_ribbon"
    autosp.param$figure.convergence.dir <- "figure_convergence"
    autosp.param$figure.spectra.dir <- "figure_spectra"
    autosp.param$figure.slope.error.dir <- "figure_slope_error"
    autosp.param$figure.skewness.dir <- "figure_skewness"
    autosp.param$figure.similarity.heatmap.dir <- "figure_similarity_heatmap"

    autosp.param$table.convergence.dir <- "table_convergence"
    autosp.param$table.spectra.dir <- "table_spectra"
    autosp.param$table.slope.error.dir <- "table_slope_error"
    autosp.param$table.skewness.dir <- "table_skewness"

  }

  # cytometer-specific parameters

  get.param.function <- get0( sprintf( "get.autospectral.param.%s", cytometer ) )

  check.critical( ! is.null( get.param.function ), "unsupported cytometer" )

  autosp.param <- get.param.function( autosp.param )

  autosp.param

}

