# get_autospectral_param_opteon.r

#' Get Autospectral Parameters for Opteon Cytometer
#'
#' This function returns parameters for running a calculation of unmixing with
#'     autospectral, without creating any figures or tables.
#'
#' @title Get Autospectral Parameters for Opteon Cytometer
#' @description Returns parameters for running a calculation of unmixing with
#'     autospectral, without creating any figures or tables.
#' @param autosp.param A list of initial autospectral parameters.
#' @return A list of autospectral parameters specific to the Opteon cytometer.
#' @export

get.autospectral.param.opteon <- function( autosp.param )
{
  # add cytometer-specific parameters

  autosp.param$cytometer <- "Opteon"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 13841000

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 8687233

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 16777216

  autosp.param$default.scatter.parameter <- c( "FSC-A", "VSSC-A" )

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 100000,
          pos = 6.68,
          neg = 0,
          width = -3000
        )

  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC", "-H", "Width" )

  autosp.param$af.channel <- "V525-A"

  autosp.param$data.step <- 3e6

  # spectral parameters

  autosp.param$plot.gate.factor <- 2

  autosp.param$af.density.threshold <- 0.75

  autosp.param$af.gate.param < list(
          density.threshold = 0.001,
          region.auto = TRUE
        )

  autosp.param$af.figure.gate.scale.expand <- 0.01

  autosp.param$af.gate.bound.density.neigh.size <- 3

  autosp.param$af.gate.bound.density.grid.n <- 100

  autosp.param$af.lower.quantile <- c( 0.01, 0.15 )
  autosp.param$af.upper.quantile <- c( 0.85, 0.99 )
  autosp.param$af.bound.margin <- 2

  autosp.param

}

