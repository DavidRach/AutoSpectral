# get_autospectral_param_aurora.r

#' Get Autospectral Parameters for Aurora Cytometer
#'
#' This function returns parameters for running a calculation of unmixing with
#'     autospectral, without creating any figures or tables.
#'
#' @title Get Autospectral Parameters for Aurora Cytometer
#' @description Returns parameters for running a calculation of unmixing with
#'     autospectral, without creating any figures or tables.
#' @importsFrom base get0
#' @param autosp.param A list of initial autospectral parameters.
#' @return A list of autospectral parameters specific to the Aurora cytometer.
#' @examples
#' \dontrun{
#' get.autospectral.param.aurora(autosp.param)
#' }


get.autospectral.param.aurora <- function( autosp.param )
{
  # add cytometer-specific parameters

  autosp.param$cytometer <- "Aurora"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 4194304

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 4194304

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 4194304

  autosp.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autosp.param$default.time.parameter <- "Time"

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 4194304,
          pos = 5.62,
          neg = 0,
          width = -1000
        )

  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC" )

  autosp.param$af.channel <- "V7-A"

  autosp.param$data.step <- 5e5

  # spectral parameters

  autosp.param$plot.gate.factor <- 2

  autosp.param$ribbon.breaks <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )

  autosp.param

}

