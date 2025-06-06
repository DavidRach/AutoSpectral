# get_autospectral_param_id7000.r

#' Get Autospectral Parameters for ID7000 Cytometer
#'
#' This function returns parameters for running a calculation of unmixing with
#'     AutoSpectral, without creating any figures or tables.
#'
#' @title Get Autospectral Parameters for ID7000 Cytometer
#'
#' @description Returns parameters for running a calculation of unmixing with
#' AutoSpectral, without creating any figures or tables.
#'
#' @param autosp.param A list of initial AutoSpectral parameters.
#'
#' @return A list of AutoSpectral parameters specific to the ID7000 cytometer.
#'
#' @export

get.autospectral.param.id7000 <- function( autosp.param )
{
  # add cytometer-specific parameters
  autosp.param$cytometer <- "ID7000"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 1048576

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 1048576

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 1000000

  autosp.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autosp.param$default.time.parameter <- "TIME"

  autosp.param$default.transformation.param <- list(
    length = 256,
    max.range = 1000000,
    pos = 5,
    neg = 0,
    width = -250
  )

  autosp.param$non.spectral.channel <- c( "TIME", "SSC", "FSC" )

  autosp.param$af.channel <- "405CH7-A"

  autosp.param$data.step <- 1e5

  # spectral parameters

  autosp.param$plot.gate.factor <- 5

  autosp.param$af.density.threshold <- 0.75

  autosp.param$af.gate.param <- list(
    density.threshold = 0.001,
    region.auto = TRUE
  )

  autosp.param$af.figure.gate.scale.expand <- 0.01

  autosp.param$af.gate.bound.density.neigh.size <- 3

  autosp.param$af.gate.bound.density.grid.n <- 200

  autosp.param$af.lower.quantile <- c( 0.01, 0.15 )
  autosp.param$af.upper.quantile <- c( 0.85, 0.99 )
  autosp.param$af.bound.margin <- 2

  # autosp.param$gate.bound.density.grid.n.cells <- 200
  # autosp.param$gate.region.density.grid.n.cells <- 200
  # autosp.param$gate.region.max.density.grid.n.cells <- 200
  # autosp.param$gate.bound.density.grid.n.beads <- 200
  # autosp.param$gate.region.density.grid.n.beads <- 200
  # autosp.param$gate.region.max.density.grid.n.beads <- 200

  autosp.param

}
