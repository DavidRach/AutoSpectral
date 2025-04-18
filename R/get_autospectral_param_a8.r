# get_autospectral_param_a8.r

#' Get Autospectral Parameters for DiscoverA8 Cytometer
#'
#' This function returns parameters for running a calculation of unmixing with
#'     autospectral, without creating any figures or tables.
#'
#' @title Get Autospectral Parameters for DiscoverA8 Cytometer
#' @description Returns parameters for running a calculation of unmixing with
#'     autospectral, without creating any figures or tables.
#' @importsFrom base get0
#' @param autosp.param A list of initial autospectral parameters.
#' @return A list of autospectral parameters specific to the DiscoverA8 cytometer.
#' @examples
#' \dontrun{
#' get.autospectral.param.a8(autosp.param)
#' }


get.autospectral.param.a8 <- function( autosp.param )
{
  # add cytometer-specific parameters

  autosp.param$cytometer <- "DiscoverA8"

  autosp.param$scatter.data.min.x <- 0

  autosp.param$scatter.data.max.x <- 24140237

  autosp.param$scatter.data.min.y <- 0

  autosp.param$scatter.data.max.y <- 16488107

  autosp.param$expr.data.min <- -111

  autosp.param$expr.data.max <- 24140237

  autosp.param$large.gate.quantile <- 0.25
  autosp.param$large.gate.scaling.x <- 2.5
  autosp.param$large.gate.scaling.y <- 6

  autosp.param$default.scatter.parameter <- c( "FSC-A", "SSC (Violet)-A" )

  autosp.param$default.transformation.param <- list(
          length = 256,
          max.range = 2147483648.0,
          pos = 8.33,
          neg = 0,
          width = -500
        )

  autosp.param$time.and.scatter <- c( "FSC-A", "FSC-H", "FSC-W",
                                      "SSC (Violet)-A", "SSC (Violet)-H", "SSC (Violet)-W",
                                      "SSC (Imaging)-A", "SSC (Imaging)-H", "SSC (Imaging)-W",
                                      "LightLoss (Imaging)-A", "LightLoss (Imaging)-H","LightLoss (Imaging)-W",
                                      "LightLoss (Violet)-A", "LightLoss (Violet)-H", "LightLoss (Violet)-W"  )

  autosp.param$non.spectral.channel <- c( "Time", "SSC", "FSC", "-H", "-W", "-T",
                                          "Delta", "Plate",
                                          "Radial", "Correlation", "Intensity",
                                          "Eccentricity", "Diffusivity", "Center",
                                          "Moment", "Size", "LightLoss",
                                          "Saturated", "Sorted", "Row", "Column",
                                          "Img", "Protocol", "EventLabel",
                                          "Region", "Gate", "Index", "Phase",
                                          "Event", "Drop", "Spectral", "Waveform",
                                          "Merged", "Flow", "Packet", "Reserved" )

  autosp.param$af.channel <- "V6 (515)-A"

  autosp.param$data.step <- 5e6

  # spectral parameters

  autosp.param$plot.gate.factor <- 5

  autosp.param$af.density.threshold <- 0.75

  autosp.param$af.gate.param <- list(
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

