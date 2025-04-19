 # fit_af_spline.r


#' Fit Spline to Autofluorescence Data
#'
#' This function fits a spline to autofluorescence data, removing extreme events
#'     and defining bounds equally far from zero. It uses robust linear modeling
#'     and identifies events within a specified number of standard deviations
#'     from the spline.
#'
#' @title Fit Spline to Autofluorescence Data
#' @description Fits a spline to autofluorescence data, removing extreme events
#'     and defining bounds equally far from zero.
#'
#' @importFrom MASS rlm
#' @importFrom tripack tri.mesh convex.hull
#' @importFrom stats predict
#' @param af.data A matrix containing the autofluorescence data.
#' @param af.gate.idx A vector of indexes indicating the gate.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#'     A list of essential parameters.
#' @return A list with the boundaries of the events within the specified number
#'     of standard deviations from the spline.
#' @export



fit.af.spline <- function( af.data, af.gate.idx, asp ){

  af.cells <- af.data[ -af.gate.idx, ]

  af.cells.upper <- af.cells[ af.cells[ , 2 ] > 0, ]
  af.cells.lower <- af.cells[ af.cells[ , 2 ] < 0, ]

  # make y-values in lower positive to allow bounds to be set equally far from zero
  af.cells.lower[ , 2 ] <- af.cells.lower[ , 2 ] * -1

  af.split <- list(
    lower = af.cells.lower,
    upper = af.cells.upper
  )

  af.remove.bounds <- lapply( af.split, function( af ) {

    # trim data to remove extreme events and low (less AF) events
    x.range <- range( af[ , 1 ] )
    x.range <- abs( x.range[ 1 ] - x.range[ 2 ] )
    y.range <- range( af[ , 2 ] )
    y.range <- abs( y.range[ 1 ] - y.range[ 2 ] )

    x.bound.low <- min( af[ , 1 ] ) + ( x.range * asp$af.spline.x.bound.factor.low )
    x.bound.high <- max( af[ , 1 ] ) * asp$af.spline.x.bound.factor.high
    y.bound.low <- min( af[ , 2 ] ) + ( y.range * asp$af.spline.y.bound.factor.low )
    y.bound.high <- max( af[ , 2 ] ) * asp$af.spline.y.bound.factor.high

    # fit a spline using rlm
    model.data <- data.frame( af[
        af[ , 1 ] > x.bound.low & af[ , 1 ] < x.bound.high &
        af[ , 2 ] > y.bound.low & af[ , 2 ] < y.bound.high , ] )

    rlm.fit <- rlm( PC2 ~ PC1, data = model.data, maxit = asp$af.spline.maxit )

    if ( !rlm.fit$converged ) {
      warning( "The IRLS algorithm employed in 'rlm' did not converge." )
    }

    # define events within n standard deviations of the spline
    predicted <- predict( rlm.fit, newdata = model.data )
    model.data$predicted <- predicted
    model.data$residuals <- model.data$PC2 - predicted

    sd.residuals <- sd( model.data$residuals )

    model.fit <- model.data[ abs( model.data$residuals ) <= asp$af.spline.sd.n * sd.residuals, ]

    # get the boundary of those events
    af.remove.boundary <- convex.hull( tri.mesh( model.fit$PC1, model.fit$PC2 ) )

    af.remove.boundary
  } )

  # re-convert lower y to negative values
  names( af.remove.bounds ) <- c( "lower", "upper" )

  af.remove.bounds$lower$y <- af.remove.bounds$lower$y * -1

  return( af.remove.bounds )

}
