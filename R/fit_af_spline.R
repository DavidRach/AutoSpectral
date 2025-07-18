# fit_af_spline.r

#' @title Fit Spline to Autofluorescence Data
#'
#' @description
#' This function fits a spline to autofluorescence data, removing extreme events
#' and defining bounds equally far from zero. It uses robust linear modeling
#' and identifies events within a specified number of standard deviations from
#' the spline.
#'
#' @importFrom MASS rlm
#' @importFrom tripack tri.mesh convex.hull
#' @importFrom stats predict
#' @importFrom sp point.in.polygon
#'
#' @param af.cells A matrix containing the autofluorescence data.
#' @param non.af.cells A matrix containing the low autofluorescence data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return The boundaries of the events within the specified number of standard
#' deviations from the spline.
#'
#' @export

fit.af.spline <- function( af.cells, non.af.cells, asp ) {

  # set up boundaries structure for two gates
  af.boundaries <- list(
    upper = NULL,
    lower = NULL
  )

  # set lower bound based on non-AF cells
  x.bound.low <- max( non.af.cells[ , 1 ] ) * ( 1 + asp$af.spline.x.bound.factor.low )
  y.bound.low <- max( non.af.cells[ , 2 ] ) * ( 1 + asp$af.spline.y.bound.factor.low )

  # fit a spline using rlm
  model.data <- data.frame( rbind( af.cells, non.af.cells ) )

  if ( nrow( model.data ) < 10 )
    stop( "Failed to identify autofluorescence" )

  rlm.fit <- rlm( PC2 ~ PC1, data = model.data, maxit = asp$af.spline.maxit )

  if ( !rlm.fit$converged )
    warning( "The IRLS algorithm employed in 'rlm' did not converge." )

  # define events within n standard deviations of the spline
  predicted <- predict( rlm.fit, newdata = model.data )
  model.data$predicted <- predicted
  model.data$residuals <- model.data$PC2 - predicted

  sd.residuals <- sd( model.data$residuals )

  model.fit <- model.data[ abs( model.data$residuals ) <= asp$af.spline.sd.n * sd.residuals, ]
  model.fit <- model.fit[ which( model.fit$PC1 > x.bound.low ), ]

  # get the boundary of those events
  af.remove.boundary <- convex.hull( tri.mesh( model.fit$PC1, model.fit$PC2 ) )

  if ( length( af.remove.boundary ) == 1 )
    stop( "Failed to identify autofluorescence" )
  else
    af.boundaries$upper <- af.remove.boundary

  if ( asp$af.remove.pop == 1 ) {
    return( af.boundaries )

  } else {
    # gate out first AF population and repeat to find a second AF population
    gate.population.pip <- point.in.polygon(
      model.data[ , 1 ], model.data[ , 2 ],
      af.remove.boundary$x, af.remove.boundary$y )

    gate.population.idx <- which( gate.population.pip == 0 )

    model.data <- model.data[ gate.population.idx, ]

    if ( nrow( model.data ) < 10 )
      return( af.boundaries )

    rlm.fit <- rlm( PC1 ~ PC2, data = model.data, maxit = asp$af.spline.maxit )

    if ( !rlm.fit$converged )
      warning( "The IRLS algorithm employed in 'rlm' did not converge." )

    predicted <- predict( rlm.fit, newdata = model.data[ , c( "PC1", "PC2" ) ] )
    model.data$predicted <- predicted
    model.data$residuals <- model.data$PC1 - predicted

    sd.residuals <- sd( model.data$residuals )

    model.fit <- model.data[ abs( model.data$residuals ) <= asp$af.spline.sd.n * sd.residuals, ]
    model.fit <- model.fit[ which( model.fit$PC2 > y.bound.low ), ]

    if ( nrow( model.fit ) < 10 )
      return( af.boundaries )

    # get the boundary of second AF pop
    af.remove.boundary <- convex.hull( tri.mesh( model.fit$PC1, model.fit$PC2 ) )

    if ( length( af.remove.boundary ) == 1 ) {
      return( af.boundaries )
    } else {
      af.boundaries$lower <- af.remove.boundary
      return( af.boundaries )
    }
  }
}
