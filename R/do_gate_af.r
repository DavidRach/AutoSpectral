# do_gate_af.r

#' @title Perform Gating on Autofluorescence Parameters
#'
#' @description
#' This function returns a vector with the indexes of events inside the initial
#' gate on autofluorescence parameters. It proceeds through several steps to
#' define the gate boundaries and identify density maxima using numerical
#' search and Voronoi tessellations.
#'
#' @importFrom deldir deldir tile.list which.tile
#' @importFrom sp point.in.polygon
#' @importFrom tripack tri.mesh convex.hull
#' @importFrom MASS kde2d bandwidth.nrd
#'
#' @param gate.data A data frame containing the gate data.
#' @param samp A sample identifier.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A vector with the indexes of events inside the initial gate.
#'
#' @export

do.gate.af <- function( gate.data, samp, asp ) {

  gate.marker <- colnames( gate.data )

  gate.bound <- NULL
  gate.region <- NULL
  gate.boundary <- NULL

  gate.data.x.min <- min( gate.data[ , 1 ] )
  gate.data.x.max <- max( gate.data[ , 1 ] )

  gate.data.y.min <- min( gate.data[ , 2 ] )
  gate.data.y.max <- max( gate.data[ , 2 ] )

  grid.n <- asp$af.gate.bound.density.grid.n

  gate.bound.density <- MASS::kde2d(
    gate.data[ , 1 ], gate.data[ , 2 ],
    asp$af.gate.density.bw.factor * apply( gate.data, 2, bandwidth.nrd ),
    n = grid.n )

  # Identify density maxima
  neighbor.n <- asp$af.gate.bound.density.neigh.size

  gate.bound.neighbor.idx <- list(
    x = - neighbor.n :neighbor.n,
    y = - neighbor.n :neighbor.n )

  gate.bound.density.max.bool <- matrix( FALSE, nrow = grid.n, ncol = grid.n )

  for ( x.idx in 1 : grid.n )
    for ( y.idx in 1 : grid.n )
      gate.bound.density.max.bool[ x.idx, y.idx ] <-
    gate.bound.density$z[ x.idx, y.idx ] >=
    max( gate.bound.density$z[
      pmax( 0, pmin( grid.n, x.idx + gate.bound.neighbor.idx$x ) ),
      pmax( 0, pmin( grid.n, y.idx + gate.bound.neighbor.idx$y ) ) ] )

  gate.bound.density.max.idx <- which( gate.bound.density.max.bool,
                                       arr.ind = TRUE )

  gate.bound.density.max.n <- nrow( gate.bound.density.max.idx )

  check.critical( gate.bound.density.max.n >= 1,
                  paste0( "gate error: no population found in sample bound", samp ) )

  gate.bound.density.max <- data.frame(
    x = gate.bound.density$x[ gate.bound.density.max.idx[ , 1 ] ],
    y = gate.bound.density$y[ gate.bound.density.max.idx[ , 2 ] ],
    z = gate.bound.density$z[ gate.bound.density.max.idx ] )

  gate.bound.density.max <- gate.bound.density.max[
    order( gate.bound.density.max$z, decreasing = TRUE ), ]

  row.names( gate.bound.density.max ) <- NULL
  gate.bound.density.max$num.label <- paste0( " ",
                                              row.names( gate.bound.density.max ) )

  # Identify the target maximum
  target.max.idx <- asp$af.gate.target.max
  target.max <- gate.bound.density.max[ target.max.idx, ]

  gate.bound.voronoi <- deldir( gate.bound.density.max,
                                rw = c( gate.data.x.min, gate.data.x.max, gate.data.y.min,
                                        gate.data.y.max ), suppressMsge = TRUE )

  tiles <- tile.list( gate.bound.voronoi )

  gate.bound <- list(
    density = gate.bound.density,
    density.max = gate.bound.density.max,
    density.max.n = gate.bound.density.max.n,
    density.max.data.idx = target.max.idx,
    density.max.target = target.max,
    voronoi = gate.bound.voronoi,
    x.low = gate.data.x.min,
    x.high = gate.data.x.max,
    y.low = gate.data.y.min,
    y.high = gate.data.y.max
  )

  # Identify points in the tile of the global maximum
  target.tile <- tiles[[ target.max.idx ]]
  tile.polygon <- cbind( target.tile$x, target.tile$y )

  gate.region.data.idx <- which( point.in.polygon(
    gate.data[, 1], gate.data[, 2],
    tile.polygon[ , 1 ], tile.polygon[ , 2 ] ) > 0 )

  check.critical( length( gate.region.data.idx ) > 0,
                  paste( "Error: No points in target Voronoi tile for sample",
                         samp ) )

  # define region boundaries
  gate.region.data <- gate.data[ gate.region.data.idx, ]

  gate.region.x.low <- min( gate.region.data[ , 1 ] )
  gate.region.x.high <- max( gate.region.data[ , 1 ]  )
  gate.region.y.low <- min( gate.region.data[ , 2 ] )
  gate.region.y.high <- max( gate.region.data[ , 2 ] )

  # get density maxima in region
  gate.region.density <- MASS::kde2d(
    gate.data[ gate.region.data.idx, 1 ], gate.data[ gate.region.data.idx, 2 ],
    asp$af.gate.density.bw.factor * apply( gate.data[ gate.region.data.idx, ], 2,
                                           bandwidth.nrd ),
    n = grid.n )

  gate.region.neighbor.idx <- list( x = - neighbor.n : neighbor.n,
                                    y = - neighbor.n : neighbor.n )

  gate.region.density.max.bool <- matrix( FALSE, nrow = grid.n, ncol = grid.n )

  for ( x.idx in 1 : grid.n )
    for ( y.idx in 1 : grid.n )
      gate.region.density.max.bool[ x.idx, y.idx ] <-
    gate.region.density$z[ x.idx, y.idx ] >=
    max( gate.region.density$z[
      pmax( 0, pmin( grid.n, x.idx + gate.region.neighbor.idx$x ) ),
      pmax( 0, pmin( grid.n, y.idx + gate.region.neighbor.idx$y ) ) ] )

  gate.region.density.max.idx <- which( gate.region.density.max.bool,
                                        arr.ind = TRUE )

  gate.region.density.max.n <- nrow( gate.region.density.max.idx )

  check.critical( gate.region.density.max.n >= 1,
                  paste( "gate error: no population found in sample region", samp ) )

  gate.region.density.max <- data.frame(
    x = gate.region.density$x[ gate.region.density.max.idx[ , 1 ] ],
    y = gate.region.density$y[ gate.region.density.max.idx[ , 2 ] ],
    z = gate.region.density$z[ gate.region.density.max.idx ] )

  gate.region.density.max <- gate.region.density.max[
    order( gate.region.density.max$z, decreasing = TRUE ), ]

  row.names( gate.region.density.max ) <- NULL
  gate.region.density.max$num.label <- paste0( " ",
                                               row.names( gate.region.density.max ) )

  if ( gate.region.density.max.n > 1 ) {

    # get voronoi tesselation for density maxima
    gate.region.voronoi <- deldir( gate.region.density.max,
                                   rw = c( gate.region.x.low, gate.region.x.high,
                                           gate.region.y.low, gate.region.y.high ),
                                   suppressMsge = TRUE )

    gate.region.tile <- tile.list( gate.region.voronoi )

    # get data in the tile of largest maximum
    gate.region.density.max.data.idx <- gate.region.data.idx[
      sapply( gate.region.data.idx, function( grdi )
        which.tile( gate.data[ grdi, 1 ], gate.data[ grdi, 2 ],
                    gate.region.tile ) == 1 )
    ]
  } else {
    gate.region.voronoi <- NULL
    gate.region.density.max.data.idx <- gate.region.data.idx

  }

  gate.region <- list(
    data.idx = gate.region.data.idx,
    density = gate.region.density,
    density.max = gate.region.density.max,
    density.max.n = gate.region.density.max.n,
    density.max.data.idx = gate.region.density.max.data.idx,
    voronoi = gate.region.voronoi,
    x.low = gate.region.x.low,
    x.high = gate.region.x.high,
    y.low = gate.region.y.low,
    y.high = gate.region.y.high
  )

  # threshold data in region around target maximum
  gate.region.max.density <- MASS::kde2d(
    gate.data[ gate.region.density.max.data.idx, 1 ],
    gate.data[ gate.region.density.max.data.idx, 2 ],
    asp$af.gate.density.bw.factor *
      apply( gate.data[ gate.region.density.max.data.idx, ], 2,
             bandwidth.nrd ),
    n = grid.n )

  gate.region.max.density.interp <- interp.surface( gate.region.max.density,
                                                    gate.data[ gate.region.density.max.data.idx, ] )

  gate.region.max.density.threshold <-
    ( 1 - asp$af.gate.param$density.threshold ) *
    min( gate.region.max.density.interp ) +
    asp$af.gate.param$density.threshold * max( gate.region.max.density.interp )

  gate.population.strict.idx <- gate.region.density.max.data.idx[
    gate.region.max.density.interp > gate.region.max.density.threshold ]

  gate.population.strict.idx <- gate.population.strict.idx[
    ! duplicated( data.frame( gate.data[ gate.population.strict.idx, ] ) ) ]

  gate.population.boundary <- convex.hull( tri.mesh(
    gate.data[ gate.population.strict.idx, 1 ],
    gate.data[ gate.population.strict.idx, 2 ] ) )

  gate.population.pip <- point.in.polygon(
    gate.data[ , 1 ], gate.data[ , 2 ],
    gate.population.boundary$x, gate.population.boundary$y )

  gate.population <- list( boundary = gate.population.boundary )

  if ( ! is.null( asp$figure.clean.control.dir ) )
    gate.af.plot( samp, gate.data, gate.bound, gate.region,
                  gate.population, asp )

  if ( asp$af.gate.bound.strict ) {

    gate.population.idx <- which( gate.population.pip != 0 )

  } else {

    gate.population.idx <- which( gate.data[ , 1 ] > gate.region.x.low &
                                    gate.data[ , 1 ] < gate.region.x.high &
                                    gate.data[ , 2 ] > gate.region.y.low &
                                    gate.data[ , 2 ] < gate.region.y.high )
  }

  return( gate.population.idx )
}
