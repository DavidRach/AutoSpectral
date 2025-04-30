# get_universal_negative.r


#' @title Get Universal Negative Control
#'
#' @description This function identifies and processes the universal negative
#'     control for a given sample, including scatter matching and plotting.
#'
#' @importFrom stats quantile median mad
#' @importFrom KernSmooth bkde2D dpik
#' @importFrom sp point.in.polygon Polygon Polygons SpatialPolygons
#' @importFrom grDevices contourLines
#'
#' @param clean.expr.data List containing cleaned expression data.
#' @param samp Sample identifier.
#' @param universal.negatives Named vector mapping samples to their universal
#'     negatives.
#' @param scatter.param Vector of scatter parameters.
#' @param peak.channels Named vector mapping samples to their peak channels.
#' @param downsample Logical indicating whether to downsample the data.
#' @param negative.n Number of negative events to select.
#' @param positive.n Number of positive events to select.
#' @param spectral.channel Vector of spectral channel names.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param control.type Named vector mapping samples to their control types.
#' @param scatter.match Logical indicating whether to perform scatter matching.
#'
#' @return A data frame containing the selected positive and scatter-matched
#'     negative events.
#' @export



get.universal.negative <- function( clean.expr.data, samp,
                                    universal.negatives,
                                    scatter.param,
                                    peak.channels, downsample,
                                    negative.n, positive.n,
                                    spectral.channel, asp,
                                    control.type,
                                    scatter.match = TRUE ){

  if ( asp$verbose )
    cat( paste( "\033[34m", "Getting universal negative for", samp, "\033[0m", "\n" ) )

  pos.control.expr <- clean.expr.data[[ samp ]]

  neg.control.expr <- clean.expr.data[[ universal.negatives[ samp ] ]]

  peak.channel <- peak.channels[ samp ]

  pos.peak.channel <- pos.control.expr[ , peak.channel ]

  neg.peak.channel <- neg.control.expr[ , peak.channel ]

  # define positive events as those above a threshold (default 99.5%) in the negative
  if ( samp == "AF" ) {
    threshold <- asp$positivity.threshold.af
  } else {
    threshold <- asp$positivity.threshold
  }

  positivity.threshold <- quantile( neg.peak.channel, threshold )

  pos.above.threshold <- pos.peak.channel[ pos.peak.channel > positivity.threshold ]

  # warn if few events in positive
  if ( length( pos.above.threshold ) < asp$min.cell.warning.n ){
    cat( paste( "\033[31m", "Warning! Fewer than",  asp$min.cell.warning.n,
    "gated events in", samp ),  "\033[0m", "\n"  )
  }

  # stop if fewer than minimum acceptable events, returning original data
  if ( length( pos.above.threshold ) < asp$min.cell.stop.n ){
    return( rbind( pos.control.expr, neg.control.expr ) )
  }

  check.critical( length( pos.above.threshold > asp$min.cell.stop.n ),
                  paste( "\033[31m", "Error! Fewer than", asp$min.cell.warning.n,
                  "events remain in", samp, "\033[0m", "\n" ) )

  # select only brightest positive.n events
  if ( length( pos.above.threshold ) >= positive.n ){
    pos.selected <- sort( pos.above.threshold, decreasing = TRUE )[ 1:positive.n ]
  } else {
    pos.selected <- pos.above.threshold
  }

  # scatter-match negative
  # recover full data
  pos.selected.expr <- pos.control.expr[ names( pos.selected ), ]

  # find scatter-matched events in the universal negative
  # if using beads, default to no matching
  sample.control.type <- control.type[[ samp ]]

  if ( sample.control.type == "beads" ) {
    scatter.match <- FALSE
  }

  if ( scatter.match ) {

    pos.scatter.coord <- pos.selected.expr[ , scatter.param ]

    bdw.x <- asp$gate.bound.density.bw.factor.cells * dpik( pos.scatter.coord[ , 1 ] )
    bdw.y <- asp$gate.bound.density.bw.factor.cells * dpik( pos.scatter.coord[ , 2 ] )

    pos.bound.density <- bkde2D( pos.scatter.coord,
                                 bandwidth = c( bdw.x, bdw.y ),
                                 gridsize = c( asp$gate.bound.density.grid.n.cells,
                                               asp$gate.bound.density.grid.n.cells ) )

    names( pos.bound.density ) <- c( "x", "y", "z" )

    contour.levels <- contourLines(pos.bound.density$x, pos.bound.density$y,
                                   pos.bound.density$z,
                                   levels = quantile(pos.bound.density$z,
                                                     probs = asp$scatter.match.threshold ))

    contour.coords <- cbind( contour.levels[[ 1 ]]$x, contour.levels[[ 1 ]]$y )

    contour.polygon <- Polygon( contour.coords )
    contour.polygon <- Polygons( list( contour.polygon ), ID = "density contour" )
    contour.polygon <- SpatialPolygons( list( contour.polygon ) )

    points.within.contour <- pos.scatter.coord[ point.in.polygon( pos.scatter.coord[ , 1 ],
                                                                  pos.scatter.coord[ , 2 ],
                                                                  contour.levels[[ 1 ]]$x,
                                                                  contour.levels[[ 1 ]]$y ) == 1, ]


    pos.scatter.gate <- convex.hull( tri.mesh(
      points.within.contour[ , 1 ],
      points.within.contour[ , 2 ] ) )

    neg.scatter.matched.pip <- point.in.polygon(
      neg.control.expr[ , scatter.param[ 1 ] ],
      neg.control.expr[ , scatter.param[ 2 ] ],
      pos.scatter.gate$x, pos.scatter.gate$y )

    neg.population.idx <- which( neg.scatter.matched.pip != 0 )

    # warn if few events in negative
    if ( length( neg.population.idx ) < asp$min.cell.warning.n ){
      cat( paste( "\033[31m", "Warning! Fewer than",  asp$min.cell.warning.n,
                  "scatter-matched negative events for", samp ),  "\033[0m", "\n"  )
    }

    # stop if fewer than minimum acceptable events, returning original negative
    if ( length( neg.population.idx ) < asp$min.cell.stop.n ){

      # downsample original negative
      if ( nrow( neg.control.expr ) > negative.n ) {
        neg.population.idx <- sample( nrow( neg.control.expr ), negative.n )
        neg.control.expr <- neg.control.expr[ neg.population.idx, ]
      }

      return( rbind( pos.selected.expr, neg.control.expr ) )
    }

    if ( length( neg.population.idx ) > negative.n )
      neg.population.idx <- sample( neg.population.idx, negative.n )

    neg.scatter.matched <- neg.control.expr[ neg.population.idx, ]

  } else {
    neg.population.idx <- sample( 1:nrow( neg.control.expr ), negative.n )
    neg.scatter.matched <- neg.control.expr[ neg.population.idx, ]
  }


  if ( !is.null( asp$figure.spectral.ribbon.dir ) ){
    spectral.ribbon.plot( pos.selected.expr, neg.scatter.matched,
                               spectral.channel, asp, samp )
  }

  if ( !is.null( asp$figure.scatter.dir.base ) ){
    scatter.match.plot( pos.selected.expr, neg.scatter.matched, samp,
                        scatter.param, asp )
  }

  rbind( pos.selected.expr, neg.scatter.matched )

}
