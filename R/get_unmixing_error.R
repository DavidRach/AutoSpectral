# get_unmixing_error.r

#' @title Calculate Unmixing Error
#'
#' @description
#' This function calculates the unmixing error for each fluorophore, including
#' intercept, coefficient and slope, using robust linear modeling.
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @param expr.data.unmix Data frame containing unmixed expression data.
#' @param fluorophores Vector of fluorophore names.
#' @param scale.untransformed Logical indicating whether to scale untransformed
#' data.
#' @param transform.inv Function to apply the inverse transformation.
#' @param flow.event.sample Vector indicating the sample for each flow event.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A list containing the unmixing correction matrices for intercept,
#' coefficient and slope.

get.unmixing.error <- function( expr.data.unmix, fluorophores,
                                scale.untransformed, transform.inv,
                                flow.event.sample, asp ){

  if ( asp$parallel ){
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  # set up spectral collection
  fluorophore.n <- length( fluorophores )

  unmixed.spillover.zero <- rep( 0, fluorophore.n )
  names( unmixed.spillover.zero ) <- fluorophores

  unmixed.matrix.corr <- lapply.function( fluorophores, function( fl ) {

    fluorophore.proper <- fl

    unmixed.spillover.corr.inte <- unmixed.spillover.zero
    unmixed.spillover.corr.coef <- unmixed.spillover.zero
    unmixed.spillover.corr.slop <- unmixed.spillover.zero

    for( channel in fluorophores ){

      if ( channel == fluorophore.proper ){

        unmixed.spillover.corr.coef[ channel ] <- 1.0
        unmixed.spillover.corr.slop[ channel ] <- 1.0

      } else {

        channel.expr.proper <- expr.data.unmix[ which( flow.event.sample == fl ),
                                                fluorophore.proper ]

        channel.expr <- expr.data.unmix[ which( flow.event.sample == fl ),
                                         channel ]

        spillover.model.result <- fit.robust.linear.model(
          channel.expr.proper, channel.expr,
          fluorophore.proper, channel, asp$rlm.iter.max )

        unmixed.spillover.corr.inte[ channel ] <- spillover.model.result[ 1 ]
        unmixed.spillover.corr.coef[ channel ] <- spillover.model.result[ 2 ]

        # get slope in untransformed scale
        if ( scale.untransformed ){

          unmixed.spillover.corr.slop[ channel ] <- unmixed.spillover.corr.coef[ channel ]

        } else {

          y1p <- min( channel.expr.proper )
          y2p <- max( channel.expr.proper )

          x1p <- unmixed.spillover.corr.inte[ channel ] +
            unmixed.spillover.corr.coef[ channel ] * y1p
          x2p <- unmixed.spillover.corr.inte[ channel ] +
            unmixed.spillover.corr.coef[ channel ] * y2p

          if ( y1p == y2p || x1p == x2p )
            unmixed.spillover.corr.slop[ channel ] <- 0
          else
          {
            y1 <- transform.inv( y1p )
            y2 <- transform.inv( y2p )

            x1 <- transform.inv( x1p )
            x2 <- transform.inv( x2p )

            unmixed.spillover.corr.slop[ channel ] <-
              unmixed.spillover.corr.coef[ channel ] *
              ( x2 - x1 ) * ( y2p - y1p ) /
              ( ( x2p - x1p ) * ( y2 - y1 ) )
          }
        }
      }
    } # channel

    c( unmixed.spillover.corr.inte, unmixed.spillover.corr.coef,
       unmixed.spillover.corr.slop )

  } )

  unmixed.matrix.corr <- do.call( rbind, unmixed.matrix.corr )
  rownames( unmixed.matrix.corr ) <- fluorophores

  unmixing.corr <- list(
    inte = unmixed.matrix.corr[ , 1 : fluorophore.n ],
    coef = unmixed.matrix.corr[ , 1 : fluorophore.n +
                                  fluorophore.n ],
    slop = unmixed.matrix.corr[ , 1 : fluorophore.n +
                                  2 * fluorophore.n ]
  )

  return( unmixing.corr )
}
