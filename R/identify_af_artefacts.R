# identify_af_artefacts.r

#' @title Identify Autofluorescence Artefacts
#'
#' @description
#' This function identifies spikes of autofluorescence based on a
#' negative control, defines gates for removal, plots removal gates, plots
#' spectra of cells with and without removal, and returns gate boundaries
#' and the components needed to transform single-stained controls.
#'
#' @importFrom stats prcomp median mad
#' @importFrom sp point.in.polygon
#' @importFrom flowWorkspace flowjo_biexp
#'
#' @param clean.expr List containing cleaned expression data.
#' @param universal.neg Name of the universal negative control.
#' @param spectral.channel Vector of spectral channel names.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A list containing the autofluorescence components and gate boundaries.
#' @export

id.af.artefacts <- function( clean.expr, universal.neg, spectral.channel,
                             asp ) {

  if ( asp$verbose )
    message( paste( "\033[34m", "Identifying autofluorescence contamination in",
                    universal.neg, "\033[0m" ) )

  expr.data.neg <- clean.expr[[ universal.neg ]][ , spectral.channel ]

  if ( nrow( expr.data.neg ) > asp$gate.downsample.n.cells ) {
    set.seed( asp$gate.downsample.seed )
    downsample.idx <- sample( nrow( expr.data.neg ), asp$gate.downsample.n.cells )
    expr.data.neg <- expr.data.neg[ downsample.idx, ]
  }

  expr.data.center <- apply( expr.data.neg, 2, median )
  expr.data.scale <- apply( expr.data.neg, 2, mad )

  scaled.data <- scale( expr.data.neg, center = expr.data.center,
                        scale = expr.data.scale )

  af.pca <- prcomp( scaled.data, center = FALSE, scale. = FALSE )

  # use the first two components to identify main AF signatures to be removed
  af.components <- apply( af.pca$rotation[ , 1:2 ], 2, function( x ) x / max( x ) )
  af.components <- t( apply( af.components, 2, function( x ) x / max( x ) ) )

  # unmix using first two components
  gate.data <- unmix.ols( expr.data.neg, af.components )

  biexp.transform <- flowjo_biexp( channelRange = asp$default.transformation.param$length,
                                   maxValue = asp$default.transformation.param$max.range,
                                   pos = asp$default.transformation.param$pos,
                                   neg = asp$default.transformation.param$neg,
                                   widthBasis = asp$default.transformation.param$width,
                                   inverse = FALSE )

  # transform for easier visualization and faster gating
  gate.data <- apply( gate.data, 2, biexp.transform )

  # get gate defining low AF region
  af.gate.idx <- do.gate.af( gate.data, universal.neg, asp )

  af.cells <- gate.data[ -af.gate.idx, ]

  # include up to 500 non-AF cells
  if ( length( af.gate.idx ) > 500 ) {
    set.seed( asp$gate.downsample.seed )
    sample.idx <- sample( af.gate.idx, 500 )
    non.af.cells <- gate.data[ sample.idx, ]
  } else {
    non.af.cells <- gate.data[ af.gate.idx, ]
  }

  af.boundaries <- fit.af.spline( af.cells, non.af.cells, asp )

  gate.population.pip.upper <- point.in.polygon(
    gate.data[ , 1 ], gate.data[ , 2 ],
    af.boundaries$upper$x, af.boundaries$upper$y )

  if ( asp$af.remove.pop != 1 & !is.null( af.boundaries$lower ) ) {
    # exclude both AF pops
    gate.population.pip.lower <- point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      af.boundaries$lower$x, af.boundaries$lower$y )

    gate.population.idx <- which( gate.population.pip.lower == 0 &
                                    gate.population.pip.upper == 0 )

  } else {
    gate.population.idx <- which( gate.population.pip.upper == 0 )
  }

  if ( asp$figures ) {

    if ( asp$af.remove.pop != 1 & !is.null( af.boundaries$lower ) ){
      gate.af.sample.plot( universal.neg, af.data = gate.data,
                           af.boundaries$lower,
                           af.boundaries$upper,
                           asp )
    } else {
      gate.af.sample.plot( universal.neg, af.data = gate.data,
                           af.boundary.lower = NULL,
                           af.boundary.upper = af.boundaries$upper,
                           asp )
    }

    spectral.ribbon.plot( expr.data.neg, expr.data.neg[ gate.population.idx, ],
                          spectral.channel, asp, universal.neg,
                          plot.prefix = asp$af.plot.filename,
                          af = TRUE,
                          removed.data = expr.data.neg[ -gate.population.idx, ] )
  }

  return( list( af.components, af.boundaries) )

}
