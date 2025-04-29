# identify_af_artefacts.r

# identifies spikes of autofluorescence based on a negative control
# defines gates for removal
# plots removal gates
# plots spectra of cells with and without removal
# returns gate boundaries and the components needed to transform single-stained controls


#' @title Identify Autofluorescence Artefacts
#'
#' @description This function identifies spikes of autofluorescence based on a
#'     negative control, defines gates for removal, plots removal gates, plots
#'     spectra of cells with and without removal, and returns gate boundaries
#'     and the components needed to transform single-stained controls.
#'
#' @importFrom stats prcomp median mad
#' @importFrom sp point.in.polygon
#'
#' @param clean.expr List containing cleaned expression data.
#' @param universal.neg Name of the universal negative control.
#' @param spectral.channel Vector of spectral channel names.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#'
#' @return A list containing the autofluorescence components and gate boundaries.
#' @export



identify.af.artefacts <- function( clean.expr, universal.neg,
                                   spectral.channel, asp ) {

  if ( asp$verbose )
    cat( paste( "\033[34m", "Identifying autofluorescence contamination in",
                universal.neg, "\033[0m", "\n" ) )

  expr.data.neg <- clean.expr[[ universal.neg ]][ , spectral.channel ]

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

  # get gate defining low AF region
  af.gate.idx <- do.gate.af( gate.data, universal.neg, asp )

  af.cells <- gate.data[ -af.gate.idx, ]

  af.boundaries <- fit.af.spline( gate.data, af.gate.idx, asp )

  gate.population.pip.lower <- point.in.polygon(
    gate.data[ , 1 ], gate.data[ , 2 ],
    af.boundaries$lower$x, af.boundaries$lower$y )

  gate.population.pip.upper <- point.in.polygon(
    gate.data[ , 1 ], gate.data[ , 2 ],
    af.boundaries$upper$x, af.boundaries$upper$y )

  if ( asp$af.remove.pop != 1 ) {
    # return both boundaries
    gate.population.idx <- which( gate.population.pip.lower == 0 &
                                    gate.population.pip.upper == 0 )

  } else {
    # return boundary that excludes more cells
    pip.lower.sum <- sum( gate.population.pip.lower )

    pip.upper.sum <- sum( gate.population.pip.upper )

    if ( pip.upper.sum >= pip.lower.sum ) {

      gate.population.idx <- which( gate.population.pip.upper == 0 )

      af.boundaries <- af.boundaries[ "upper" ]

    } else {

      gate.population.idx <- which( gate.population.pip.lower == 0 )

      af.boundaries <- af.boundaries[ "lower" ]
      # rename for ease of tracking later
      names( af.boundaries ) <- gsub( "lower", "upper", names( af.boundaries ) )

    }
  }

  if ( !is.null( asp$figure.clean.control.dir ) ) {

    if ( length( af.boundaries ) == 2 ){

      plot.gate.af.sample( universal.neg, af.data = gate.data,
                           af.boundaries$lower, af.boundaries$upper,
                           asp )

    } else {

      plot.gate.af.sample( universal.neg, af.data = gate.data,
                           af.boundary.lower = NULL,
                           af.boundary.upper = af.boundaries$upper,
                           asp )

      }

  }

  if ( !is.null( asp$figure.clean.control.dir ) ) {

    spectral.ribbon.plot( expr.data.neg, expr.data.neg[ gate.population.idx, ],
                          spectral.channel, asp, universal.neg,
                          plot.prefix = asp$af.plot.filename,
                          af = TRUE,
                          removed.data = expr.data.neg[ -gate.population.idx, ] )
  }

  return( list( af.components, af.boundaries) )

}
