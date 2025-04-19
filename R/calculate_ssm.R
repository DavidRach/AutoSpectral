# calculate_ssm.r

#' @title calculate.ssm
#'
#' @description Calculates a spillover spread matrix (ssm) for a set of fluorophores.
#'
#' @importFrom MASS ginv
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and 1,
#'     with fluorophores in rows and detectors (channels) in columns. More than one
#'     spectra are needed to calculate anything here.
#' @return The spillover spread matrix (ssm) for the
#'     spectral flow cytometry panel. Spread from row fluor into column fluor.
#' @export

calculate.ssm <- function( spectra ) {

  # calculate SSM
  unmixing.matrix <- ginv( t( spectra ) )

  fluorophores <- rownames( spectra )
  fluorophore.n <- length( fluorophores )

  # initialize empty ssm
  ssm <- matrix( 0, nrow = fluorophore.n, ncol = fluorophore.n,
                 dimnames = list( fluorophores,
                                  fluorophores ) )

  # identify the peak channel for each fluorophore
  peak.channel <- apply( unmixing.matrix, 1, function( x ) which.max( x ) )

  # compute the SSM
  for( i in 1:fluorophore.n ) {
    for( k in 1:fluorophore.n ) {
      if( i != k ) {
        ssm[ i, k ] <- ( unmixing.matrix[ i, peak.channel[ k ] ] /
                           unmixing.matrix[ k, peak.channel[ k ] ] )^2
      }
    }
  }

  return( ssm )
}
