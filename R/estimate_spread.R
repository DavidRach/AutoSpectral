# estimate_spread.r

#' @title Estimate Spread
#'
#' @description
#' A crude approach to predicting spillover spread in the unmixed space based
#' on the spectral matrix.
#'
#' @importFrom stats runif
#'
#' @param spectra The spectral matrix, dimensions fluorophores x detectors.
#'
#' @return The estimated spread matrix
#'
#' @export

estimate.spread <- function( spectra ) {
  fluorophore.n <- nrow( spectra )
  detector.n <- ncol( spectra )

  spectra.t <- t( spectra )
  signal.variance <- runif( fluorophore.n, min = 1e3, max = 1e5 )

  detector.variance <- ( spectra.t^2 ) *
    matrix( signal.variance,
           nrow = detector.n,
           ncol = fluorophore.n,
           byrow = TRUE )

  unmixing.matrix <- solve( t( spectra.t) %*% spectra.t ) %*% t( spectra.t )

  total.detector.variance <- rowSums( detector.variance )
  variance.matrix <- diag( total.detector.variance )

  fluorophore.variance.matrix <- unmixing.matrix %*% variance.matrix %*% t(unmixing.matrix)

  predicted.spread.matrix <- sqrt( abs( fluorophore.variance.matrix ) )
  diag( predicted.spread.matrix ) <- 0

  normalization.matrix <- matrix( signal.variance,
                                 nrow = fluorophore.n,
                                 ncol = fluorophore.n,
                                 byrow = FALSE )

  estimated.spread.matrix <- predicted.spread.matrix / normalization.matrix

  return( estimated.spread.matrix )
}
