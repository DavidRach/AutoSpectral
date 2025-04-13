# unmix_wls.r

#' unmix.wls
#'
#' Performs spectral unmixing using weighted least squares.
#' Weighting is by channel power
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and detectors
#' (raw channels) in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and 1,
#' with fluorophores in rows and detectors (channels) in columns.
#' @return Unmixed data with cells in rows and fluorophores in columns.
#' @export


unmix.wls <- function( raw.data, spectra ) {

  spectra <- t( spectra )

  channel.var <- apply( raw.data, 2, sd )

  # weights are inverse of channel variances
  channel.weights <- 1 / ( channel.var + 1e-6 )

  W <- diag( channel.weights )

  # Weighted LS solution: (M^T W M)^{-1} M^T W
  unmixing.matrix <- solve( t( spectra ) %*% W %*% spectra ) %*%
    ( t( spectra ) %*% W )

  unmixed.data <- raw.data %*% t( unmixing.matrix )

  colnames( unmixed.data ) <- colnames( spectra )

  return( unmixed.data )

}
