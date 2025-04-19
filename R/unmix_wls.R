# unmix_wls.r
#'
#' @title Unmix Using Weighted Least Squares
#'
#' @description This function performs unmixing of raw data using weighted least
#'    squares (WLS) based on the provided spectra. Weighting is by channel power.
#'
#' @importFrom stats sd
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and detectors
#'     (raw channels) in columns. Columns should be fluorescent data only and must
#'     match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and 1,
#'     with fluorophores in rows and detectors (channels) in columns.
#' @return A matrix containing unnmixed data with cells in rows and fluorophores in columns.
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
