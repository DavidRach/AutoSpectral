# unmix_ols.r

#' @title unmix.ols
#'
#' @description
#' Performs spectral unmixing using ordinary least squares
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.ols <- function( raw.data, spectra ) {

  spectra <- t( spectra )

  unmixing.matrix <- solve( crossprod( spectra ) ) %*% t( spectra )

  unmixed.data <- raw.data %*% t( unmixing.matrix )

  colnames( unmixed.data ) <- colnames( spectra )

  return( unmixed.data )

}
