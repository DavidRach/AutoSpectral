# calculate_complexity_index.r

#' @title Calculate Complexity Index
#'
#' @description
#' Calculates the complexity index (condition number) for a set of
#' fluorophores in a spectral panel.
#'
#' @param spectra A matrix of spectral signatures of fluorophores, normalized
#' between 0 and 1.
#' Fluorophores should be in rows and detectors (channels) in columns.
#' More than one spectrum is required to perform the calculation.
#'
#' @return A numeric value representing the complexity index (condition number)
#' for the spectral flow cytometry panel.
#'
#' @export


calculate.complexity.index <- function( spectra ) {

  svd.result <- svd( spectra )

  singular.values <- svd.result$d

  complexity.index <- max( singular.values ) / min( singular.values )

  complexity.index <- round( complexity.index, 2 )

  return( complexity.index )
}
