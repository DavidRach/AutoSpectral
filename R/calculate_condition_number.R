# calculate_condition.number.r

#' @title Calculate Condition Number
#'
#' @description
#' Calculates the mixing matrix condition number ("Complexity IndexTM") for a
#' set of fluorophores in a spectral panel.
#'
#' @param spectra A matrix of spectral signatures of fluorophores, normalized
#' between 0 and 1.
#' Fluorophores should be in rows and detectors (channels) in columns.
#' More than one spectrum is required to perform the calculation.
#'
#' @return A numeric value representing the condition number for the spectral
#' flow cytometry panel.
#'
#' @export


calculate.condition.number <- function( spectra ) {

  svd.result <- svd( spectra )

  singular.values <- svd.result$d

  condition.number <- max( singular.values ) / min( singular.values )

  condition.number <- round( condition.number, 2 )

  return( condition.number )
}
