# calculate_complexity_index.r

#' @title calculate.complexity.index
#'
#' @description Calculates the complexity index for a set of fluorophores (a spectral panel).
#'
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and 1,
#'     with fluorophores in rows and detectors (channels) in columns. More than one
#'     spectra are needed to calculate anything here.
#' @return The complexity index (condition number) for the
#'     spectral flow cytometry panel.
#' @export

calculate.complexity.index <- function( spectra ) {

  svd.result <- svd( spectra )

  singular.values <- svd.result$d

  complexity.index <- max( singular.values ) / min( singular.values )

  complexity.index <- round( complexity.index, 2 )

  return( complexity.index )
}
