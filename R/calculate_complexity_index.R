# calculate_complexity_index.r

#' @title calculate.complexity.index
#'
#' @description Calculates the complexity index for a set of fluorophores (a spectral panel).
#'
#' @importFrom lsa cosine
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and 1,
#' with fluorophores in rows and detectors (channels) in columns. More than one
#' spectra are needed to calculate anything here.
#' @return The complexity index (condition number) for the
#' spectral flow cytometry panel.
#' @export

calculate.complexity.index <- function( spectra ) {

  similarity.matrix <- cosine( t( spectra ) )

  complexity.index <- sum( similarity.matrix[ lower.tri( similarity.matrix ) ] )
  complexity.index <- round( complexity.index, 2 )

  return( complexity.index )
}
