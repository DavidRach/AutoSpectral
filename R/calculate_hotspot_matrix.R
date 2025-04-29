# calculate_hotspot_matrix.r

#' @title Calculate Hotspot Matrix
#' @description Calculates the hotspot (unmixing spread) matrix
#'
#' @importFrom MASS ginv
#' @param spectra Fluorophore spectra as fluorophores x detectors
#' @return The hotspot matrix
#' @export


calculate.hotspot.matrix <- function( spectra ) {

  similarity.matrix <- cosine.similarity( t( spectra ) )

  hotspot.matrix <- sqrt( abs( MASS::ginv( similarity.matrix ) ) )

  rownames( hotspot.matrix ) <- rownames( spectra )
  colnames( hotspot.matrix ) <- rownames( spectra )

  return( hotspot.matrix )
}
