# calculate_similarity_matrix.r

#' @title calculate.similarity.matrix
#'
#' @description Calculates the similarity matrix (set of cosine similarity values)
#'     for a set of fluorophores (a spectral panel).
#'
#' @importFrom lsa cosine
#' @param spectra Spectral signatures of fluorophores, normalized between 0 and 1,
#'     with fluorophores in rows and detectors (channels) in columns. More than one
#'     spectra are needed to calculate anything here.
#' @return The similarity matrix for the spectral flow cytometry panel.
#' @export

calculate.similarity.matrix <- function( spectra ) {

  similarity.matrix <- cosine( t( spectra ) )

  return( similarity.matrix )
}
