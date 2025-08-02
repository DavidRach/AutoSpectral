# cosine_similarity.r

#' @title Calculate cosine similarity
#' @description
#' Calculates the cosine similarity between columns of the input matrix (spectra).
#'
#' @param spectra The spectral matrix (or dataframe), represented as
#' fluorophores in rows and detectors in columns.
#'
#' @return The cosine similarity matrix in fluorophores x fluorophores.
#' @export

cosine.similarity <- function( spectra ) {
  spectra.t <- as.matrix( t( spectra ) )
  euclidean.norm <- sqrt( colSums( spectra.t^2 ) + 1e-9 )
  dot.product <- t( spectra.t ) %*% spectra.t
  similarity.matrix <- dot.product / outer( euclidean.norm, euclidean.norm )
  return( similarity.matrix )
}
