# cosine_similarity.r

#' @title Calculate cosine similarity
#' @description
#' Calculates the cosine similarity between columns of the input matrix (spectra).
#'
#' @param spectra.t The transposed spectral matrix (or dataframe), represented as
#'     detectors in rows and fluorophores in columns.
#' @return The cosine similarity matrix in fluorophores x fluorophores.
#' @export

cosine.similarity <- function( spectra.t ) {
  spectra.t <- as.matrix( spectra.t )
  euclidean.norm <- sqrt( colSums( spectra.t^2 ) )
  dot.product <- t( spectra.t ) %*% spectra.t
  similarity.matrix <- dot.product / outer( euclidean.norm, euclidean.norm )
  return( similarity.matrix )
}
