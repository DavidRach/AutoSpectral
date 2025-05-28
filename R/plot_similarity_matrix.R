# plot_similarity_matrix.r

#' @title Similarity Matrix Plot
#'
#' @description This function plots the similarity matrix as a heatmap and saves
#'     it as a JPEG file. It also calculates and displays the complexity index.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_minimal
#' @importFrom ggplot2 coord_fixed element_text labs ggsave
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr pivot_longer
#'
#' @param spectra Data frame containing spectral data.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param plot.prefix Optional prefix for the plot filename.
#'
#' @return Saves the heatmap plot as a JPEG file in the specified directory.
#' @export

similarity.matrix.plot <- function( spectra, asp, plot.prefix = NULL ){

  if ( !is.null( plot.prefix ) ){
    similarity.heatmap.filename <- paste0( plot.prefix, " ",
                                           asp$similarity.heatmap.file.name,
                                           ".jpg" )
  } else {
    similarity.heatmap.filename <- paste0( asp$similarity.heatmap.file.name,
                                           ".jpg" )
  }

  # similarity matrix
  similarity.matrix <- cosine.similarity( t( spectra ) )

  complexity.index <- calculate.complexity.index( spectra )

  similarity.df <- as.data.frame( similarity.matrix )

  similarity.df <- similarity.df %>%
    mutate( Fluor1 = factor( rownames( similarity.df ),
                             levels = rownames( similarity.df ) ) ) %>%
    tidyr::pivot_longer( cols = -Fluor1, names_to = "Fluor2", values_to = "value" ) %>%
    mutate( Fluor2 = factor( Fluor2, levels = rev( colnames( similarity.matrix ) ) ) )

  similarity.heatmap <- ggplot( similarity.df, aes( Fluor1, Fluor2, fill = value ) ) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    coord_fixed( ratio = 1 ) +
    theme( axis.text.x = element_text( angle = 90, hjust = 1 ) ) +
    labs( x = paste( "Complexity Index", complexity.index ),
          y = NULL, fill = "Cosine Similarity" )

  ggsave( filename = file.path( asp$figure.similarity.heatmap.dir,
                            similarity.heatmap.filename ),
          plot = similarity.heatmap,
          width = asp$figure.width, height = asp$figure.height )

}
