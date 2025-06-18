# similarity_matrix_plot.r

#' @title Similarity Matrix Plot
#'
#' @description
#' This function plots the similarity matrix (cosine similarity) as a heatmap
#' and saves it as a JPEG file. It also calculates and displays the complexity
#' index (condition number) of the matrix.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_minimal
#' @importFrom ggplot2 coord_fixed element_text labs ggsave
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr pivot_longer
#'
#' @param spectra Data frame or matrix containing spectral data.
#' @param filename Character string for the output file. Default is
#' `autospectral_similarity_matrix`.
#' @param plot.prefix Optional prefix for the plot filename. Default is `NULL`
#' @param output.dir File path where the plot will be created. Default is
#' `figure_similarity_heatmap`. The directory will be created if it does not
#' already exist.
#' @param figure.width Numeric. Width of output plot. Default is `8`.
#' @param figure.height Numeric. Height of output plot. Default is `6`.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `viridis`. Options
#' are the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#'
#' @return Saves the heatmap plot as a JPEG file in the similarity directory.
#'
#' @export

similarity.matrix.plot <- function( spectra,
                                    filename = "autospectral_similarity_matrix",
                                    plot.prefix = NULL,
                                    output.dir = "figure_similarity_heatmap",
                                    figure.width = 8, figure.height = 6,
                                    color.palette = "viridis" ){

  if ( !is.null( plot.prefix ) ){
    similarity.heatmap.filename <- paste0( plot.prefix, " ", filename, ".jpg" )
  } else {
    similarity.heatmap.filename <- paste0( filename, ".jpg" )
  }

  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

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
    scale_fill_viridis_c( option = color.palette ) +
    theme_minimal() +
    coord_fixed( ratio = 1 ) +
    theme( axis.text.x = element_text( angle = 90, hjust = 1 ) ) +
    labs( x = paste( "Complexity Index", complexity.index ),
          y = NULL, fill = "Cosine Similarity" )

  ggsave( filename = file.path( output.dir, similarity.heatmap.filename ),
          plot = similarity.heatmap,
          width = figure.width, height = figure.height )

}
