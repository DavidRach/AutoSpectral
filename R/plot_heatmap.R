# plot_heatmap.r


#' @title Plot Heatmap
#'
#' @description This function plots a matrix as a heatmap and saves it as a
#'     JPEG file.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_minimal
#' @importFrom ggplot2 coord_fixed element_text labs ggsave
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr pivot_longer
#'
#' @param matrix Matrix or dataframe containing spectral data.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param plot.prefix Optional prefix for the plot filename.
#' @param number.labels Logical indicating whether to add number labels to the heatmap.
#'
#' @return Saves the heatmap plot as a JPEG file and the SSM data as a CSV file
#'     in the specified directory.
#' @export

plot.heatmap <- function( matrix, asp, plot.prefix = NULL, number.labels ) {

  if( !is.null( plot.prefix ) ){
    heatmap.filename <- paste( plot.prefix, "heatmap.jpg" )
  } else {
    heatmap.filename <- "heatmap.jpg"
  }

  # plot and save
  heatmap.df <- as.data.frame( matrix )

  heatmap.df <- heatmap.df %>%
    mutate( Fluor1 = factor( rownames( heatmap.df ),
                             levels = rownames( heatmap.df ) ) ) %>%
    pivot_longer( cols = -Fluor1, names_to = "Fluor2", values_to = "value" ) %>%
    mutate( Fluor2 = factor( Fluor2, levels = rev( colnames( heatmap.df ) ) ) )

  heatmap.plot <- ggplot( heatmap.df, aes( Fluor1, Fluor2, fill = value ) ) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    coord_fixed( ratio = 1 ) +
    theme( axis.text.x = element_text( angle = 90, hjust = 1 ) ) +
    labs( x = NULL,
          y = NULL, fill = plot.prefix )

  if( number.labels ){
    heatmap.plot <- heatmap.plot + geom_text( aes( label = round( heatmap.df$value, 2 ) ),
                                          color = "white", size = 3 )
  }

  ggsave( filename = file.path( asp$figure.similarity.heatmap.dir,
                                heatmap.filename ),
          plot = heatmap.plot,
          width = asp$figure.width, height = asp$figure.height )

}
