# create_heatmap.r


#' @title Create Heatmap Plot
#'
#' @description
#' This function plots a matrix as a heatmap and saves it as a JPEG file.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_minimal
#' @importFrom ggplot2 coord_fixed element_text labs ggsave
#' @importFrom dplyr mutate %>% filter
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#'
#' @param matrix Matrix or dataframe containing spectral data.
#' @param number.labels Logical indicating whether to add number labels to
#' the heatmap.
#' @param plot.prefix Optional prefix for the plot filename. Default is `NULL`,
#' in which case the file will just be called `heatmap.jpg`
#' @param legend.label Character string that will appear on the heatmap legend.
#' Default is `heatmap`
#' @param triangular Logical. Plot the lower triangle of the matrix only,
#' diagonal included. Default is `FALSE`.
#' @param output.dir Optional output directory. Default is `NULL`, in which case
#' the working directory will be used.
#' @param fixed.scale Logical, determines whether to use an externally supplied
#' fixed scale (min and max) for the heatmap color scale. Useful for putting
#' multiple plots on the same scale. Default is `FALSE`
#' @param scale.min Optional numeric. Minimum for the fixed color scale.
#' Default is `NULL`, for no fixed scaling.
#' @param scale.max Optional numeric. Maximum for the fixed color scale.
#' Default is `NULL`, for no fixed scaling.
#' @param color.palette Optional character string defining the viridis color
#' palette to be used for the fluorophore traces. Default is `viridis`. Options
#' are the viridis color options: `magma`, `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#' @param figure.width Numeric. Width of the heatmap figure. Default is `8`.
#' @param figure.height Numeric. Height of the heatmap figure. Default is `6`.
#'
#' @return Saves the heatmap plot as a JPEG file and the SSM data as a CSV file
#' in the specified directory.
#'
#' @export

create.heatmap <- function( matrix, number.labels, plot.prefix = NULL,
                            legend.label = "heatmap",
                            triangular = FALSE,
                            output.dir = NULL,
                            fixed.scale = FALSE,
                            scale.min = NULL, scale.max = NULL,
                            color.palette = "viridis",
                            figure.width = 8, figure.height = 6 ) {

  if ( !is.null( plot.prefix ) ) {
    heatmap.filename <- paste( plot.prefix, "heatmap.jpg" )
  } else {
    heatmap.filename <- "heatmap.jpg"
  }

  if ( is.null( output.dir ) )
    output.dir <- getwd()

  # rearrange data
  heatmap.df <- as.data.frame( matrix )
  row.levels <- rownames( heatmap.df )
  col.levels <- colnames( heatmap.df )
  heatmap.df$Fluor1 <- row.levels

  heatmap.long <- heatmap.df %>%
    pivot_longer( cols = -Fluor1, names_to = "Fluor2", values_to = "value" ) %>%
    mutate( Fluor1 = factor( Fluor1, levels = row.levels ),
            Fluor2 = factor( Fluor2, levels = col.levels ) )


  if ( triangular ){
    heatmap.long <- heatmap.long %>%
      dplyr::filter( as.integer( Fluor1 ) >= as.integer( Fluor2 ) ) %>%
      mutate( Fluor1 = factor( Fluor1, levels = rev( row.levels ) ) )
  } else {
    heatmap.long <- heatmap.long %>%
      mutate( Fluor1 = factor( Fluor1, levels = rev( row.levels ) ) )
  }

  # plot and save
  heatmap.plot <- ggplot( heatmap.long, aes( Fluor1, Fluor2, fill = value ) ) +
    geom_tile() +
    theme_classic() +
    coord_fixed( ratio = 1) +
    theme( axis.text.x = element_text( angle = 90, hjust = 1 ) ) +
    labs( x = NULL, y = NULL, fill = legend.label )

  if ( fixed.scale ) {
    heatmap.plot <- heatmap.plot +
      scale_fill_viridis_c( option = color.palette,
                            limits = c( scale.min, scale.max ) )
  } else {
    heatmap.plot <- heatmap.plot +
      scale_fill_viridis_c( option = color.palette )
  }

  if ( number.labels ) {
    heatmap.plot <- heatmap.plot +
      geom_text( aes( label = round( .data$value, 2 ) ),
                 color = "white", size = 3 )
  }

  ggsave( filename = file.path( output.dir, heatmap.filename ),
          plot = heatmap.plot,
          width = figure.width, height = figure.height )

}
