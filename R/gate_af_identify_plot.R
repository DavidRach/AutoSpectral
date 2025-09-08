# gate_af_identify_plot.r

#' @title Plot Autofluorescence Identification Gate
#'
#' @description
#' This function plots the sample being used to identify intrusive
#' autofluorescence in the single-stained controls. The input data are expected
#' to be PCA projections of the unstained sample with an accompanying region to
#' identify the low-autofluorescence cell region.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_gradientn theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 guide_colorbar
#' @importFrom scattermore geom_scattermore
#' @importFrom fields interp.surface
#' @importFrom rlang .data
#'
#' @param gate.data Matrix containing autofluorescence data points.
#' @param samp Sample identifier.
#' @param gate.region Dataframe containing region boundary information.
#' @param gate.bound.density Density (e.g., from MASS:kde2d) information
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return Saves the plot as a JPEG file in the specified directory.
#'
#' @export

gate.af.identify.plot <- function( gate.data, samp, gate.region,
                                   gate.bound.density, asp ) {

  gate.data.ggp <- data.frame(
    x = gate.data[ , 1 ],
    y = gate.data[ , 2 ],
    z = interp.surface( gate.bound.density, gate.data )
    )

  density.palette <- get.density.palette( gate.data.ggp$z, asp )

  # get axis labels
  axes.labels <- colnames( gate.data )

  # get data range & step
  x.min <- min( gate.data[ , 1 ] )
  x.max <- max( gate.data[ , 1 ] )

  y.min <- min( gate.data[ , 2 ] )
  y.max <- max( gate.data[ , 2 ] )

  x.breaks <- round( seq( x.min, x.max, length.out = 10 ) )
  y.breaks <- round( seq( y.min, y.max, length.out = 10 ) )

  gate.plot <- ggplot( gate.data.ggp, aes( .data$x, .data$y,
                                         color = .data$z ) ) +
    scale_x_continuous(
      name = axes.labels[ 1 ],
      breaks = x.breaks,
      labels = x.breaks,
      limits = c( x.min, x.max ),
      expand = expansion( asp$af.figure.gate.scale.expand ) ) +
    scale_y_continuous(
      name = axes.labels[ 2 ],
      breaks = y.breaks,
      labels = y.breaks,
      limits = c( y.min, y.max ),
      expand = expansion( asp$af.figure.gate.scale.expand ) ) +
    geom_scattermore( pointsize = 1.2 * asp$figure.gate.point.size,
                      alpha = 1, na.rm = TRUE ) +
    scale_color_gradientn( "", labels = NULL, colors = density.palette,
                           guide = guide_colorbar( barwidth = asp$figure.gate.bar.width,
                                                   barheight = asp$figure.gate.bar.height ) ) +
    theme_bw() +
    theme( plot.margin = margin( asp$figure.margin, asp$figure.margin,
                                 asp$figure.margin, asp$figure.margin ),
           legend.margin = margin( asp$figure.gate.bar.margin,
                                   asp$figure.gate.bar.margin, asp$figure.gate.bar.margin,
                                   asp$figure.gate.bar.margin ),
           legend.position = "none",
           axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
           axis.text = element_text( size = asp$figure.axis.text.size ),
           axis.title = element_text( size = asp$figure.axis.title.size ),
           panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() )

  gate.plot <- gate.plot +
    geom_path( aes( .data$x, .data$y, color = NULL ),
               data = gate.region, linewidth = asp$figure.gate.line.size )

  ggsave( file.path( asp$figure.clean.control.dir,
                     paste( asp$af.plot.define.filename, samp, ".jpg" ) ),
          plot = gate.plot, width = asp$figure.width,
          height = asp$figure.height )

}
