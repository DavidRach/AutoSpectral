# gate_sample_plot.r

#' @title Plot Pre-defined Gate on Sample
#'
#' @description
#' This function plots a pre-defined gate on a sample, using ggplot2 and other
#' necessary packages.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_gradientn theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 guide_colorbar geom_text geom_point
#' @importFrom scattermore geom_scattermore
#' @importFrom MASS kde2d bandwidth.nrd
#' @importFrom fields interp.surface
#' @importFrom rlang .data
#'
#' @param samp Sample identifier.
#' @param gate.data Matrix containing gate data points.
#' @param gate.marker Vector containing gate marker names.
#' @param gate.boundary List containing gate boundary information.
#' @param scatter.and.channel.label Named vector mapping scatter and
#' channel labels.
#' @param control.type Type of control: `beads` or `cells`
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return Saves the plot as a JPEG file in the specified directory.
#'
#' @export

gate.sample.plot <- function( samp, gate.data, gate.marker, gate.boundary,
                              scatter.and.channel.label, control.type, asp )
{

  if ( control.type == "beads" ){

    gate.bound.density.bw.factor <- asp$gate.bound.density.bw.factor.beads
    gate.bound.density.grid.n <- asp$plot.gate.factor * asp$gate.bound.density.grid.n.beads

  } else {

    gate.bound.density.bw.factor <- asp$gate.bound.density.bw.factor.cells
    gate.bound.density.grid.n <- asp$plot.gate.factor * asp$gate.bound.density.grid.n.cells

  }

  gate.bound.density <- MASS::kde2d( gate.data[ , 1 ], gate.data[ , 2 ],
                                     gate.bound.density.bw.factor *
                                       apply( gate.data, 2, bandwidth.nrd ),
                                     n = gate.bound.density.grid.n )

  gate.data.ggp <- data.frame(
    x = gate.data[ , 1 ],
    y = gate.data[ , 2 ],
    z = interp.surface( gate.bound.density, gate.data ) )

  gate.boundary$x[ gate.boundary$x > asp$scatter.data.max.x ] <- asp$scatter.data.max.x
  gate.boundary$y[ gate.boundary$y > asp$scatter.data.max.y ] <- asp$scatter.data.max.y

  gate.boundary.ggp <- data.frame(
    x = c( gate.boundary$x,
           gate.boundary$x[ 1 ] ),
    y = c( gate.boundary$y,
           gate.boundary$y[ 1 ] )
  )

  density.palette <- get.density.palette( gate.data.ggp$z, asp )

  x.lab.idx <- which( scatter.and.channel.label == gate.marker[ 1 ] )
  x.lab <- names( scatter.and.channel.label[ x.lab.idx ] )
  y.lab.idx <- which( scatter.and.channel.label == gate.marker[ 2 ] )
  y.lab <- names( scatter.and.channel.label[ y.lab.idx ] )

  gate.plot <- ggplot( gate.data.ggp, aes( .data$x, .data$y,
                                           color = .data$z ) ) +
    scale_x_continuous(
      name = x.lab,
      breaks = seq( asp$scatter.data.min.x,
                    asp$scatter.data.max.x, asp$data.step ),
      labels = paste0( round( seq( asp$scatter.data.min.x, asp$scatter.data.max.x,
                                   asp$data.step ) / 1e6, 1 ), "e6" ),
      limits = c( asp$scatter.data.min.x,
                  asp$scatter.data.max.x ),
      expand = expansion( asp$figure.gate.scale.expand ) ) +
    scale_y_continuous(
      name = y.lab,
      breaks = seq( asp$scatter.data.min.y,
                    asp$scatter.data.max.y, asp$data.step ),
      labels = paste0( round( seq( asp$scatter.data.min.y, asp$scatter.data.max.y,
                                   asp$data.step ) / 1e6, 1 ), "e6" ),
      limits = c( asp$scatter.data.min.y,
                  asp$scatter.data.max.y ),
      expand = expansion( asp$figure.gate.scale.expand ) ) +
    geom_scattermore( pointsize = asp$figure.gate.point.size,
                      stroke = 0.1 * asp$figure.gate.point.size, alpha = 1, na.rm = TRUE ) +
    scale_color_gradientn( "", labels = NULL, colors = density.palette,
                           guide = guide_colorbar( barwidth = asp$figure.gate.bar.width,
                                                   barheight = asp$figure.gate.bar.height ) ) +
    geom_path( aes( .data$x, .data$y, color = NULL ),
               data = gate.boundary.ggp, linewidth = asp$figure.gate.line.size ) +
    theme_bw() +
    theme( plot.margin = margin( asp$figure.margin, asp$figure.margin,
                                 asp$figure.margin, asp$figure.margin ),
           legend.position = "none",
           axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
           axis.text = element_text( size = asp$figure.axis.text.size ),
           axis.text.x = element_text( angle = 45, hjust = 1 ),
           axis.title = element_text( size = asp$figure.axis.title.size ),
           panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() )

  ggsave( file.path( asp$figure.gate.dir, sprintf( "%s.jpg", samp ) ),
          plot = gate.plot, width = asp$figure.width,
          height = asp$figure.height )

}
