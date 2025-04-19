# plot_gate_af_sample.r


#' @title Plot Gate Autofluorescence Sample
#'
#' @description This function plots the gate autofluorescence sample,
#'     including upper and lower boundaries, using ggplot2 and other necessary
#'     packages.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_gradientn theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 guide_colorbar
#' @importFrom scattermore geom_scattermore
#' @importFrom KernSmooth bkde2D dpik
#' @importFrom fields interp.surface
#' @importFrom rlang .data
#'
#' @param samp Sample identifier.
#' @param af.data Matrix containing autofluorescence data points.
#' @param af.boundary.lower Matrix containing lower boundary information.
#' @param af.boundary.upper Matrix containing upper boundary information.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#'
#' @return Saves the plot as a JPEG file in the specified directory.
#' @export



plot.gate.af.sample <- function( samp, af.data,
                                 af.boundary.lower, af.boundary.upper,
                                 asp ) {

  bandwidth.x <- asp$af.plot.bw.factor * dpik( af.data[ , 1 ] )
  bandwidth.y <- asp$af.plot.bw.factor * dpik( af.data[ , 2 ] )


  af.bound.density <- suppressWarnings(
    bkde2D( af.data, bandwidth = c( bandwidth.x, bandwidth.y ),
      gridsize = c( asp$af.plot.density.grid.n, asp$af.plot.density.grid.n ) )
  )

  names( af.bound.density ) <- c( "x", "y", "z" )

  af.data.ggp <- data.frame(
    x = af.data[ , 1 ],
    y = af.data[ , 2 ],
    z = interp.surface( af.bound.density, af.data ) )

  if( !is.null( af.boundary.lower ) ){
    af.boundary.lower.ggp <- data.frame(
      x = c( af.boundary.lower$x,
             af.boundary.lower$x[ 1 ] ),
      y = c( af.boundary.lower$y,
             af.boundary.lower$y[ 1 ] )
    )
  }

  af.boundary.upper.ggp <- data.frame(
    x = c( af.boundary.upper$x,
           af.boundary.upper$x[ 1 ] ),
    y = c( af.boundary.upper$y,
           af.boundary.upper$y[ 1 ] )
  )

  density.palette <- get.density.palette( af.data.ggp$z, asp )

  # get axis labels
  axes.labels <- colnames( af.data )

  # get data range & step
  x.min <- min( af.data[ , 1 ] )
  x.max <- max( af.data[ , 1 ] )

  y.min <- min( af.data[ , 2 ] )
  y.max <- max( af.data[ , 2 ] )

  x.breaks <- round( seq( x.min, x.max, length.out = 10 ) )
  y.breaks <- round( seq( y.min, y.max, length.out = 10 ) )

  gate.plot <- ggplot( af.data.ggp, aes( .data$x, .data$y,
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
                      stroke = 0.1 * asp$figure.gate.point.size, alpha = 1 ) +
    scale_color_gradientn( "", labels = NULL, colors = density.palette,
                           guide = guide_colorbar( barwidth = asp$figure.gate.bar.width,
                                                   barheight = asp$figure.gate.bar.height ) ) +
    theme_bw() +
    theme( plot.margin = margin( asp$figure.margin, asp$figure.margin,
                                 asp$figure.margin, asp$figure.margin ),
           legend.margin = margin( asp$figure.gate.bar.margin,
                                   asp$figure.gate.bar.margin, asp$figure.gate.bar.margin,
                                   asp$figure.gate.bar.margin ),
           axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
           axis.text = element_text( size = asp$figure.axis.text.size ),
           axis.title = element_text( size = asp$figure.axis.title.size ),
           panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() )

  gate.plot <- gate.plot +
    geom_path( aes( .data$x, .data$y, color = NULL ),
               data = af.boundary.upper.ggp, linewidth = asp$figure.gate.line.size )

  if( !is.null( af.boundary.lower ) ){
    gate.plot <- gate.plot +
      geom_path( aes( .data$x, .data$y, color = NULL ),
                 data = af.boundary.lower.ggp, linewidth = asp$figure.gate.line.size )
  }

  ggsave( file.path( asp$figure.clean.control.dir,
                     paste( asp$af.plot.filename, samp, ".jpg" ) ),
          plot = gate.plot, width = asp$figure.width,
          height = asp$figure.height )

}
