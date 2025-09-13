# gate_af_plot.r

#' @title Plot Autofluorescence Gate
#'
#' @description
#' This function plots the AF gate, including intermediate steps, using ggplot2
#' and other necessary packages.
#'
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_color_gradientn theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin expansion ggsave
#' @importFrom ggplot2 guide_colorbar
#' @importFrom scattermore geom_scattermore
#' @importFrom fields interp.surface
#' @importFrom rlang .data
#'
#' @param samp Sample identifier.
#' @param gate.data Data frame containing gate data points.
#' @param gate.bound List containing gate boundary information.
#' @param gate.region List containing gate region information.
#' @param gate.population List containing gate population information.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return Saves the plot as a JPEG file in the specified directory.
#' @export

gate.af.plot <- function( samp, gate.data, gate.bound, gate.region,
                          gate.population, asp ) {

    gate.data.ggp <- data.frame(
          x = gate.data[ , 1 ],
          y = gate.data[ , 2 ],
          z = interp.surface( gate.bound$density, gate.data ) )

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
            stroke = 0.1 * asp$figure.gate.point.size, alpha = 1, na.rm = TRUE ) +
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
            axis.text.x = element_text( angle = 45, hjust = 1 ),
            axis.title = element_text( size = asp$figure.axis.title.size ),
            panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank() )

    gate.bound.ggp <- data.frame(
            x = c(
                gate.bound$x.low,
                gate.bound$x.high,
                gate.bound$x.high,
                gate.bound$x.low,
                gate.bound$x.low
            ),
            y = c(
                gate.bound$y.low,
                gate.bound$y.low,
                gate.bound$y.high,
                gate.bound$y.high,
                gate.bound$y.low
            )
        )

     gate.region.ggp <- data.frame(
            x = c(
                gate.region$x.low,
                gate.region$x.high,
                gate.region$x.high,
                gate.region$x.low,
                gate.region$x.low
            ),
            y = c(
                gate.region$y.low,
                gate.region$y.low,
                gate.region$y.high,
                gate.region$y.high,
                gate.region$y.low )
        )

    gate.plot <- gate.plot +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                data = gate.region.ggp, linewidth = asp$figure.gate.line.size )

   gate.boundary.ggp <- data.frame(
            x = c( gate.population$boundary$x,
                gate.population$boundary$x[ 1 ] ),
            y = c( gate.population$boundary$y,
                gate.population$boundary$y[ 1 ] )
        )

    gate.plot <- gate.plot +
            geom_path( aes( .data$x, .data$y, color = NULL ),
                data = gate.boundary.ggp, linewidth = asp$figure.gate.line.size )

    ggsave( file.path( asp$figure.clean.control.dir,
                       paste( samp, asp$af.plot.define.filename, ".jpg" ) ),
            plot = gate.plot, width = asp$figure.width,
            height = asp$figure.height )

}

