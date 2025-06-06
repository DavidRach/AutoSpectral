# scatter_match_plot.r
# plotting function for scatter-matching of universal negative

#' @title Plot Scatter-Matching of Universal Negative
#'
#' @description This function generates scatter plots for matching positive
#' and negative expression data based on specified scatter parameters.
#'
#' @importFrom KernSmooth bkde2D dpik
#' @importFrom ggplot2 ggplot aes scale_color_gradientn facet_wrap
#' @importFrom ggplot2 xlab ylab scale_x_continuous scale_y_continuous theme_bw
#' @importFrom ggplot2 theme element_rect element_text margin ggsave guide_colorbar
#' @importFrom rlang .data
#' @importFrom scattermore geom_scattermore
#'
#' @param pos.expr.data A matrix containing the positive expression data.
#' @param neg.expr.data A matrix containing the negative expression data.
#' @param fluor.name A character string specifying the fluorophore name.
#' @param scatter.param A character vector specifying the scatter parameters.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return None. The function saves the generated scatter plot to a file.
#'
#' @export


scatter.match.plot <- function( pos.expr.data, neg.expr.data, fluor.name,
                                scatter.param, asp ){

  pos.scatter.plot <- data.frame( pos.expr.data[ , scatter.param ],
                                  check.names = FALSE )

  neg.scatter.plot  <- data.frame( neg.expr.data[ , scatter.param ],
                                   check.names = FALSE )

  pos.scatter.plot$group <- fluor.name
  neg.scatter.plot$group <- "Negative"

  scatter.plot.data <- rbind( pos.scatter.plot, neg.scatter.plot )

  bandwidth.x <- asp$gate.bound.density.bw.factor.cells * dpik( scatter.plot.data[ , 1 ] )
  bandwidth.y <- asp$gate.bound.density.bw.factor.cells * dpik( scatter.plot.data[ , 2 ] )

  scatter.density <- suppressMessages( suppressWarnings(
    bkde2D( scatter.plot.data, bandwidth = c( bandwidth.x, bandwidth.y ),
    gridsize = c( asp$gate.bound.density.grid.n.cells,
                  asp$gate.bound.density.grid.n.cells ) ) ) )

  names( scatter.density ) <- c( "x", "y", "z" )

  data.ggp <- data.frame(
    x = scatter.plot.data[ , 1 ],
    y = scatter.plot.data[ , 2 ],
    z = interp.surface( scatter.density, scatter.plot.data ),
    group = scatter.plot.data$group )

  density.palette <- get.density.palette( data.ggp$z, asp )

  ggplot( data.ggp, aes( .data$x, .data$y, color = .data$z, group ) ) +
    geom_scattermore( pointsize = asp$figure.gate.point.size,
                      stroke = asp$figure.gate.point.size, alpha = 1, na.rm = TRUE ) +

    scale_color_gradientn( "", labels = NULL, colors = density.palette,
          guide = guide_colorbar( barwidth = asp$figure.gate.bar.width,
          barheight = asp$figure.gate.bar.height ) ) +
    facet_wrap( ~ group, ncol = 2 ) +
    xlab( scatter.param[ 1 ] ) +
    ylab( scatter.param[ 2 ] ) +
    scale_x_continuous( breaks = seq( asp$scatter.data.min.x,
                                      asp$scatter.data.max.x, asp$data.step ),
                        labels = paste0( round( seq( asp$scatter.data.min.x,
                                                     asp$scatter.data.max.x,
                                                     asp$data.step ) / 1e6, 1 ),
                                         "e6" ),
                        limits = c( asp$scatter.data.min.x,
                                    asp$scatter.data.max.x ) ) +
    scale_y_continuous( breaks = seq( asp$scatter.data.min.y,
                                      asp$scatter.data.max.y, asp$data.step ),
                        labels = paste0( round( seq( asp$scatter.data.min.y,
                                                     asp$scatter.data.max.y,
                                                     asp$data.step ) / 1e6, 1 ),
                                         "e6" ),
                        limits = c( asp$scatter.data.min.y,
                                    asp$scatter.data.max.y ) ) +
    theme_bw() +
    theme( plot.margin = margin( asp$figure.margin, asp$figure.margin,
                                 asp$figure.margin, asp$figure.margin ),
           legend.position = "none",
           strip.background = element_rect( fill = "white" ),
           strip.text = element_text( size = asp$scatter.match.plot.text.size,
                                      face = asp$scatter.match.plot.text.face ),
           axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
           axis.text = element_text( size = asp$figure.axis.text.size ),
           axis.title = element_text( size = asp$figure.axis.title.size ),
           panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() )

  scatter.plot.filename <- paste( fluor.name, asp$scatter.match.plot.filename )

  ggsave( scatter.plot.filename, path = asp$figure.scatter.dir.base,
          width = asp$scatter.match.plot.width, height = asp$scatter.match.plot.height )

}
