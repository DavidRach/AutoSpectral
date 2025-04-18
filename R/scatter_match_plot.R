# plotting function for scatter-matching of universal negative

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
  
  
  scatter.density <- bkde2D(
    scatter.plot.data,
    bandwidth = c( bandwidth.x, bandwidth.y ),
    gridsize = c( asp$gate.bound.density.grid.n.cells,
                  asp$gate.bound.density.grid.n.cells ) )
  
  names( scatter.density ) <- c( "x", "y", "z" )
  
  data.ggp <- data.frame(
    x = scatter.plot.data[ , 1 ],
    y = scatter.plot.data[ , 2 ],
    z = interp.surface( scatter.density, scatter.plot.data ),
    group = scatter.plot.data$group )
  
  density.palette <- get.density.palette( data.ggp$z, asp )
  
  ggplot( data.ggp, aes( .data$x, .data$y, color = .data$z, group ) ) +
    geom_scattermore( pointsize = asp$figure.gate.point.size,
                      stroke = asp$figure.gate.point.size, alpha = 1 ) +
    
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
           strip.text = element_text( size = 15, face = "bold" ),
           axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
           axis.text = element_text( size = asp$figure.axis.text.size ),
           axis.title = element_text( size = asp$figure.axis.title.size ),
           panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() )
  
  
  scatter.plot.filename <- paste( fluor.name, "universal negative scatter plot.jpg" )
  
  ggsave( scatter.plot.filename, path = asp$figure.scatter.dir.base,
          width = 12, height = 6 )
  
}