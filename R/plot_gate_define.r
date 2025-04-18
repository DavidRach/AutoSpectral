# plot_gate_define.r
#
# Copyright (c) 2020 VIB (Belgium) & Babraham Institute (United Kingdom)
#
# Software written by Carlos P. Roca, as research funded by the European Union.
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.


# Plots gate, including boundary and region.

plot.gate.define <- function( samp, gate.data, gate.marker, gate.bound,
    gate.region, gate.population, scatter.and.channel.label, asp )
{
  
  gate.data.ggp <- data.frame(
      x = gate.data[ , 1 ],
      y = gate.data[ , 2 ],
      z = interp.surface( gate.bound$density, gate.data ) )
  
  gate.data.ggp$z[ is.na( gate.data.ggp$z ) ] <- min( gate.bound$density$z )
    
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
    
    gate.boundary.ggp <- data.frame(
      x = c( gate.population$boundary$x,
             gate.population$boundary$x[ 1 ] ),
      y = c( gate.population$boundary$y,
             gate.population$boundary$y[ 1 ] )
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
            stroke = 0.1 * asp$figure.gate.point.size, alpha = 1 ) +
        scale_color_gradientn( "", labels = NULL, colors = density.palette,
            guide = guide_colorbar( barwidth = asp$figure.gate.bar.width,
                barheight = asp$figure.gate.bar.height ) ) +
        geom_path( aes( .data$x, .data$y, color = NULL ),
                 data = gate.bound.ggp, linewidth = asp$figure.gate.line.size,
                 linetype = "dashed" ) +
        geom_path( aes( .data$x, .data$y, color = NULL ),
                 data = gate.region.ggp, linewidth = asp$figure.gate.line.size ) +
        geom_path( aes( .data$x, .data$y, color = NULL ),
                 data = gate.boundary.ggp, linewidth = asp$figure.gate.line.size ) +
        geom_point( data = gate.bound$density.max,
                  size = 1.9 * asp$figure.gate.point.size,
                  stroke = 0.1 * asp$figure.gate.point.size,
                  color = asp$gate.tesselation.color ) +
        geom_text( data = gate.bound$density.max,
                 aes( label = .data$num.label ),
                 hjust = 0, vjust = 0, size = asp$figure.axis.text.size / 2.5,
                 color = asp$gate.tesselation.color ) +
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

    ggsave( file.path( asp$figure.gate.dir, sprintf( "%s.jpg", samp ) ),
            plot = gate.plot, width = asp$figure.width,
            height = asp$figure.height )
    
}

