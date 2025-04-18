# plot_af_heatmap.r
# AF error heatmaps

plot.af.heatmap <- function( error.matrix, metric,
                            x.var, y.var, x.var.names, y.var.names, asp,
                            plot.prefix = NULL,
                            plot.width = NULL, plot.height = NULL, 
                            number.labels = TRUE ){
  
  if( is.null( plot.prefix ) ){
    plot.title <- paste( metric, "heatmap.jpg" )
  } else {
    plot.title <- paste( plot.prefix, metric, "heatmap.jpg"  )
  }
  
  if( is.null( plot.width ) ){
    plot.width <- length( x.var.names ) * 0.5 + 3
  }
  
  if( is.null( plot.height ) ){
    plot.height <- length( y.var.names ) * 0.5 + 1
  }
  
  error.matrix <- as.data.frame( error.matrix )
  
  error.matrix[[ y.var ]] <- rownames( error.matrix )

  error.matrix$row.sum <- rowSums( error.matrix[ setdiff( names( error.matrix), y.var)])

  error.long <- error.matrix %>%
    pivot_longer(
      cols = -c( all_of( y.var ), "row.sum" ), 
      names_to = x.var, 
      values_to = metric
    ) %>%
    mutate(
      !!sym( y.var ) := factor( .data[[ y.var ]], levels = rev( y.var.names ) ),
      !!sym( x.var ) := factor( .data[[ x.var ]], levels = x.var.names )
    )
  
  x.offset <- length( unique( error.long[[ x.var ]] ) ) + 1

  af.heatmap <- ggplot( error.long, aes( x = !!sym( x.var ), y = !!sym( y.var ), 
                                fill = !!sym( metric ) ) ) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs( title = metric ) +
    theme_minimal() +
    coord_fixed( ratio = 1, clip = "off" ) +
    expand_limits( x = x.offset ) +
    theme(
      axis.text.x = element_text( angle = 90, hjust = 1 ),
      legend.position = "left",
      plot.margin = margin( 10, 30, 10, 10 ), 
      panel.background = element_rect( fill = "transparent", colour = NA ),
      plot.background = element_rect( fill = "transparent", colour = NA ),
      panel.border = element_blank()
    )
  
  if( number.labels ){
    af.heatmap <- af.heatmap + geom_text( aes( label = round( !!sym( metric ), 2)), 
                   color = "white", size = 3 )
  }
  
  af.heatmap <- af.heatmap + geom_text( data = error.long, 
                     aes( x = x.offset, y = !!sym( y.var ), 
                          label = round( row.sum, 2 ) ),
                     size = 3, hjust = 0 )
  
  ggsave( file = file.path( asp$figure.af.dir, plot.title ),
          plot = af.heatmap,
          width = plot.width, height = plot.height,
          limitsize = FALSE )
  
}