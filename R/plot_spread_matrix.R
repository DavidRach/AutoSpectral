# plot_spread_matrix.r

plot.spread.matrix <- function( spectra, asp, plot.prefix = NULL ) {
 
  if( !is.null( plot.prefix ) ){
    ssm.heatmap.filename <- paste0( plot.prefix, " ",
                                           asp$ssm.heatmap.file.name, 
                                           ".jpg" )
  } else {
    ssm.heatmap.filename <- paste0( asp$ssm.heatmap.file.name, 
                                           ".jpg" )
  }
  
  # calculate SSM
  ssm <- calculate.ssm( spectra )
  
  # plot and save
  ssm.df <- as.data.frame( ssm )
  
  ssm.df <- ssm.df %>%
    mutate( Fluor1 = factor( rownames( ssm.df ), 
                             levels = rownames( ssm.df ) ) ) %>%
    pivot_longer( cols = -Fluor1, names_to = "Fluor2", values_to = "value" ) %>%
    mutate( Fluor2 = factor( Fluor2, levels = rev( colnames( ssm.df ) ) ) )
  
  ssm.heatmap <- ggplot( ssm.df, aes( Fluor1, Fluor2, fill = value ) ) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    coord_fixed( ratio = 1 ) +
    theme( axis.text.x = element_text( angle = 90, hjust = 1 ) ) +
    labs( x = NULL, 
          y = NULL, fill = "Spillover Spread" )
  
  ggsave( file = file.path( asp$figure.similarity.heatmap.dir,
                            ssm.heatmap.filename ),
          plot = ssm.heatmap,
          width = asp$figure.width, height = asp$figure.height )
  
  # save table
  ssm.heatmap.filename <- sub( ".jpg", ".csv", ssm.heatmap.filename )
  write.csv( file = file.path( asp$figure.similarity.heatmap.dir,
                               ssm.heatmap.filename ), ssm )
  
}