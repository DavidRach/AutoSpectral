# plot_similarity_matrix.r

plot.similarity.matrix <- function( spectra, asp, plot.prefix = NULL ){
  
  if( !is.null( plot.prefix ) ){
    similarity.heatmap.filename <- paste0( plot.prefix, " ",
                                           asp$similarity.heatmap.file.name, 
                                           ".jpg" )
  } else {
    similarity.heatmap.filename <- paste0( asp$similarity.heatmap.file.name, 
                                           ".jpg" )
  }
  
  # similarity matrix
  similarity.matrix <- cosine( t( spectra ) )
  
  complexity.index <- sum( similarity.matrix[ lower.tri( similarity.matrix ) ] )
  complexity.index <- round( complexity.index, 2 )
  
  similarity.df <- as.data.frame( similarity.matrix )
  
  similarity.df <- similarity.df %>%
    mutate( Fluor1 = factor( rownames( similarity.df ), 
                             levels = rownames( similarity.df ) ) ) %>%
    pivot_longer( cols = -Fluor1, names_to = "Fluor2", values_to = "value" ) %>%
    mutate( Fluor2 = factor( Fluor2, levels = rev( colnames( similarity.matrix ) ) ) )
  
  similarity.heatmap <- ggplot( similarity.df, aes( Fluor1, Fluor2, fill = value ) ) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    coord_fixed( ratio = 1 ) +
    theme( axis.text.x = element_text( angle = 90, hjust = 1 ) ) +
    labs( x = paste( "Complexity Index", complexity.index ), 
          y = NULL, fill = "Cosine Similarity" )
  
  ggsave( file = file.path( asp$figure.similarity.heatmap.dir,
                            similarity.heatmap.filename ),
          plot = similarity.heatmap,
          width = asp$figure.width, height = asp$figure.height )
  
}