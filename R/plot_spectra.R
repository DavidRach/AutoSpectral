# plotting function for spectra

# needs ifelse for case spectra.n = 1

plot.spectra <- function( spectral.matrix, flow.control,
                         plot.title = "Fluorophore Spectra",
                         plot.dir, split.lasers = TRUE ){
  
  colnames( spectral.matrix ) <- flow.control$spectral.channel
  fluor.spectra.plotting <- data.frame( spectral.matrix, check.names = FALSE )
  fluor.spectra.plotting$Fluorophore <- rownames( fluor.spectra.plotting )
  
  
  # get excitation laser
  laser.order <- c( "UV", "Violet","Blue", "YellowGreen", "Red" )
  
  fluorophore.database <- read.csv( file.path( asp$database.dir, "fluorophore_database.csv" ) )
  fluorophore.database[ fluorophore.database == "" ] <- NA
  lasers <- setNames( fluorophore.database$excitation.laser, fluorophore.database$fluorophore )
  
  laser.idx <- match( fluor.spectra.plotting$Fluorophore, names( lasers ) )
  
  fluor.spectra.plotting$Laser <- lasers[ laser.idx ]
  fluor.spectra.plotting$Laser[ is.na( fluor.spectra.plotting$Laser ) ] <- "Violet"
  fluor.spectra.plotting$Laser <- factor( fluor.spectra.plotting$Laser, levels = laser.order )
  
  plot.width <- ( ncol( fluor.spectra.plotting ) - 1 ) / 64 * 12
  plot.height <- 5 + round( nrow( fluor.spectra.plotting ) / 8, 0 )
  
  fluor.spectra.long <- pivot_longer( fluor.spectra.plotting, -c( Fluorophore, Laser ),
                                      names_to = "Detector", 
                                      values_to = "Intensity" )
  
  fluor.spectra.long$Detector <-  factor( fluor.spectra.long$Detector, 
                                         levels = unique( fluor.spectra.long$Detector ),
                                         ordered = TRUE )
  
  spectra.plot <- ggplot( fluor.spectra.long, 
                         aes( x = Detector, y = Intensity, 
                             group = Fluorophore, color = Fluorophore ) ) +
    geom_path( linewidth = asp$figure.spectra.line.size ) +
    geom_point( size = asp$figure.spectra.point.size ) +
    labs( title = plot.title,
         x = "Detector",
         y = "Normalized Intensity" ) +
    theme_minimal() +
    theme( axis.text.x = element_text( angle = 45, hjust = 1 )  ) +
    theme( legend.position = "bottom" )
  
  ggsave( file.path( plot.dir, sprintf( "%s.jpg", plot.title )), 
          spectra.plot,
          width = plot.width, height = plot.height,
          limitsize = FALSE )
  
  if( split.lasers ){
    # get number of lasers used
    laser.n <- length( unique( fluor.spectra.plotting$Laser ) )
    
    plot.height <- ( plot.height - 1 ) * laser.n
    
    spectra.plot.split <- ggplot( fluor.spectra.long, 
                                  aes( x = Detector, y = Intensity, 
                                       group = Fluorophore, color = Fluorophore ) ) +
      geom_path( linewidth = asp$figure.spectra.line.size ) +
      geom_point( size = asp$figure.spectra.point.size ) +
      facet_wrap( ~ Laser, nrow = laser.n ) +
      labs( title = plot.title,
            x = "Detector",
            y = "Normalized Intensity" ) +
      theme_minimal() +
      theme( axis.text.x = element_text( angle = 45, hjust = 1 )  ) +
      theme( legend.position = "bottom" )
    
    ggsave( file.path( plot.dir, sprintf( "%s by laser.jpg", plot.title )), 
            spectra.plot.split,
            width = plot.width, height = plot.height,
            limitsize = FALSE )
  }
  
  
}