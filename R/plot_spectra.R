# plotting function for spectra

# needs ifelse for case spectra.n = 1


#' @title Plot Fluorophore Spectra
#'
#' @description This function plots the fluorophore spectra, optionally splitting
#'     by excitation lasers, and saves the plots as JPEG files.
#'
#' @importFrom ggplot2 ggplot aes geom_path geom_point labs theme_minimal
#' @importFrom ggplot2 element_text facet_wrap ggsave theme
#' @importFrom tidyr pivot_longer
#' @importFrom utils read.csv
#' @importFrom stats setNames
#'
#' @param spectral.matrix Matrix containing spectral data.
#' @param flow.control List containing flow control information, including
#'     spectral channels.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param plot.title Title for the plot.
#' @param plot.dir Directory to save the plot files.
#' @param split.lasers Logical indicating whether to split the plot by
#'     excitation lasers.
#'
#' @return Saves the plot(s) as JPEG files in the specified directory.
#' @export


plot.spectra <- function( spectral.matrix, flow.control, asp,
                         plot.title = "Fluorophore Spectra",
                         plot.dir, split.lasers = TRUE ){

  colnames( spectral.matrix ) <- flow.control$spectral.channel
  fluor.spectra.plotting <- data.frame( spectral.matrix, check.names = FALSE )
  fluor.spectra.plotting$Fluorophore <- rownames( fluor.spectra.plotting )

  # get excitation laser
  laser.order <- c( "UV", "Violet","Blue", "YellowGreen", "Red" )

  data.path <- system.file( "extdata", "fluorophore_database.csv",
                            package = "AutoSpectral" )
  fluorophore.database <- read.csv( data.path )
  fluorophore.database[ fluorophore.database == "" ] <- NA
  lasers <- setNames( fluorophore.database$excitation.laser, fluorophore.database$fluorophore )

  laser.idx <- match( fluor.spectra.plotting$Fluorophore, names( lasers ) )

  fluor.spectra.plotting$Laser <- lasers[ laser.idx ]
  fluor.spectra.plotting$Laser[ is.na( fluor.spectra.plotting$Laser ) ] <- "Violet"
  fluor.spectra.plotting$Laser <- factor( fluor.spectra.plotting$Laser, levels = laser.order )

  plot.width <- ( ncol( fluor.spectra.plotting ) - 1 ) / 64 * 12
  plot.height <- 5 + round( nrow( fluor.spectra.plotting ) / 8, 0 )

  fluor.spectra.long <- tidyr::pivot_longer( fluor.spectra.plotting, -c( Fluorophore, Laser ),
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

  if ( split.lasers ){
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
