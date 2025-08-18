# spectral_ribbon_plot.r

#' @title Spectral Ribbon Plot
#'
#' @description
#' This function generates spectral ribbon plots for positive and negative
#' expression data.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes scale_y_continuous geom_bin2d facet_wrap xlab
#' @importFrom ggplot2 ylab scale_fill_gradientn theme_minimal theme element_text
#' @importFrom ggplot2 element_blank ggsave
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom scales trans_new
#'
#' @param pos.expr.data A matrix containing the positive expression data.
#' @param neg.expr.data A matrix containing the negative expression data.
#' @param spectral.channel A character vector specifying the spectral channels.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param fluor.name A character string specifying the fluorophore name.
#' @param title A character string to prefix the plot file name.
#' Default is `NULL`
#' @param af A logical value indicating whether autofluorescence removal is
#' being performed. Default is `FALSE`
#' @param removed.data A matrix containing the removed data, if applicable.
#' Default is `NULL`. If omitted, only two groups (facets) are plotted.
#' @param figure.dir Output folder where the figures will be created. Default is
#' `NULL`, enabling automatic selection inside AutoSpectral.
#' @param factor.names Optional titles for the facets on the plot. Default is
#' `NULL`, enabling automatic selection inside AutoSpectral. A character vector
#' containing three elements should be provided unless `af` is set to `TRUE` and
#' `removed.data` is `NULL`, in which case two labels must be provided.
#'
#' @return None. The function saves the generated spectral ribbon plot to
#' a file.
#'
#' @export

spectral.ribbon.plot <- function( pos.expr.data, neg.expr.data,
                                  spectral.channel, asp, fluor.name,
                                  title = NULL, af = FALSE,
                                  removed.data = NULL, figure.dir = NULL,
                                  factor.names = NULL ){

  if ( !af ) {

    if ( is.null( title ) )
      title <- "Scatter match"

    if ( is.null( figure.dir ) )
      figure.dir <- asp$figure.spectral.ribbon.dir

    neg.mfi <- apply( neg.expr.data[ , spectral.channel ], 2, median )

    pos.background.subtracted <- data.frame( pos.expr.data[ , spectral.channel ],
                                             check.names = FALSE )
    pos.background.subtracted <- sweep( pos.background.subtracted, 2, neg.mfi,
                                        FUN = "-" )

    pos.data.plot <- data.frame( pos.expr.data[ , spectral.channel ],
                                 check.names = FALSE )

    neg.data  <- data.frame( neg.expr.data[ , spectral.channel ],
                             check.names = FALSE )

    if ( is.null( factor.names ) ) {
      pos.background.subtracted$group <- fluor.name
      pos.data.plot$group <- paste( "Raw", fluor.name )
      neg.data$group <- "Negative"

      ribbon.plot.data <- rbind( pos.background.subtracted, pos.data.plot, neg.data )

      ribbon.plot.data$group <- factor( ribbon.plot.data$group,
                                        levels = c( "Negative",
                                                    paste( "Raw", fluor.name ),
                                                    fluor.name ) )
    } else {

      if ( length( factor.names ) != 3 )
        stop( "Three labels must be provided via `factor.names` if used." )

      pos.background.subtracted$group <- factor.names[ 1 ]
      pos.data.plot$group <- factor.names[ 2 ]
      neg.data$group <- factor.names[ 3 ]

      ribbon.plot.data <- rbind( pos.background.subtracted, pos.data.plot, neg.data )

      ribbon.plot.data$group <- factor( ribbon.plot.data$group,
                                        levels = factor.names )
    }

  } else {

    if ( is.null( title ) )
      title <- "AF removal"

    if ( is.null( figure.dir ) )
      figure.dir <- asp$figure.clean.control.dir

    original.data <- data.frame( pos.expr.data[ , spectral.channel ],
                                 check.names = FALSE )

    cleaned.data <- data.frame( neg.expr.data[ , spectral.channel ],
                                check.names = FALSE )

    if ( !is.null( removed.data ) )
      removed.data <- data.frame( removed.data[ , spectral.channel ],
                                  check.names = FALSE )

    if ( is.null( factor.names ) ) {
      original.data$group <- paste( "Original", fluor.name )
      cleaned.data$group <- paste( "Cleaned", fluor.name )
      removed.data$group <- "Removed events"

      ribbon.plot.data <- rbind( original.data, cleaned.data, removed.data )

      ribbon.plot.data$group <- factor( ribbon.plot.data$group,
                                        levels = c( paste( "Original", fluor.name ),
                                                    paste( "Cleaned", fluor.name ),
                                                    "Removed events" ) )
    } else {

      if ( !is.null( removed.data ) & length( factor.names ) != 3 )
        stop( "Three labels must be provided via `factor.names` if used with 3 groups." )
      if ( is.null( removed.data ) & length( factor.names ) != 2 )
        stop( "Two labels must be provided via `factor.names` if used with 2 groups." )

      original.data$group <- factor.names[ 1 ]
      cleaned.data$group <- factor.names[ 2 ]

      if ( !is.null( removed.data ) ) {
        removed.data$group <- factor.names[ 3 ]
        ribbon.plot.data <- rbind( original.data, cleaned.data, removed.data )
      } else {
        ribbon.plot.data <- rbind( original.data, cleaned.data )
      }

      ribbon.plot.data$group <- factor( ribbon.plot.data$group,
                                        levels = factor.names )
    }

  }

  ribbon.plot.long <- tidyr::pivot_longer( ribbon.plot.data,
                                           cols = -group,
                                           names_to = "channel",
                                           values_to = "value" )

  ribbon.breaks <- asp$ribbon.breaks
  ribbon.labels <- sapply( ribbon.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  ribbon.limits <- c( asp$ribbon.plot.min, asp$expr.data.max )

  # biexponential transform for scale
  biexp.transform <- flowjo_biexp( channelRange = asp$default.transformation.param$length,
                                  maxValue = asp$default.transformation.param$max.range,
                                  pos = asp$default.transformation.param$pos,
                                  neg = asp$default.transformation.param$neg,
                                  widthBasis = asp$default.transformation.param$width,
                                  inverse = FALSE )

  biexp.inverse <- flowjo_biexp( channelRange = asp$default.transformation.param$length,
                                maxValue = asp$default.transformation.param$max.range,
                                pos = asp$default.transformation.param$pos,
                                neg = asp$default.transformation.param$neg,
                                widthBasis = asp$default.transformation.param$width,
                                inverse = TRUE )

  plot.biexp.transform <- trans_new(
    name = "biexp",
    transform = biexp.transform,
    inverse = biexp.inverse
  )

  ribbon.plot <- suppressWarnings( ggplot( ribbon.plot.long,
          aes( factor( channel, levels = unique( channel ) ), value ) ) +
    scale_y_continuous( trans = plot.biexp.transform,
                        breaks = ribbon.breaks,
                        limits = ribbon.limits,
                        labels = ribbon.labels ) +
    geom_bin2d( bins = asp$ribbon.bins ) +
    facet_wrap( ~ group, ncol = 1 ) +
    xlab( "Detector" ) +
    ylab( "Intensity" ) +
    scale_fill_gradientn( colours = asp$density.palette.base.color,
                          values = asp$ribbon.scale.values ) +
    theme_minimal() +
    theme( axis.text.x = element_text( angle = asp$ribbon.plot.axis.text.angle,
                                       vjust = 1, hjust = 1 ),
           panel.grid.minor = element_blank(),
           legend.position = "none",
           strip.text = element_text( size = asp$ribbon.plot.strip.text.size,
                                      face = asp$ribbon.plot.strip.text.face ) )
  )

  ribbon.plot.filename <- paste( title, fluor.name, asp$ribbon.plot.filename )

  suppressWarnings(
    ggsave( ribbon.plot.filename, plot = ribbon.plot,
            path = figure.dir,
            width = asp$ribbon.plot.width, height = asp$ribbon.plot.height )
  )
}
