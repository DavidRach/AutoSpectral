# refine_unmixing_log_plot_together.r

#' @title Plot Density of Log Values Together
#'
#' @description
#' This function plots the density of log values for initial and final steps
#' of AutoSpectral, and calculation with positive and negative populations,
#' using ggplot2.
#'
#' @importFrom ggplot2 ggplot aes geom_density scale_color_manual scale_fill_manual
#' @importFrom ggplot2 scale_linetype_manual labs theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin ggsave
#' @importFrom dplyr mutate
#' @importFrom rlang sym .data
#'
#' @param x.table List of data frames containing values to be plotted.
#' @param x.label Label for the x-axis.
#' @param plot.file.path Path to save the plot file.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return Saves the plot as a JPEG file in the specified directory.

refine.unmixing.log.plot.together <- function( x.table, x.label, plot.file.path,
                                               asp ) {

  # get values avoiding matrix diagonal
  x.data <- lapply( x.table, function( xt ) {
    xt.n <- nrow( xt )
    xt[ -( 1 + 0:(xt.n-1) * (xt.n+1) ) ]
  } )
  x.data.n <- length( x.data )

  # remove other zero values
  x.data <- lapply( x.data, function( xd )
    xd[ xd != 0 ]
  )

  x.ggdata <- lapply( 1 : x.data.n, function( xd.idx ) {
    data.frame( log.abs.x = log10( abs( x.data[[ xd.idx ]] ) ),
                sign = sign( x.data[[ xd.idx ]] ) * xd.idx )
  } )
  x.ggdata <- do.call( rbind, x.ggdata )
  x.ggdata$sign <- factor( x.ggdata$sign,
                           levels = c( 0, rep( 1:x.data.n, each = 2 ) * c( -1, 1 ) ) )

  scale.color.value <- c(
    "-1" = asp$density.color.initial, "1" = asp$density.color.initial,
    "-2" = asp$density.color.final, "2" = asp$density.color.final,
    "-3" = asp$density.color.posnegpop, "3" = asp$density.color.posnegpop
  )

  scale.linetype.value <- c(
    "-1" = "dashed", "1" = "solid",
    "-2" = "dashed", "2" = "solid",
    "-3" = "dashed" , "3"= "solid"
  )

  x.ggplot <- ggplot( x.ggdata, aes( x = .data$log.abs.x, color = .data$sign,
                                     fill = .data$sign, linetype = .data$sign ) ) +
    geom_density( alpha = 0.1, linewidth = asp$figure.density.line.size ) +
    scale_color_manual( values = scale.color.value ) +
    scale_fill_manual( values = scale.color.value  ) +
    scale_linetype_manual( values = scale.linetype.value ) +
    labs( x = sprintf( "log10( Absolute %s )", x.label ), y = "Density" ) +
    theme_bw() +
    theme( plot.margin = margin( asp$figure.margin,
                                 asp$figure.margin, asp$figure.margin,
                                 asp$figure.margin ),
           axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
           axis.text = element_text( size = asp$figure.axis.text.size ),
           axis.title = element_text( size = asp$figure.axis.title.size ),
           panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.position = "none" )

  ggsave( plot.file.path, plot = x.ggplot,
          width = asp$figure.width, height = asp$figure.height )
}
