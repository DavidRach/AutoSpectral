# plot_density_log.r

# Plots density of log values segregated by sign, ignoring diagonal
# coefficients and zeroes.

#' @title Plot Density of Log Values
#'
#' @description This function plots the density of log values, excluding matrix
#'     diagonal and zero values, using ggplot2.
#'
#' @importFrom ggplot2 ggplot aes geom_density scale_color_manual scale_fill_manual
#' @importFrom ggplot2 scale_linetype_manual labs theme_bw theme element_line
#' @importFrom ggplot2 element_text element_rect margin ggsave
#' @importFrom dplyr mutate
#' @importFrom rlang sym .data
#'
#' @param x.data Data frame containing values to be plotted.
#' @param x.label Label for the x-axis.
#' @param plot.file.path Path to save the plot file.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#'
#' @return Saves the plot as a JPEG file in the specified directory.
#' @export



refine.unmixing.log.plot <- function( x.data, x.label, plot.file.path, asp )
{
    # get values avoiding matrix diagonal
    x.data.n <- nrow( x.data )
    x.dat <- x.data[ -( 1 + 0:(x.data.n-1) * (x.data.n+1) ) ]

    # remove other zero values
    x.dat <- x.dat[ x.dat != 0 ]

    x.ggdata <- data.frame( log.abs.x = log10( abs( x.dat ) ),
        sign = factor( sign( x.dat ), levels = c( -1, 0, 1 ) ) )

    scale.color.value <- c( "-1" = asp$density.color.single,
        "1" = asp$density.color.single )

    scale.linetype.value <- c( "-1" = "dashed", "1" = "solid" )

    x.ggplot <- ggplot( x.ggdata, aes( x = .data$log.abs.x, color = .data$sign,
            fill = .data$sign, linetype = .data$sign ) ) +
        geom_density( alpha = 0.1, linewidth = asp$figure.density.line.size ) +
        scale_color_manual( values = scale.color.value ) +
        scale_fill_manual( values = scale.color.value ) +
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
