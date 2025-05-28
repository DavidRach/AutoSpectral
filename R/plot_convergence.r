# plot_convergence.r

#' @title Convergence Plot
#' @description Plots convergence of iterative refinement of the spillover matrix.
#'
#' @importFrom ggplot2 ggplot aes labs geom_hline geom_point scale_x_continuous
#' @importFrom ggplot2 scale_y_log10 scale_shape_manual theme_bw theme
#' @importFrom ggplot2 element_line element_text element_rect margin ggsave
#' @param convergence.log Dataframe with convergence data of AutoSpectral.
#' @param asp List with AutoSpectral parameters.
#' @return None. The function saves the generated convergence plot to a file.
#' @export

convergence.plot <- function( convergence.log, asp )
{
    convergence.ggdata <- convergence.log
    convergence.ggdata$delta.change <- abs( convergence.ggdata$delta.change )

    plot.xaxis.max  <- ceiling( max( convergence.log$iter ) / 10 ) * 10
    if ( plot.xaxis.max == 10 )
        plot.xaxis.step <- 5
    else if ( plot.xaxis.max <= 50 )
        plot.xaxis.step <- 10
    else
        plot.xaxis.step <- 25

    convergence.ggplot <- ggplot( convergence.ggdata, aes( x = .data$iter,
            shape = .data$scale ) ) +
        labs( x = "Iteration", y = "Absolute spillover error" ) +
        geom_hline( yintercept = asp$rs.delta.threshold.untr,
            linewidth = asp$figure.convergence.line.size,
            linetype = "dashed" ) +
        geom_hline( yintercept = asp$rs.delta.threshold.tran,
            linewidth = asp$figure.convergence.line.size,
            linetype = "dashed" ) +
        geom_hline( yintercept = asp$rs.delta.threshold.change,
            linewidth = asp$figure.convergence.line.size,
            linetype = "dashed" ) +
        geom_point( aes( y = .data$delta ),
            size = asp$figure.convergence.point.size,
            color = asp$convergence.color.delta ) +
        geom_point( aes( y = .data$delta.max ),
            size = asp$figure.convergence.point.size,
            color = asp$convergence.color.delta.max ) +
        geom_point( aes( y = .data$delta.change ),
            size = asp$figure.convergence.point.size,
            color = asp$convergence.color.delta.change ) +
        scale_x_continuous( limits = c( 0, plot.xaxis.max ),
            breaks = seq( 0, plot.xaxis.max, plot.xaxis.step ) ) +
        scale_y_log10() +
        scale_shape_manual( values = c(
            "linear" = asp$convergence.shape.linear,
            "bi-exp" = asp$convergence.shape.biexp,
            "posnegpop" = asp$convergence.shape.posnegpop
        ) ) +
        theme_bw() +
        theme( legend.position = "none",
            plot.margin = margin( asp$figure.margin, asp$figure.margin,
                asp$figure.margin, asp$figure.margin ),
            axis.ticks = element_line( linewidth = asp$figure.panel.line.size ),
            axis.text = element_text( size = asp$figure.axis.text.size ),
            axis.title = element_text( size = asp$figure.axis.title.size ),
            panel.border = element_rect( linewidth = asp$figure.panel.line.size ),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank() )

    ggsave( file.path( asp$figure.convergence.dir,
        sprintf( "%s.jpg", asp$convergence.file.name ) ),
        plot = convergence.ggplot,
        width = asp$figure.width, height = asp$figure.height )
}

