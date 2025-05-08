# plot_unmix_fix.r

#' @title Plot Unmixing Fixing
#'
#' @description Plots the data before and after automated unmixing corrections.
#'
#' @importFrom ggplot2 ggplot scale_x_continuous scale_y_continuous aes facet_wrap
#' @importFrom ggplot2 scale_color_gradientn guide_colorbar geom_path theme_classic
#' @importFrom ggplot2 labeller margin ggsave
#' @importFrom scattermore geom_scattermore
#' @importFrom dplyr %>% mutate
#' @importFrom KernSmooth bkde2D dpik
#' @importFrom rlang .data
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom scales trans_new
#'
#' @param original.data The fully stained sample with the original unmixing.
#' @param compensated.data The fully stained sample with compensation corrections
#'     applied.
#' @param plot.idx The index of the biplot comparisons to be plotted.
#' @param asp The AutoSpectral parameter list.
#' @param unstained.thresholds Numeric between 0 and 1, default 0.99. The threshold
#'     used to determine positivity in the fully stained sample (the percentile
#'     on the unstained sample in that channel).
#' @param unstained.margin Numeric, default 1.5. The fudge factor above the
#'     unstained.threshold.
#' @param spread.estimate Numeric matrix. The estimated spillover spread for the
#'     panel.
#' @param output.dir The folder where the plot will be saved.
#'
#' @return No returns. Creates a jpeg plot.
#'
#' @export
#'


plot.unmix.fix <- function( original.data, compensated.data, plot.idx,
                            asp, unstained.thresholds, unstained.margin,
                            spread.estimate, output.dir ) {

  # pull out marker names for plotting
  x.var <- colnames( compensated.data )[ plot.idx[ 1 ] ]
  y.var <- colnames( compensated.data )[ plot.idx[ 2 ] ]

  # add source column and combine data for plotting
  original.data <- original.data %>%
    data.frame( check.names = FALSE ) %>%
    mutate( source = "original" )
  fixed.data <- compensated.data %>%
    data.frame( check.names = FALSE ) %>%
    mutate( source = "fixed" )

  combined.data <- rbind( original.data, fixed.data )
  combined.data$source <- factor( combined.data$source,
                                  levels = c( "original", "fixed" ) )

  # get density for color scheme
  combined.data.density <- combined.data[ , c( x.var, y.var ) ]
  bandwidth.x <- suppressWarnings( asp$fix.unmix.bw.factor * dpik( combined.data.density[ , 1 ],
                                                                gridsize = asp$fix.unmix.grid.n )
  )
  bandwidth.y <- suppressWarnings( asp$fix.unmix.bw.factor * dpik( combined.data.density[ , 2 ],
                                                                gridsize = asp$fix.unmix.grid.n )
  )

  data.density <- suppressMessages( suppressWarnings( bkde2D(
    combined.data.density,
    bandwidth = c( bandwidth.x, bandwidth.y ),
    gridsize = c( asp$fix.unmix.grid.n,
                  asp$fix.unmix.grid.n ) ) ) )

  names( data.density ) <- c( "x", "y", "z" )

  # arrange data for plotting
  data.ggp <- data.frame(
    x = combined.data.density[ , 1 ],
    y = combined.data.density[ , 2 ],
    source = combined.data$source,
    z = interp.surface( data.density, combined.data.density ) )

  colnames( data.ggp )[ 1:2 ] <- c( x.var, y.var )

  # get color palette
  density.palette <- get.density.palette( data.ggp$z, asp )

  # move to asp
  spread.breaks <- asp$fix.unmix.plot.breaks

  # get spread coefficents and positivity thresholds
  spread.val.x <- spread.estimate[ x.var, y.var ]
  y.thr <- unstained.thresholds[ y.var ]

  spread.boundary.x <- data.frame(
    x = spread.breaks,
    y = y.thr*unstained.margin + abs( spread.breaks ) * spread.val.x
  )

  spread.val.y <- spread.estimate[ y.var, x.var ]
  x.thr <- unstained.thresholds[ x.var ]

  spread.boundary.y <- data.frame(
    x = x.thr*unstained.margin + abs( spread.breaks ) * spread.val.y,
    y = spread.breaks
  )

  biexp.transform <- flowjo_biexp( channelRange = 256,
                                   maxValue = asp$default.transformation.param$max.range,
                                   pos = asp$default.transformation.param$pos,
                                   neg = asp$default.transformation.param$neg,
                                   widthBasis = asp$default.transformation.param$width,
                                   inverse = FALSE )

  biexp.inverse <- flowjo_biexp( channelRange = 256,
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

  ribbon.labels <- sapply( asp$ribbon.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  ribbon.limits <- c( asp$default.transformation.param$width * 10, asp$expr.data.max )

  unmix.fix.plot <- ggplot( data.ggp,
                            aes( .data[[ x.var ]], .data[[ y.var ]],
                                 color = .data$z ) ) +
    scale_x_continuous( trans = plot.biexp.transform,
                        breaks = asp$ribbon.breaks,
                        limits = ribbon.limits,
                        labels = ribbon.labels ) +
    scale_y_continuous( trans = plot.biexp.transform,
                        breaks = asp$ribbon.breaks,
                        limits = ribbon.limits,
                        labels = ribbon.labels ) +
    geom_scattermore( pointsize = asp$figure.gate.point.size, alpha = 1,
                      na.rm = TRUE ) +
    scale_color_gradientn( "", labels = NULL, colors = density.palette,
                           guide = guide_colorbar( barwidth = asp$figure.gate.bar.width,
                                                   barheight = asp$figure.gate.bar.height ) ) +
    geom_path( aes( .data$x, .data$y, color = NULL ),
               data = spread.boundary.x, linewidth = asp$figure.gate.line.size ) +
    geom_path( aes( .data$x, .data$y, color = NULL ),
               data = spread.boundary.y, linewidth = asp$figure.gate.line.size ) +
    theme_classic() +
    facet_wrap( ~source, labeller = labeller( source = c( "original" = "Original",
                                                          "fixed" = "Fixed" ) ) ) +
    theme( plot.margin = margin( asp$figure.margin, asp$figure.margin,
                                 asp$figure.margin, asp$figure.margin ),
           legend.position = "none" )

  ggsave( filename = file.path( output.dir, asp$fix.unmixing.plot ),
          unmix.fix.plot, width = asp$fix.unmix.figure.width,
          height = asp$fix.unmix.figure.height )

}
