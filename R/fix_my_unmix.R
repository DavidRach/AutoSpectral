# fix_my_unmix.r

#' @title Fix My Unmix
#' @description Attempt to automatically fix unmixing errors in fully stained
#'     unmixed samples.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom flowCore read.FCS
#' @importFrom sp point.in.polygon
#' @importFrom utils write.csv
#' @importFrom flowWorkspace flowjo_biexp
#'
#' @param spectra The spectral matrix, fluorophores x detectors. In this case,
#'     it does not need to be a perfect match for the fluorophores in the fully
#'     stained sample and may, for instance, be pulled from a database for the
#'     same specification of cytometer.
#' @param unstained.sample File path and name for a raw unstained sample, ideally
#'     acquired the same day, and ideally matching the autofluorescence of the
#'     fully stained sample.
#' @param fully.stained.sample File path and name for a raw fully stained sample.
#' @param flow.control The flow.control list.
#' @param asp The AutoSpectral parameter list.
#' @param large.gate Logical, whether to use a large gate for identifying cells.
#'     Default is TRUE and will work best for any panels containing any non-
#'     lymphocyte markers.
#' @param max.iter Numeric. Limits the number of iterations performed by the algorithm
#'     to avoid overfitting.
#' @param downsample Logical/numeric. Controls downsampling to speed up the algorithm.
#'     If FALSE, no downsampling will be performed. If a numeric, that number of
#'     cells will be used for downsampling. For rare markers, use larger numbers.
#' @param unstained.threshold Numeric between 0 and 1, default 0.99. The threshold
#'     used to determine positivity in the fully stained sample (the percentile
#'     on the unstained sample in that channel).
#' @param unstained.margin Numeric, default 1.3. The fudge factor above the
#'     unstained.threshold.
#' @param convergence.threshold Numeric, default 0.01. If the algorithm reaches
#'     this minimal level of spillover error, it will stop.
#' @param biexp Logical indicating whether to apply biexponential transformation
#'     (default is TRUE).
#'
#' @return A list containing the corrected spillover matrix (fluorophores x fluorophores),
#'     the corrected compensation matrix and the corrected spectra for unmixing
#'     (fluorophores x detectors).
#'
#' @export


fix.my.unmix <- function( spectra, unstained.sample, fully.stained.sample,
                          flow.control, asp, large.gate = TRUE,
                          max.iter = 20, downsample = 20000,
                          unstained.threshold = 0.99, unstained.margin = 1.3,
                          convergence.threshold = 0.03, biexp = TRUE ) {

  if ( ! dir.exists( asp$fix.unmixing.dir ) )
    dir.create( asp$fix.unmixing.dir )

  if ( biexp ) {
    biexp.transform <- flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue = asp$default.transformation.param$max.range,
      pos = asp$default.transformation.param$pos,
      neg = asp$default.transformation.param$neg,
      widthBasis = asp$default.transformation.param$width,
      inverse = FALSE
    )
  }

  fluorophores <- rownames( spectra )
  fluorophore.n <- length( fluorophores )

  spectra <- as.matrix( spectra )

  # calculate similarity matrix
  similarity.matrix <- cosine.similarity( t( spectra ) )

  # calculate unmixing matrix
  unmixing.matrix <- solve( spectra %*% t( spectra ) ) %*% spectra

  # estimate spread matrix
  spread.estimate <- estimate.spread( spectra )

  # import unstained sample
  unstained.raw <- suppressWarnings(
    flowCore::read.FCS( unstained.sample, transformation = FALSE,
                        truncate_max_range = FALSE, emptyValue = FALSE )
  )
  unstained.raw <- flowCore::exprs( unstained.raw )[ , flow.control$scatter.and.channel.spectral ]

  message( "\033[32mGating raw data.\033[0m" )
  gate.data <- unstained.raw[ , flow.control$scatter.parameter ]

  # define large gate
  gate <- do.gate( gate.data, viability.gate = FALSE, large.gate = large.gate,
                   samp = "unstained raw",
                   scatter.and.channel.label = flow.control$scatter.and.channel.label,
                   control.type = "cells", asp )

  # get cells in gate
  gate.population.pip <- point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ],
                                           gate$x, gate$y )
  gate.population.idx <- which( gate.population.pip != 0 )
  unstained.raw <- unstained.raw[ gate.population.idx, flow.control$spectral.channel ]

  # unmix
  message( "\033[32mUnmixing unstained sample.\033[0m" )
  unstained.unmixed <- unmix.ols( unstained.raw, spectra )

  if ( biexp )
    unstained.unmixed <- apply( unstained.unmixed, 2, biexp.transform )

  # calculate nth percentiles across all channels
  unstained.thresholds <- apply( unstained.unmixed, 2, function( col )
    quantile( col , probs = unstained.threshold ) )

  # import fully stained sample
  message( "\033[32mReading fully stained raw data.\033[0m" )
  fully.stained.raw <- suppressWarnings(
    flowCore::read.FCS( fully.stained.sample, transformation = FALSE,
                        truncate_max_range = FALSE, emptyValue = FALSE )
  )
  fully.stained.raw <- flowCore::exprs( fully.stained.raw )[ , flow.control$scatter.and.channel.spectral ]

  # get cells in gate
  gate.data <- fully.stained.raw[ , flow.control$scatter.parameter ]
  gate.population.pip <- point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ],
                                           gate$x, gate$y )
  gate.population.idx <- which( gate.population.pip != 0 )
  fully.stained.raw <- fully.stained.raw[ gate.population.idx, flow.control$spectral.channel ]

  # unmix
  message( "\033[32mUnmixing fully stained sample.\033[0m" )
  fully.stained.unmixed <- unmix.ols( fully.stained.raw, spectra )
  # rm( fully.stained.raw )
  # rm( unstained.raw )

  downsample.idx <- sample( 1:nrow( fully.stained.unmixed ), downsample )
  fully.stained.downsampled <- fully.stained.unmixed[ downsample.idx, ]

  if ( biexp )
    fully.stained.downsampled <- apply( fully.stained.downsampled, 2, biexp.transform )

  # set initial values for iteration variables
  rs.convergence <- FALSE
  rs.exit <- FALSE
  rs.iter <- 0
  rs.iter.width <- floor( log10( max.iter ) ) + 1

  rs.delta <- -1.0
  rs.delta.history <- rep( -1, asp$rs.delta.history.n )

  rs.convergence.log <- data.frame(
    iter = numeric(),
    delta = numeric(),
    delta.max = numeric(),
    delta.change = numeric(),
    stringsAsFactors = FALSE
  )

  # set up compensation matrix zero
  spectra.zero <- rep( 0, length( fluorophores ) )
  names( spectra.zero ) <- fluorophores

  # set initial values for spillover calculation
  spillover.curr <- diag( fluorophore.n )
  rownames( spillover.curr ) <- fluorophores
  colnames( spillover.curr ) <- fluorophores
  spillover.update <- spillover.curr - diag( fluorophore.n )

  # main loop
  while ( ! rs.exit ) {

    # update spillover matrix and calculate compensation matrix
    spillover.curr <- spillover.curr + spillover.update
    spillover.curr <- sweep( spillover.curr, 1, diag( spillover.curr ), "/" )

    compensation.curr <- solve( spillover.curr )

    # apply compensation
    unmixed.comp <- fully.stained.downsampled %*% compensation.curr

    # get spillover coefficients
    marker.spillover <- lapply( fluorophores, function( fl ) {

      peak.channel.expr <- unmixed.comp[ , fl ]

      fluor.spectra.coef <- spectra.zero

      for ( channel in fluorophores ) {
        if ( channel == fl ) {
          fluor.spectra.coef[ channel ] <- 1.0
        } else {
          channel.expr <- unmixed.comp[ , channel ]

          # select negative events using estimated spread
          channel.thr <- unstained.thresholds[ fl ] * unstained.margin
          channel.spread <- spread.estimate[ channel, fl ]

          channel.expr.idx <- which( channel.expr <
                                       ( channel.thr + peak.channel.expr * channel.spread ) )

          channel.expr <- channel.expr[ channel.expr.idx ]
          peak.channel.data <- peak.channel.expr[ channel.expr.idx ]

          # fit robust linear model
          spectra.model.result <- fit.robust.linear.model(
            peak.channel.data, channel.expr, fl, channel, asp, fix.unmix = TRUE )

          fluor.spectra.coef[ channel ] <- spectra.model.result[ 2 ]
        }
      }

      # normalize fluor.spectra.coef
      fluor.spectra.coef <- fluor.spectra.coef / max( fluor.spectra.coef )

      fluor.spectra.coef
    })

    marker.spillover <- do.call( rbind, marker.spillover )
    rownames( marker.spillover ) <- fluorophores

    # ratio by similarity matrix
    ratioed.coefficients <- marker.spillover * similarity.matrix

    # select worst pair for plotting on first run
    if ( rs.iter == 0 ) {
      plot.select <- ratioed.coefficients
      diag( plot.select ) <- 0
      plot.select <- abs( plot.select )
      plot.idx <- which( plot.select == max( plot.select ), arr.ind = TRUE )
    }

    # update spillover matrix & calculate compensation matrix
    compensation.matrix <- solve( ratioed.coefficients )

    # get slope error and update delta variables
    slope.error <- ratioed.coefficients - diag( fluorophore.n )

    rs.delta.prev <- rs.delta
    rs.delta <- sd( slope.error )
    rs.delta.max <- max( abs( slope.error ) )

    if ( rs.delta.prev >= 0 )
      rs.delta.history[ rs.iter %% asp$rs.delta.history.n + 1 ] <-
      rs.delta - rs.delta.prev
    else
      rs.delta.history[ rs.iter %% asp$rs.delta.history.n + 1 ] <- -1

    rs.delta.change <- mean( rs.delta.history )

    rs.convergence.log[ rs.iter + 1, ] <- list( rs.iter, rs.delta, rs.delta.max,
                                                rs.delta.change )

    if ( asp$verbose )
    {
      message( sprintf( "iter %0*d, delta %g, delta.max %g, delta.change %g",
        rs.iter.width, rs.iter, rs.delta, rs.delta.max, rs.delta.change ) )
    }

    # detect convergence
    rs.convergence.now <- ( rs.delta.max < convergence.threshold )
    rs.stalled <- ( rs.delta.change > -asp$rs.delta.threshold.change )

    # determine if exit conditions have been met
    if ( rs.convergence.now && rs.convergence ) {
      rs.exit <- TRUE
      message( "\033[34mConverged \033[0m" )
    } else if ( rs.iter >= max.iter ) {
      rs.exit <- TRUE
      message( "\033[33m Reached iteration limit \033[0m" )
    } else if ( rs.stalled ) {
      rs.exit <- TRUE
      message( "\033[31m Refinement stalled, stopping \033[0m" )
    }

    # save convergence state (allows one extra pass)
    rs.convergence <- rs.convergence.now

    # update spillover matrix
    spillover.update <- slope.error %*% spillover.curr * 0.1

    rs.iter <- rs.iter + 1
  }

  # heatmap of coefficients
  create.heatmap( spillover.curr, asp, number.labels = FALSE,
                plot.prefix = asp$fix.unmixing.heatmap,
                legend.label = "Spillover", triangular = TRUE,
                output.dir = asp$fix.unmixing.dir )

  # write spillover as csv
  write.csv( spillover.curr, file = file.path( asp$fix.unmixing.dir,
                                               asp$fix.spillover.filename ) )

  # write compensation as csv
  compensation.matrix <- solve( spillover.curr )
  write.csv( compensation.matrix, file = file.path( asp$fix.unmixing.dir,
                                                    asp$fix.compensation.filename ) )

  # write compensation as flowjo.mtx
  write.flowjo.mtx( spillover.curr,
                    file.path( asp$fix.unmixing.dir, asp$fix.compensation.flowjo ) )

  if ( biexp ) {
    inverse.transform <- flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue = asp$default.transformation.param$max.range,
      pos = asp$default.transformation.param$pos,
      neg = asp$default.transformation.param$neg,
      widthBasis = asp$default.transformation.param$width,
      inverse = TRUE
    )

    fully.stained.downsampled <- apply( fully.stained.downsampled, 2, inverse.transform )
    unmixed.comp <- fully.stained.downsampled %*% compensation.matrix
    unstained.unmixed <- apply( unstained.unmixed, 2, inverse.transform )
    unstained.thresholds <- apply( unstained.unmixed, 2, function( col )
      quantile( col , probs = unstained.threshold ) )
  }

  # plot example of fixing on marker pair with biggest error
  unmix.fix.plot( fully.stained.downsampled, unmixed.comp, plot.idx, asp,
                  unstained.thresholds, unstained.margin,
                  spread.estimate, output.dir = asp$fix.unmixing.dir )

  # back convert to unmixing spectra
  spectra.update.reverted <- spillover.curr %*% spectra

  spectra.update.reverted <- t( apply( spectra.update.reverted,
                                       1, function( x ) x/max( x ) ) )
  write.csv( spectra.update.reverted, file = file.path( asp$fix.unmixing.dir,
                                                    asp$fix.spectra.filename ) )

  # plot original spectra and adjusted spectra
  spectral.trace( spectra, flow.control, asp, plot.title = "Original spectra",
                plot.dir = asp$fix.unmixing.dir, split.lasers = TRUE )
  spectral.trace( spectra.update.reverted, flow.control, asp,
                plot.title = "Adjusted spectra",
                plot.dir = asp$fix.unmixing.dir, split.lasers = TRUE )

  # return updated spectra, spillover & compensation
  return( list(
    Spillover = spillover.curr,
    Compensation = compensation.matrix,
    Spectra = spectra.update.reverted
  ) )

}
