# fix_my_unmix.r

#' @title Fix My Unmix
#' @description
#' Attempt to automatically fix unmixing errors in fully stained unmixed samples.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom flowCore read.FCS
#' @importFrom sp point.in.polygon
#' @importFrom utils write.csv
#' @importFrom flowWorkspace flowjo_biexp
#'
#' @param spectra The spectral matrix, fluorophores x detectors. In this case,
#' it does not need to be a perfect match for the fluorophores in the fully
#' stained sample and may, for instance, be pulled from a database for the
#' same specification of cytometer.
#' @param unstained.sample File path and name for a raw unstained sample,
#' ideally acquired the same day, and ideally matching the autofluorescence of
#' the fully stained sample.
#' @param fully.stained.sample File path and name for a raw fully stained sample.
#' @param asp The AutoSpectral parameter list.
#' @param verbose Logical. Turns messages on or off. Default is `TRUE`.
#' @param output.dir Character. Directory where plots and tables will be created.
#' Directory will be created if it does not exist. Default is `./FixMyUnmix`.
#' @param large.gate Logical, whether to use a large gate for identifying cells.
#' Default is `FALSE`. For panels containing any non-lymphocyte markers,
#' use `TRUE`.
#' @param max.iter Numeric. Limits the number of iterations performed by the
#' algorithm to avoid overfitting. Default is `100`.
#' @param gate.downsample.n Logical/numeric. Controls downsampling for the
#' definition of the gate to speed up the processing. If `FALSE`, no downsampling
#' will be performed for gating. If a numeric, that number of cells will be
#' used. Does not affect subsequent steps. Default is `1e5`.
#' @param model.downsample.n Logical/numeric. Controls downsampling for the
#' number of cells to be used in the algorithm itself. Default is `2e5`. For
#' samples with rare markers, use a larger number, which will be somehwat slower.
#' @param pos.threshold Numeric between 0 and 1, default `0.99` The threshold
#' used to determine positivity in the fully stained sample (the percentile
#' on the unstained sample in that channel).
#' @param neg.threshold Numeric between 0 and 1, default `0.8`. The threshold
#' used to determine negativity in the fully stained sample (the percentile
#' on the unstained sample in that channel). Deprecated.
#' @param unstained.margin Numeric, default `1.3`. The fudge factor above the
#' unstained.threshold.
#' @param convergence.threshold Numeric, default `0.005` corresponding to half
#' a percent. If the algorithm reaches this minimal level of spillover error,
#' it will stop.
#' @param biexp Logical indicating whether to apply biexponential
#' transformation. Default is `FALSE`.
#' @param rs.lambda Numeric, default 0.1. Controls how rapidly the spillover
#' corrections are applied in order to reach convergence.
#' @param rs.delta.history.n Numeric. Number of iterations to track in order to
#' determine convergence. Default is `10`.
#' @param rs.delta.threshold.change Numeric. Threshold to determine convergence
#' based on a lack of change between iterations. Default is `1e-6`.
#' @param min.neg.n Numeric. The minimum number of cells (points) in the
#' negative region required to calculate corrections. Default is `50`. Channels
#' with fewer events will not be considered.
#' @param max.neg.n Numeric. The maximum number of negative events to use for
#' the calculation. Default is `500`. Larger numbers may improve accuracy in
#' messy samples, at the cost of increased calculation time.
#' @param min.pos.n Numeric. The minimum number of cells (points) in the
#' positive region required to calculate corrections. Default is `50`. Channels
#' with fewer events will not be considered.
#' @param max.pos.n Numeric. The maximum number of positive events to use for
#' the calculation. Default is `100`. The brightest n events in each channel
#' will be used, so this controls selection of the positive events and
#' processing time.
#' @param rlm.max.iter Numeric. The maximum number of iterations for the
#' internal robust linear model call. Default is `100`.
#' @param scatter.channels The names of the FSC and SSC channels to be used for
#' gating. Must match parameter names from the FCS files. Default is
#' `flow.control$scatter.parameter`, which requires the use of `AutoSpectral`,
#' either `define.flow.control` or `reload.flow.control`.
#' @param spectral.channels The names of the fluorescent detector channels.
#' Must match parameter names from the FCS files. Default is
#' `flow.control$spectral.channel`, which requires the use of `AutoSpectral`,
#' either `define.flow.control` or `reload.flow.control`.
#' @param scatter.and.spectral The names of the scatter and fluorescent detector
#' channels. Must match parameter names from the FCS files. Default is
#' `flow.control$scatter.and.channel.spectral`, which requires the use of
#' `AutoSpectral`, either `define.flow.control` or `reload.flow.control`.
#' @param channel.labels Labels for the scatter channels and peak channels for
#' the fluorophores. Used for gating plots. Default is
#' `flow.control$scatter.and.channel.label`, which requires the use of
#' `AutoSpectral`, either `define.flow.control` or `reload.flow.control`.
#' @param heatmap.title Character. Title for the saved heatmap showing identified
#' unmixing errors. Default is `Fixing_unmixing`.
#' @param spillover.filename Character. Name of the CSV file carrying the
#' spillover (fluorophore x fluorophore in the unmixed space) values. Default
#' is `Fixed_spillover.csv`.
#' @param spectra.filename Character. Name of the CSV file carrying the fixed
#' spectra (fluorophore x detectors in the raw space) values. Default
#' is `Fixed_spectra.csv`.
#' @param compensation.filename Character. Name of the CSV file carrying the
#' compensation matrix (fluorophore x fluorophore in the unmixed space). Default
#' is `Fixed_compensation.csv`.
#' @param flowjo.mtx.filename Character. Name of the XML .mtx file carrying the
#' compensation matrix for use in FlowJo. Default is `Fixed_compensation.mtx`.
#' @param heatmap.number.labels Logical. Whether to include numbers on each
#' square tile of the heatmap of unmixing errors showing the spillover value.
#' Default is `FALSE` for better aesthetics.
#' @param triangular.heatmap Logical. Whether to produce a triangular `TRUE` or
#' square `FALSE` heatmap. Default is `TRUE`.
#' @param heatmap.label Character. Title to put in the legend of the heatmap of
#' identified unmixing errors. Default is `Spillover`.
#' @param heatmap.color.palette Optional character string defining the viridis
#' color palette to be used for the heatmap. Default is `viridis`. Options are
#' the viridis color options: `magma`,  `inferno`, `plasma`, `viridis`,
#' `cividis`, `rocket`, `mako` and `turbo`.
#' @param heatmap.fig.width Numeric. Width of the heatmap figure.
#' Default is `8`.
#' @param heatmap.fig.height Numeric. Height of the heatmap figure.
#' Default is `6`.
#' @param spectra.color.palette Optional character string defining the viridis
#' palette to be used for the fluorophore traces. Default is `NULL`, in which
#' case default R Brewer colors will be assigned automatically. Options are the
#' viridis color options: `magma`, `inferno`, `plasma`, `viridis`, `cividis`,
#' `rocket`, `mako` and `turbo`.
#'
#' @return A list containing the corrected spillover matrix (fluorophores x
#' fluorophores), the corrected compensation matrix and the corrected spectra
#' for unmixing (fluorophores x detectors).
#'
#' @export


fix.my.unmix <- function( spectra, unstained.sample, fully.stained.sample,
                          asp,
                          verbose = TRUE,
                          output.dir = "./FixMyUnmix",
                          large.gate = FALSE,
                          max.iter = 100,
                          gate.downsample.n = 1e5, model.downsample.n = 2e5,
                          pos.threshold = 0.99, neg.threshold = 0.8,
                          unstained.margin = 1.3,
                          convergence.threshold = 0.005, biexp = FALSE,
                          rs.lambda = 0.1, rs.delta.history.n = 10,
                          rs.delta.threshold.change = 1e-6,
                          min.neg.n = 50, max.neg.n = 500,
                          min.pos.n = 50, max.pos.n = 100,
                          rlm.max.iter = 100,
                          scatter.channels = flow.control$scatter.parameter,
                          spectral.channels = flow.control$spectral.channel,
                          scatter.and.spectral = flow.control$scatter.and.channel.spectral,
                          channel.labels = flow.control$scatter.and.channel.label,
                          heatmap.title = "Fixing_unmixing",
                          spillover.filename = "Fixed_spillover.csv",
                          spectra.filename = "Fixed_spectra.csv",
                          compensation.filename = "Fixed_compensation.csv",
                          flowjo.mtx.filename = "Fixed_compensation.mtx",
                          heatmap.number.labels = FALSE,
                          triangular.heatmap = TRUE,
                          heatmap.label = "Spillover",
                          heatmap.color.palette = "plasma",
                          heatmap.fig.width = 8,
                          heatmap.fig.height = 6,
                          spectra.color.palette = "magma"
                          ) {

  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

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
  similarity.matrix[ similarity.matrix <= 0 ] <- 1e-6

  # calculate unmixing matrix
  unmixing.matrix <- solve( spectra %*% t( spectra ) ) %*% spectra

  # estimate spread matrix
  spread.estimate <- estimate.spread( spectra )

  # import unstained sample
  unstained.raw <- suppressWarnings(
    flowCore::read.FCS( unstained.sample, transformation = FALSE,
                        truncate_max_range = FALSE, emptyValue = FALSE )
  )
  unstained.raw <- flowCore::exprs( unstained.raw )[ , scatter.and.spectral ]

  if ( verbose )
    message( "\033[32mGating raw data.\033[0m" )

  # set up downsampling better with set.seed, sample, idx
  if ( nrow( unstained.raw ) > gate.downsample.n ) {
    set.seed( 42 )
    gate.idx <- sample( 1:nrow( unstained.raw ), gate.downsample.n )
    gate.data <- unstained.raw[ gate.idx, scatter.channels ]
  } else {
    gate.data <- unstained.raw[ , scatter.channels ]
  }

  # define large gate
  gate <- do.gate( gate.data, viability.gate = FALSE, large.gate = large.gate,
                   samp = "unstained raw",
                   scatter.and.channel.label = channel.labels,
                   control.type = "cells", asp )

  # get cells in gate
  gate.population.pip <- point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ],
                                           gate$x, gate$y )
  gate.population.idx <- which( gate.population.pip != 0 )
  unstained.raw <- unstained.raw[ gate.population.idx, spectral.channels ]

  # unmix
  if ( verbose )
    message( "\033[32mUnmixing unstained sample.\033[0m" )

  unstained.unmixed <- unmix.ols( unstained.raw, spectra )

  if ( biexp )
    unstained.unmixed <- apply( unstained.unmixed, 2, biexp.transform )

  # calculate nth percentiles across all channels
  pos.thresholds <- apply( unstained.unmixed, 2, function( col )
    quantile( col, probs = pos.threshold ) )
  neg.thresholds <- apply( unstained.unmixed, 2, function( col )
    quantile( col, probs = neg.threshold ) )


  # import fully stained sample
  if ( verbose )
    message( "\033[32mReading fully stained raw data.\033[0m" )

  fully.stained.raw <- suppressWarnings(
    flowCore::read.FCS( fully.stained.sample, transformation = FALSE,
                        truncate_max_range = FALSE, emptyValue = FALSE )
  )
  fully.stained.raw <- flowCore::exprs( fully.stained.raw )[ , scatter.and.spectral ]

  # get cells in gate
  gate.data <- fully.stained.raw[ , scatter.channels ]
  gate.population.pip <- point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ],
                                           gate$x, gate$y )
  gate.population.idx <- which( gate.population.pip != 0 )
  fully.stained.raw <- fully.stained.raw[ gate.population.idx, spectral.channels ]

  # unmix
  if ( verbose )
    message( "\033[32mUnmixing fully stained sample.\033[0m" )

  if ( nrow( fully.stained.raw ) > model.downsample.n ) {
    set.seed( 42 )
    downsample.idx <- sample( 1:nrow( fully.stained.raw ), model.downsample.n )
    fully.stained.raw <- fully.stained.raw[ downsample.idx, ]
  }

  fully.stained.unmixed <- unmix.ols( fully.stained.raw, spectra )
  rm( fully.stained.raw )
  rm( unstained.raw )

  if ( biexp )
    fully.stained.unmixed <- apply( fully.stained.unmixed, 2, biexp.transform )

  # set initial values for iteration variables
  rs.convergence <- FALSE
  rs.exit <- FALSE
  rs.iter <- 0
  rs.iter.width <- floor( log10( max.iter ) ) + 1

  rs.delta <- -1.0
  rs.delta.history <- rep( -1, rs.delta.history.n )

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
    unmixed.comp <- fully.stained.unmixed %*% compensation.curr

    # get spillover coefficients
    marker.spillover <- lapply( fluorophores, function( fl ) {

      # select peak channel data
      peak.channel.expr <- unmixed.comp[ , fl ]

      fluor.spectra.coef <- spectra.zero

      for ( channel in fluorophores ) {

        if ( channel == fl ) {
          fluor.spectra.coef[ channel ] <- 1.0

        } else {
          # get data for spillover channel
          channel.expr <- unmixed.comp[ , channel ]

          # select negative events in spillover channel using estimated spread
          #channel.thr <- neg.thresholds[ fl ] * unstained.margin
          channel.thr <- pos.thresholds[ channel ] * unstained.margin
          channel.spread <- spread.estimate[ channel, fl ]

          channel.expr.idx <- which( unmixed.comp[ , channel ] <
                                      ( unmixed.comp[ , fl ] * channel.spread + channel.thr ) )

          # if everything is positive in this channel, skip
          if ( length( channel.expr.idx ) < min.neg.n ) {
            if ( verbose )
              message( paste( "Unable to identify enough negative events for",
                              fl, "versus", channel ) )

            fluor.spectra.coef[ channel ] <- 0
            next

          } else {
            # select peak channel data within the negative region for the spillover channel
            peak.channel.data <- peak.channel.expr[ channel.expr.idx ]

            # check that some of these are above the positivity threshold for the peak channel
            #peak.channel.thr <- neg.thresholds[ channel ] * unstained.margin
            peak.channel.thr <- pos.thresholds[ fl ] * unstained.margin
            rel.peak.ch.pos.idx <- which( peak.channel.data > peak.channel.thr )
            abs.peak.ch.pos.idx <- channel.expr.idx[ rel.peak.ch.pos.idx ]
            pos.n <- length( abs.peak.ch.pos.idx )

            # check that there are enough positives to work with
            if ( pos.n < 50 ) {
              if ( verbose )
                message( paste( "Unable to identify enough positive events for",
                                fl, "versus", channel ) )

              fluor.spectra.coef[ channel ] <- 0
              next

            } else {
              # select peak of positive events to use for model
              peak.pos.values <- peak.channel.expr[ abs.peak.ch.pos.idx ]
              peak.pos.ord <- order( peak.pos.values, decreasing = TRUE )[ 1:min( max.pos.n, pos.n ) ]
              peak.pos.idx <- abs.peak.ch.pos.idx[ peak.pos.ord ]

              #abs.peak.ch.neg.idx <- setdiff( channel.expr.idx, abs.peak.ch.pos.idx )
              rel.neg.idx <- which( unmixed.comp[ channel.expr.idx, fl ] < pos.thresholds[ fl ] )
              abs.peak.ch.neg.idx <- channel.expr.idx[ rel.neg.idx ]

              neg.n <- length( abs.peak.ch.neg.idx )

              if ( neg.n < min.neg.n ) {
                if ( verbose )
                  message( paste( "Not enough negatives for", fl, "vs", channel ) )

                fluor.spectra.coef[ channel ] <- 0
                next

              }
              # select max number of negative events to use for model
              ### TBD: random sampling ###
              peak.neg.idx <- abs.peak.ch.neg.idx[ 1:min( max.neg.n, neg.n ) ]

              # select positive data
              peak.channel.data <- peak.channel.expr[ c( peak.pos.idx, peak.neg.idx ) ]
              # select negatives from spillover channel
              channel.expr.data <- channel.expr[ c( peak.pos.idx, peak.neg.idx ) ]
            }
          }

          # fit robust linear model
          spectra.model.result <- fit.robust.linear.model(
            peak.channel.data, channel.expr.data, fl, channel,
            rlm.max.iter, fix.unmix = TRUE )

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
    similarity.matrix[ similarity.matrix <= 0 ] <- 1e-6
    ratioed.coefficients <- marker.spillover * sqrt( similarity.matrix )

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

    if ( rs.delta.prev >= 0 ) {
      rs.delta.history[ rs.iter %% rs.delta.history.n + 1 ] <-
        rs.delta - rs.delta.prev
    } else {
      rs.delta.history[ rs.iter %% rs.delta.history.n + 1 ] <- -1
    }

    rs.delta.change <- mean( rs.delta.history )

    rs.convergence.log[ rs.iter + 1, ] <- list( rs.iter, rs.delta, rs.delta.max,
                                                rs.delta.change )

    if ( verbose ) {
      message( sprintf( "iter %0*d, delta %g, delta.max %g, delta.change %g",
        rs.iter.width, rs.iter, rs.delta, rs.delta.max, rs.delta.change ) )
    }

    # detect convergence
    rs.convergence.now <- ( rs.delta.max < convergence.threshold )
    ### check sign here ###
    rs.stalled <- ( rs.delta.change > -rs.delta.threshold.change )

    # determine if exit conditions have been met
    if ( rs.convergence.now && rs.convergence ) {
      rs.exit <- TRUE
      if ( verbose )
        message( "\033[34mConverged \033[0m" )
    } else if ( rs.iter >= max.iter ) {
      rs.exit <- TRUE
      if ( verbose )
        message( "\033[33m Reached iteration limit \033[0m" )
    } else if ( rs.stalled ) {
      rs.exit <- TRUE
      if ( verbose )
        message( "\033[31m Refinement stalled, stopping \033[0m" )
    }

    # save convergence state (allows one extra pass)
    rs.convergence <- rs.convergence.now

    # update spillover matrix
    spillover.update <- slope.error %*% spillover.curr * rs.lambda

    rs.iter <- rs.iter + 1
  }

  # heatmap of coefficients
  create.heatmap( spillover.curr, number.labels = heatmap.number.labels,
                  plot.prefix = heatmap.title,
                  legend.label = heatmap.label, triangular = triangular.heatmap,
                  output.dir = output.dir, color.palette = heatmap.color.palette,
                  figure.width = heatmap.fig.width,
                  figure.height = heatmap.fig.height )

  # write spillover as csv
  write.csv( spillover.curr, file = file.path( output.dir, spillover.filename ) )

  # write compensation as csv
  compensation.matrix <- solve( spillover.curr )
  write.csv( compensation.matrix, file = file.path( output.dir,
                                                    compensation.filename ) )

  # write compensation as flowjo.mtx
  write.flowjo.mtx( spillover.curr,
                    file.path( output.dir, flowjo.mtx.filename ) )

  if ( biexp ) {
    inverse.transform <- flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue = asp$default.transformation.param$max.range,
      pos = asp$default.transformation.param$pos,
      neg = asp$default.transformation.param$neg,
      widthBasis = asp$default.transformation.param$width,
      inverse = TRUE
    )

    fully.stained.unmixed <- apply( fully.stained.unmixed, 2, inverse.transform )
    unmixed.comp <- fully.stained.unmixed %*% compensation.matrix
    unstained.unmixed <- apply( unstained.unmixed, 2, inverse.transform )
    neg.thresholds <- apply( unstained.unmixed, 2, function( col )
      quantile( col , probs = neg.threshold ) )
  }

  # plot example of fixing on marker pair with biggest error
  unmix.fix.plot( fully.stained.unmixed, unmixed.comp, plot.idx, asp,
                  neg.thresholds, unstained.margin,
                  spread.estimate, output.dir = output.dir )

  # back convert to unmixing spectra
  M <- t( spectra )
  unmixing.matrix <- solve( t( M ) %*% M ) %*% t( M )

  unmixing.error.pseudo.inv <- solve( t( spillover.curr )
                                      %*% spillover.curr ) %*% t( spillover.curr )

  unmixing.matrix.update <- t( t( unmixing.matrix ) %*% unmixing.error.pseudo.inv )

  spectra.update.reverted <- solve( unmixing.matrix.update %*%
                                      t( unmixing.matrix.update) ) %*% unmixing.matrix.update

  # spectra.update.reverted <- spillover.curr %*% spectra

  spectra.update.reverted <- t( apply( spectra.update.reverted,
                                       1, function( x ) x/max( x ) ) )

  write.csv( spectra.update.reverted, file = file.path( output.dir, spectra.filename ) )

  # plot original spectra and adjusted spectra
  spectral.trace( spectral.matrix = spectra,
                  plot.title = "Original spectra",
                  plot.dir = output.dir, split.lasers = TRUE,
                  asp$figure.spectra.line.size,
                  asp$figure.spectra.point.size,
                  color.palette = spectra.color.palette )
  spectral.trace( spectral.matrix = spectra.update.reverted,
                  plot.title = "Adjusted spectra",
                  plot.dir = output.dir, split.lasers = TRUE,
                  asp$figure.spectra.line.size,
                  asp$figure.spectra.point.size,
                  color.palette = spectra.color.palette )

  # return updated spectra, spillover & compensation
  return( list(
    Spillover = spillover.curr,
    Compensation = compensation.matrix,
    Spectra = spectra.update.reverted
  ) )

}
