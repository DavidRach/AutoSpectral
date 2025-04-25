# refine_unmixing.r

#' @title Refine Unmixing
#'
#' @description This function refines the unmixing process for spectral data,
#'     iteratively updating the spectra matrix and calculating unmixing errors.
#'
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom stats sd
#' @importFrom utils write.csv
#'
#' @param spectra.initial Initial spectra matrix.
#' @param flow.control List containing flow control information.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param clean.expr Logical indicating whether to clean expression data.
#' @param plot.prefix Optional prefix for the plot title.
#'
#' @return The refined spectra matrix.
#' @export


refine.unmixing <- function( spectra.initial, flow.control, asp,
                                  clean.expr = FALSE, plot.prefix = "Refined" )
{

  if ( !is.null( plot.prefix ) ) {
    plot.title <- paste( plot.prefix, asp$spectra.file.name )
  } else {
    plot.title <- asp$spectra.file.name
  }

    fluorophores <- rownames( spectra.initial )

    if ( "AF" %in% fluorophores ) {
      af.control <- "AF"
    }

    peak.channels <- apply( spectra.initial, 1, which.max )

    biexp.transform <- flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue = asp$default.transformation.param$max.range,
      pos = asp$default.transformation.param$pos,
      neg = asp$default.transformation.param$neg,
      widthBasis = asp$default.transformation.param$width,
      inverse = FALSE
    )

    transform.inv <- flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue = asp$default.transformation.param$max.range,
      pos = asp$default.transformation.param$pos,
      neg = asp$default.transformation.param$neg,
      widthBasis = asp$default.transformation.param$width,
      inverse = TRUE
    )

    # set initial values for iteration variables
    rs.convergence <- FALSE
    rs.exit <- FALSE

    rs.iter <- 0
    rs.iter.last <- FALSE
    rs.iter.width <- floor( log10( asp$rs.iter.max ) ) + 1

    rs.lambda <- asp$rs.lambda.coarse

    rs.delta <- -1.0
    rs.delta.threshold <- asp$rs.delta.threshold.untr

    rs.delta.history <- rep( -1, asp$rs.delta.history.n )

    rs.scale.untransformed <- TRUE

    rs.convergence.log <- data.frame(
        iter = numeric(),
        scale = character(),
        lambda = numeric(),
        delta = numeric(),
        delta.max = numeric(),
        delta.change = numeric(),
        stringsAsFactors = FALSE
    )

    # set initial values for spectra calculation
    spectra.zero <- matrix( 0, nrow = length( fluorophores ),
                            ncol = flow.control$spectral.channel.n )
    rownames( spectra.zero ) <- fluorophores
    colnames( spectra.zero ) <- flow.control$spectral.channel

    for ( name in names( peak.channels ) ) {
      row <- which( rownames( spectra.zero ) == name )
      col <- peak.channels[ name ]
      spectra.zero[ row, col ] <- 1
    }

    spectra.curr <- spectra.zero
    spectra.update <- spectra.zero

    spectra.update <- spectra.initial - spectra.curr

    # get raw expression data and spectra
    if ( clean.expr ) {

      expr.data <- flow.control$clean.expr[ , flow.control$spectral.channel ]

      flow.event.sample <- flow.control$clean.event.sample[ flow.control$clean.event.sample
                                                            %in% fluorophores ]
      flow.event.sample <- factor( flow.event.sample, levels = unique( flow.event.sample ) )

    } else {

      expr.data <- flow.control$expr.data[ flow.control$event.sample %in% fluorophores,
                                           flow.control$spectral.channel ]

      # merge data
      flow.event.sample <- flow.control$event.sample[ flow.control$event.sample %in% fluorophores ]
      flow.event.sample <- factor( flow.event.sample, levels = unique( flow.event.sample ) )

    }

    while ( ! rs.exit )
    {
        # update spectra matrix and calculate unmixing matrix
        spectra.curr <- spectra.curr + spectra.update

        spectra.curr.original <- spectra.curr

        # get unmixed expression data
        expr.data.unmix <- unmix.ols( expr.data, spectra.curr.original )

        if ( ! rs.scale.untransformed ){
          expr.data.unmix <- apply( expr.data.unmix, 2, biexp.transform )
        }

        check.critical(
          identical( colnames( expr.data.unmix ),
                     rownames( spectra.curr.original ) ),
          "internal error: inconsistent dye names in unmixed data"
        )

        # get unmixing error in unmixed data
        unmixing.error <- get.unmixing.error(
          expr.data.unmix, fluorophores,
          rs.scale.untransformed, transform.inv,
          flow.event.sample, asp
        )

        # backconvert error to raw space
        M <- t( spectra.curr.original )
        unmixing.matrix.curr <- solve(t(M) %*% M) %*% t(M)

        unmixing.error.pseudo.inv <- solve( t( unmixing.error$slop )
                                            %*% unmixing.error$slop ) %*% t( unmixing.error$slop )


        unmixing.matrix.update <- t( t( unmixing.matrix.curr ) %*% unmixing.error.pseudo.inv )


        spectra.update.reverted <- solve( unmixing.matrix.update %*%
                                            t( unmixing.matrix.update) ) %*% unmixing.matrix.update

        spectra.update.reverted <- t( apply( spectra.update.reverted,
                                             1, function( x ) x/max( x ) ) )

        # get slope error and update delta variables
        slope.error <- spectra.update.reverted - spectra.curr.original

        rs.delta.prev <- rs.delta
        rs.delta <- sd( slope.error )

        rs.delta.max <- max( abs( slope.error ) )

        if ( rs.delta.prev >= 0 )
            rs.delta.history[ rs.iter %% asp$rs.delta.history.n + 1 ] <-
            rs.delta - rs.delta.prev
        else
            rs.delta.history[ rs.iter %% asp$rs.delta.history.n + 1 ] <- -1

        rs.delta.change <- mean( rs.delta.history )

        rs.convergence.log[ rs.iter + 1, ] <- list( rs.iter,
            ifelse( rs.scale.untransformed, "linear", "bi-exp" ),
            rs.lambda, rs.delta, rs.delta.max, rs.delta.change )

        if ( asp$verbose )
        {
            cat( sprintf(
                "iter %0*d, %s scale, lambda %.1f, delta %g, delta.max %g, delta.change %g\n",
                rs.iter.width, rs.iter,
                ifelse( rs.scale.untransformed, "linear", "bi-exp" ),
                rs.lambda, rs.delta, rs.delta.max, rs.delta.change ) )
        }

        # update iteration variables
        if ( rs.scale.untransformed && rs.delta.max < rs.delta.threshold )
        {
            # switch to bi-exponential scale and reset lambda and delta history
            rs.scale.untransformed <- FALSE
            rs.delta.threshold <- asp$rs.delta.threshold.tran
            rs.lambda <- asp$rs.lambda.coarse
            rs.delta <- -1.0
            rs.delta.history <- rep( -1, asp$rs.delta.history.n )
            rs.delta.change <- -1
        }

        if ( rs.delta.change > - asp$rs.delta.threshold.change &&
                rs.lambda == asp$rs.lambda.coarse )
        {
            # reduce lambda and reset delta history
            rs.lambda <- asp$rs.lambda.fine
            rs.delta <- -1.0
            rs.delta.history <- rep( -1, asp$rs.delta.history.n )
            rs.delta.change <- -1
        }

        rs.convergence <- ! rs.scale.untransformed &&
            ( rs.delta.max < rs.delta.threshold ||
                    rs.delta.change > - asp$rs.delta.threshold.change )

        rs.exit <- ( rs.convergence && rs.iter.last ) ||
            ( rs.delta.change > - asp$rs.delta.threshold.change &&
                    rs.scale.untransformed ) ||
            ( ! rs.convergence && rs.iter == asp$rs.iter.max ) ||
            rs.iter > asp$rs.iter.max

        rs.iter.last <- rs.convergence

        rs.iter <- rs.iter + 1

        # update spectra matrix
        spectra.update <- rs.lambda * slope.error

    }


    # save and plot convergence and error metrics
    # save and plot slope error
    if ( ! is.null( asp$table.slope.error.dir ) )
    {
      table.slope.error.file.name <- ifelse( rs.iter.last,
                                             sprintf( "%s.csv", asp$slope.error.file.name ),
                                             sprintf( "%s_%0*d.csv", asp$slope.error.file.name,
                                                      rs.iter.width, rs.iter ) )

      write.csv( slope.error,
                 file = file.path( asp$table.slope.error.dir,
                                   table.slope.error.file.name ) )

    }

    if ( ! is.null( asp$figure.slope.error.dir ) )
    {
      figure.slope.error.file.name <- sprintf( "%s%s.jpg",
                                               asp$slope.error.file.name,
                                               ifelse( rs.iter.last, "",
                                                       sprintf( "_%0*d", rs.iter.width, rs.iter ) ) )

      plot.density.log( slope.error, "unmixing error",
                        file.path( asp$figure.slope.error.dir,
                                   figure.slope.error.file.name ),
                        asp )
    }

    # save and plot skewness
    if ( ! is.null( asp$table.skewness.dir ) )
    {
      table.skewness.file.name <- ifelse( rs.iter.last,
                                          sprintf( "%s.csv", asp$skewness.file.name ),
                                          sprintf( "%s_%0*d.csv", asp$skewness.file.name,
                                                   rs.iter.width, rs.iter ) )

      write.csv( unmixing.error$skew,
                 file = file.path( asp$table.skewness.dir,
                                   table.skewness.file.name ) )
    }

    if ( ! is.null( asp$figure.skewness.dir ) )
    {
      if ( ! is.null( af.control ) )
        spillover.skewness <- unmixing.error$skew[
          - which( rownames( unmixing.error$skew ) == af.control ),
          - which( colnames( unmixing.error$skew ) == af.control ) ]
      else
        spillover.skewness <- unmixing.error$skew

      figure.skewness.file.name <- sprintf( "%s%s.jpg",
                                            asp$skewness.file.name,
                                            ifelse( rs.iter.last, "",
                                                    sprintf( "_%0*d", rs.iter.width, rs.iter ) ) )

      plot.density.log( spillover.skewness, "spillover skewness",
                        file.path( asp$figure.skewness.dir,
                                   figure.skewness.file.name ),
                        asp )
    }

    if ( ! is.null( asp$table.convergence.dir ) ){
      write.csv( rs.convergence.log,
                 file = file.path( asp$table.convergence.dir,
                 sprintf( "%s.csv", asp$convergence.file.name ) ),
                 row.names = FALSE )
    }


    if ( ! is.null( asp$figure.convergence.dir ) ){

      plot.convergence( rs.convergence.log, asp )
    }

    # save spectral matrix
    if ( ! is.null( asp$table.spectra.dir ) )
    {
      table.spectra.file.name <- sprintf( "%s.csv", plot.title )

      write.csv( spectra.update.reverted,
                 file = file.path( asp$table.spectra.dir,
                                   table.spectra.file.name ) )
    }

    # plot spectra
    if ( ! is.null( asp$figure.similarity.heatmap.dir ) ){
      plot.spectra( spectra.update.reverted, flow.control, asp, plot.title,
                    asp$figure.spectra.dir )
    }

    # plot similarity matrix heatmap
    if ( ! is.null( asp$figure.spectra.dir ) ){
      plot.similarity.matrix( spectra.update.reverted, asp, plot.prefix )
    }

    # check convergence
    check.critical( rs.convergence,
        "no convergence in refinement of spectra matrix" )

    return( spectra.update.reverted )
}

