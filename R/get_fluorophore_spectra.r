# get_fluorophore_spectra.r


#' Get Fluorophore Spectra
#'
#' This function retrieves the fluorophore spectra for flow cytometry data,
#'     optionally using cleaned expression data and biexponential transformation.
#'     It also plots and saves the spectra, and performs cosine similarity QC for
#'     controls.
#'
#' @title Get Fluorophore Spectra
#' @description Retrieves the fluorophore spectra for flow cytometry data,
#'     optionally using cleaned expression data and biexponential transformation.
#'
#' @importFrom flowWorkspace flowjo_biexp
#' @param flow.control A list containing flow cytometry control data.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param use.clean.expr Logical indicating whether to use cleaned expression data
#'     (default is FALSE).
#' @param biexp Logical indicating whether to apply biexponential transformation
#'     (default is FALSE).
#' @param af.spectra Optional autofluorescence spectra to include.
#' @param plot.prefix Optional prefix for plot titles (default is "Initial").
#'
#' @return A matrix with the fluorophore spectra.
#' @export


get.fluorophore.spectra <- function( flow.control, asp, use.clean.expr = FALSE,
                                         biexp = FALSE, af.spectra = NULL,
                                         plot.prefix = "Initial" )
{
  spectra.zero <- rep( 0, flow.control$spectral.channel.n )
  names( spectra.zero ) <- flow.control$spectral.channel

  if ( !is.null( plot.prefix ) ) {
    plot.title <- paste( plot.prefix, asp$spectra.file.name )
  } else {
    plot.title <- asp$spectra.file.name
  }

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

  # iterate only over non-negative samples
  fluorophore.samples <- flow.control$fluorophore[ ! grepl( "negative",
                                                          flow.control$fluorophore,
                                                          ignore.case = TRUE ) ]

  fluorophore.channels <- flow.control$channel[ ! grepl( "negative",
                                                       flow.control$fluorophore,
                                                       ignore.case = TRUE ) ]

  if ( !use.clean.expr ) {
    fluorophore.event.samples <- flow.control$event.sample[ flow.control$event.sample %in%
                                                             fluorophore.samples ]
    fluorophore.event.samples <- droplevels( fluorophore.event.samples )

    expr.data <- flow.control$expr.data[ flow.control$event.sample %in%
                                           fluorophore.event.samples, ]

    if ( biexp )
      expr.data <- apply( expr.data, 2, biexp.transform )

    marker.spectra <- lapply( fluorophore.samples, function( samp ) {

      message( paste( "\033[32m", "Processing", samp, "\033[0m" ) )

      peak.channel <- fluorophore.channels[ fluorophore.samples == samp ]

      peak.channel.expr <- expr.data[
        which( fluorophore.event.samples == samp ),
        peak.channel ]

      fluor.spectra.coef <- spectra.zero

      for ( channel in flow.control$spectral.channel ) {
        if ( channel == peak.channel ) {
          fluor.spectra.coef[ channel ] <- 1.0
        } else {
          channel.expr <- expr.data[
            which( fluorophore.event.samples == samp ),
            channel ]

          # fit robust linear model
          spectra.model.result <- fit.robust.linear.model(
            peak.channel.expr, channel.expr,
            peak.channel, channel, asp )

          fluor.spectra.coef[ channel ] <- spectra.model.result[ 2 ]
        }
      }

      # normalize fluor.spectra.coef
      fluor.spectra.coef <- fluor.spectra.coef / max( fluor.spectra.coef )

      fluor.spectra.coef

    } )

    marker.spectra <- do.call( rbind, marker.spectra )
    rownames( marker.spectra ) <- fluorophore.samples

  } else {
    fluorophore.event.samples <- flow.control$clean.event.sample[ flow.control$clean.event.sample %in%
                                                                   fluorophore.samples ]
    fluorophore.event.samples <- droplevels( fluorophore.event.samples )

    expr.data <- flow.control$clean.expr[ flow.control$clean.event.sample %in%
                                            fluorophore.event.samples, ]

    if ( biexp ) {
      expr.data <- apply( expr.data, 2, biexp.transform )
    }

    marker.spectra <- lapply( fluorophore.samples, function( samp ) {

      message( paste("\033[32m", "Processing", samp, "\033[0m" ) )

      peak.channel <- fluorophore.channels[ fluorophore.samples == samp ]

      peak.channel.expr <- expr.data[
        which( fluorophore.event.samples == samp ),
        peak.channel ]

      fluor.spectra.coef <- spectra.zero

      for ( channel in flow.control$spectral.channel ) {
        if ( channel == peak.channel ) {
          fluor.spectra.coef[ channel ] <- 1.0
        } else {
          channel.expr <- expr.data[
            which( fluorophore.event.samples == samp ),
            channel ]

          # fit robust linear model
          spectra.model.result <- fit.robust.linear.model(
            peak.channel.expr, channel.expr,
            peak.channel, channel, asp )

          fluor.spectra.coef[ channel ] <- spectra.model.result[ 2 ]
        }
      }

      # normalize fluor.spectra.coef
      fluor.spectra.coef <- fluor.spectra.coef / max( fluor.spectra.coef )

      fluor.spectra.coef

    } )

    marker.spectra <- do.call( rbind, marker.spectra )
    rownames( marker.spectra ) <- fluorophore.samples
  }

  # plot spectra
  if ( asp$figures ) {
    fluorophore.spectra.plot <- marker.spectra

    if ( !is.null( af.spectra ) ) {
      fluorophore.spectra.plot <- rbind( fluorophore.spectra.plot, af.spectra )
    }

    plot.spectra( fluorophore.spectra.plot, flow.control, asp, plot.title,
                 asp$figure.spectra.dir )
  }

  if ( !is.null( asp$table.spectra.dir ) ) {
    write.csv( fluorophore.spectra.plot,
              file = file.path( asp$table.spectra.dir,
                               sprintf( "%s.csv", plot.title ) ) )
  }

  # cosine similarity QC for controls
  if ( asp$figures )
    plot.similarity.matrix( fluorophore.spectra.plot, asp, plot.prefix )

  similarity.matrix <- cosine.similarity( t( fluorophore.spectra.plot ) )

  unique.similarity <- similarity.matrix * lower.tri( similarity.matrix )

  similarity.idx <- which( unique.similarity > asp$similarity.warning.n,
                           arr.ind = TRUE )

  similarity.error <- nrow( similarity.idx ) > 0

  if ( similarity.error ) {
    fluor1 <- rownames( similarity.matrix )[ similarity.idx[ , 1 ] ]
    fluor2 <- colnames( similarity.matrix )[ similarity.idx[ , 2 ] ]
    similarity.values <- similarity.matrix[ similarity.idx ]

    similarity.qc <- data.frame( Fluor1 = fluor1, Fluor2 = fluor2,
                                Similarity = similarity.values )

    print( similarity.qc )

    message( "\033[31m Similarity over 0.95 detected for one or more pairs of fluorophores.

    Check the table below for problematic combinations.
    If both Fluor1 and Fluor2 are fluorophores,
    manually inspect the controls to confirm they have been prepared correctly.
    Check the fcs_control_table to be sure you have set it up properly.

    If one of the pair is AF, the other likely has minimal signal.
    In this case, run clean.controls and set use.clean.expr to TRUE.
    If you have already done that, manually inspect the control for real signal.
         \033[0m" )
  }

  marker.spectra
}

