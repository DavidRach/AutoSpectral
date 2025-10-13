# get_fluor_variants.r

#' @title Get Fluorophore Variants
#'
#' @description
#' Assesses variation in the spectral signature of a single-stained flow
#' cytometry control sample. Uses SOM-based clustering on the brightest positive
#' events in the file.
#'
#' @importFrom flowCore read.FCS exprs
#'
#' @param fluor The name of the fluorophore.
#' @param file.name A named vector of file names for the samples.
#' @param control.dir The directory containing the control files.
#' @param spectra A matrix containing the spectral data. Fluorophores in rows,
#' detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param spectral.channel A vector of spectral channels.
#' @param universal.negative A named vector of unstained negative samples, with
#' names corresponding to the fluorophores.
#' @param control.type Character, either "beads" or "cells". Determines the type
#' of control sample being used and the subsequent processing steps.
#' @param raw.thresholds A named vector of numerical values corresponding to
#' the threshold for positivity in each raw detector channel. Determined by the
#' 99.5th percentile on the unstained sample, typically.
#' @param unmixed.thresholds A named vector of numerical values corresponding to
#' the threshold for positivity in each unmixed channel. Determined by the
#' 99.5th percentile on the unstained sample, typically after single-cell AF
#' unmixing.
#' @param flow.channel A named vector of peak raw channels, one per fluorophore.
#' @param som.dim Numeric, default `10`. Number of x and y dimensions to use in
#' the SOM for clustering the spectral variation.
#' @param n.cells Numeric, default `2000`. Number of cells to use for defining
#' the variation in spectra. Up to `n.cells` cells will be selected as positive
#' events in the peak channel for each fluorophore, above the `pos.quantile` in
#' the unstained sample.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#' @param output.dir File path to whether the figures and .rds data file will be
#' saved. Default is `NULL`, in which case `asp$variant.dir` will be used.
#' @param sim.threshold Numeric, default `0.98`. Threshold for cosine similarity-
#' based exclusion of spectral variants. Any variant less than `sim.threshold`
#' by `cosine.similarity` from the optimized spectrum for that fluorophore (from
#' `spectra`) will be excluded from output. This helps to exclude autofluorescence
#' contamination.
#' @param figures Logical, controls whether the variation in spectra for each
#' fluorophore is plotted in `output.dir`. Default is `TRUE`.
#'
#' @return A matrix with the flow expression data.

get.fluor.variants <- function( fluor, file.name, control.dir, spectra, af.spectra,
                                spectral.channel, universal.negative, control.type,
                                raw.thresholds, unmixed.thresholds, flow.channel,
                                som.dim, n.cells, asp, verbose, output.dir,
                                sim.threshold, figures ) {

  if ( verbose )
    message( paste( "\033[34m", "Getting spectral variants for", fluor, "\033[0m" ) )

  pos.data <- suppressWarnings(
    flowCore::read.FCS( file.path( control.dir, file.name[ fluor ] ),
                        transformation = NULL,
                        truncate_max_range = FALSE,
                        emptyValue = FALSE ) )

  # read exprs for spectral channels only
  pos.data <- flowCore::exprs( pos.data )[ , spectral.channel ]

  # get data above threshold in peak channel
  # restrict to top n events
  peak.channel <- flow.channel[ fluor ]
  raw.idx <- which( pos.data[ , peak.channel ] > raw.thresholds[ peak.channel ] )
  neg.idx <- setdiff( seq_len( nrow( pos.data ) ), raw.idx )

  if ( length( raw.idx ) > n.cells * 2 ) {
    sorted.idx <- order( pos.data[ raw.idx, peak.channel ],
                         decreasing = TRUE )[ 1:( n.cells * 2 ) ]
    raw.idx <- raw.idx[ sorted.idx ]

  }

  if ( "AF" %in% rownames( spectra ) )
    no.af.spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]
  else
    no.af.spectra <- spectra

  if ( control.type[ fluor ] == "cells" ) {

    # autospectral per-cell AF unmixing
    pos.unmixed <- unmix.ols( pos.data[ raw.idx, ], no.af.spectra )

    fluorophores <- rownames( no.af.spectra )
    af.n <- nrow( af.spectra )
    fluorophore.n <- nrow( no.af.spectra )
    detector.n <- ncol( no.af.spectra )
    combined.spectra <- matrix( NA_real_, nrow = fluorophore.n + 1, ncol = detector.n )
    colnames( combined.spectra ) <- colnames( no.af.spectra )
    fluors.af <- c( fluorophores, "AF" )
    rownames( combined.spectra ) <- fluors.af
    combined.spectra[ 1:fluorophore.n, ] <- no.af.spectra
    initial.af <- matrix( 0, nrow = length( raw.idx ), ncol = 1 )
    colnames( initial.af ) <- c( "AF" )
    pos.unmixed <- cbind( pos.unmixed, initial.af )
    fitted.af <- matrix( 0, nrow = length( raw.idx ), ncol = ncol( no.af.spectra ) )
    error <- rowSums( abs( pos.unmixed[ , fluorophores, drop = FALSE ] ) )

    for ( af in seq_len( af.n ) ) {
      combined.spectra[ fluorophore.n + 1, ] <- af.spectra[ af, , drop = FALSE ]
      unmixed.af <- unmix.ols( pos.data[ raw.idx, ], combined.spectra )

      error.af <- rowSums( abs( unmixed.af[ , fluorophores, drop = FALSE ] ) )
      improved <- which( error.af < error )

      error[ improved ] <- error.af[ improved ]
      pos.unmixed[ improved, fluors.af ] <- unmixed.af[ improved, ]

      fitted.af[ improved, ] <- unmixed.af[ improved, "AF", drop = FALSE ] %*%
        af.spectra[ af, , drop = FALSE ]
    }

    # check for data above threshold in unmixed fluor channel
    pos.idx <- which( pos.unmixed[ , fluor ] > unmixed.thresholds[ fluor ]*2 )

    # subtract fitted af component from raw data
    remaining.raw <- pos.data[ raw.idx[ pos.idx ], ] - fitted.af[ pos.idx, ]

    # check that we still have data; if not, return original spectrum
    if ( length( pos.idx ) < asp$min.cell.stop.n )
      return( spectra[ fluor, , drop = FALSE ] )

    # cluster
    som.input <- cbind( pos.unmixed[ pos.idx, ], remaining.raw )
    set.seed( asp$variant.seed )
    map <- EmbedSOM::SOM( som.input, xdim = som.dim, ydim = som.dim )

    # get spectra
    variant.spectra <- t( apply( map$codes[ , spectral.channel ], 1,
                                 function( x ) x / max( x ) ) )
    variant.spectra <- as.matrix( na.omit( variant.spectra ) )
    rownames( variant.spectra ) <- paste0( fluor, "_", 1:nrow( variant.spectra ) )

  } else { # get bead control data and spectral variation
    # get negative background to subtract
    # check for universal negative, if none, use internal negative
    if ( universal.negative[ fluor ] != FALSE ) {
      neg.data <- suppressWarnings(
        flowCore::read.FCS( file.path( control.dir, universal.negative[ fluor ] ),
                  transformation = NULL,
                  truncate_max_range = FALSE,
                  emptyValue = FALSE ) )
      neg.data <- flowCore::exprs( neg.data )[ , spectral.channel ]

      # get background on up to 10k events
      if ( nrow( neg.data ) > asp$gate.downsample.n.beads ) {
        set.seed( asp$variant.seed )
        neg.idx <- sample( nrow( neg.data ), asp$gate.downsample.n.beads )
        background <- apply( neg.data[ neg.idx, spectral.channel ], 2, median )
      } else {
        background <- apply( neg.data[ , spectral.channel ], 2, median )
      }

    } else {
      # select events below unstained raw threshold

      if ( length( neg.idx ) > asp$gate.downsample.n.beads ) {
        # downsample if lots of events
        set.seed( asp$variant.seed )
        neg.idx <- sample( neg.idx, asp$gate.downsample.n.beads )
        background <- apply( pos.data[ neg.idx, spectral.channel ], 2, median )

      } else if ( length( neg.idx ) > asp$min.cell.warning.n ) {
        # use selected data below threshold if moderate numbers of events
        background <- apply( pos.data[ neg.idx, spectral.channel ], 2, median )
      } else if ( nrow( pos.data ) < asp$min.cell.stop.n ) {
        warning( paste0( "Minimal data present in sample: ", fluor,
                         "Variation assessment not possible for this fluorophore." ) )
        return( spectra[ fluor, , drop = FALSE ] )
      } else {
        # low event sample--take lower half of distribution
        neg.selected <- sort( pos.data[ , peak.channel ],
                              decreasing = FALSE )[ 1:nrow( pos.data ) / 2 ]
        background <- apply( neg.selected, 2, median )
      }
    }

    # unmix without AF
    pos.unmixed <- unmix.ols( pos.data[ raw.idx, ], no.af.spectra )

    # cluster
    som.input <- cbind( pos.unmixed, pos.data[ raw.idx, ] )
    set.seed( asp$variant.seed )
    map <- EmbedSOM::SOM( som.input, xdim = som.dim, ydim = som.dim )

    # get spectra, subtracting background from unstained/negative
    variant.spectra <- sweep( map$codes[ , spectral.channel ], 2, background, FUN = "-" )
    variant.spectra <- t( apply( variant.spectra, 1, function( x ) x / max( x ) ) )
    variant.spectra <- as.matrix( na.omit( variant.spectra ) )
    rownames( variant.spectra ) <- paste0( fluor, "_", 1:nrow( variant.spectra ) )
  }

  # qc to remove dissimilar spectral variants (usually AF contamination)
  original.spectrum <- spectra[ fluor, ]

  similar <- sapply( seq_len( nrow( variant.spectra ) ), function( sp ) {
    sim <- cosine.similarity( rbind( original.spectrum, variant.spectra[ sp, ] ) )
    sim <- sim[ lower.tri( sim ) ]
    sim > sim.threshold
  } )

  if ( ! any( similar ) )
    variant.spectra <- spectra[ fluor, , drop = FALSE ]
  else
    variant.spectra <- variant.spectra[ similar, , drop = FALSE ]

  if ( figures )
    spectral.variant.plot( variant.spectra, original.spectrum,
                           title = paste0( fluor, "_variants" ),
                           save = TRUE,
                           plot.dir = output.dir )

  return( variant.spectra )
}

