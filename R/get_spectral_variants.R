# get_spectral_variants.r

#' @title Get Spectral Variations for Fluorophores
#'
#' @description
#' This function cycles through all the fluorophores defined in `flow.control`,
#' identifying variations in spectral profiles. It does this by performing SOM
#' clustering on the positive events in the cleaned control data. The output is
#' saved as an .rds file, and figures summarizing the variation are saved, if
#' desired.
#'
#' @importFrom EmbedSOM SOM
#'
#' @param flow.control A list containing flow cytometry control data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param spectra A matrix containing the spectral data. Fluorophores in rows,
#' detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param n.cells Numeric, default `2000`. Number of cells to use for defining
#' the variation in spectra. Up to `n.cells` cells will be selected as positive
#' events in the peak channel for each fluorophore, above the `pos.quantile` in
#' the unstained sample.
#' @param pos.quantile Numeric, default `0.995`. Threshold for positivity. This
#' quantile will be used to define the maximum extent of the unstained sample in
#' the unmixed space. Anything above that will be considered positive in the
#' single-stained controls.
#' @param som.dim Numeric, default `10`. Number of x and y dimensions to use in
#' the SOM for clustering the spectral variation.
#' @param sim.threshold Numeric, default `0.98`. Threshold for cosine similarity-
#' based exclusion of spectral variants. Any variant less than `sim.threshold`
#' by `cosine.similarity` from the optimized spectrum for that fluorophore (from
#' `spectra`) will be excluded from output. This helps to exclude autofluorescence
#' contamination.
#' @param plot.dir File path to whether the figures will be saved. Default is
#' `NULL`, in which case `asp$variant.dir` will be used.
#'
#' @return A vector with the indexes of events inside the initial gate.
#'
#' @export

get.spectral.variants <- function( flow.control, asp, spectra, af.spectra,
                                   n.cells = 2000, pos.quantile = 0.995,
                                   som.dim = 10, sim.threshold = 0.98,
                                   plot.dir = NULL ) {

  # check for clean.expr.data
  if ( is.null( flow.control$clean.expr ) )
    stop( "Run clean.controls prior to getting spectral variants" )

  if ( asp$figures ) {
    if ( is.null( plot.dir ) )
      plot.dir <- asp$variant.dir
    if ( !dir.exists( plot.dir ) )
      dir.create( plot.dir )
  }

  fluorophores <- rownames( spectra )[ rownames( spectra ) != "AF" ]
  spectral.channel <- colnames( spectra )
  af.n <- nrow( af.spectra )
  fluorophore.n <- nrow( spectra )
  detector.n <- ncol( spectra )

  # get thresholds for positivity
  # use flow.control$expr.data == "AF"
  # can possibly use per-cell AF unmix, but higher threshold will be faster and impact likely minimal
  unstained.exprs <- flow.control$clean.expr[
    which( flow.control$clean.event.sample == "AF" ), spectral.channel ]

  if ( nrow( unstained.exprs ) > n.cells )
    unstained.exprs <- unstained.exprs[ 1:n.cells, ]

  background <- apply( unstained.exprs, 2, median )

  unstained.unmixed <- unmix.ols( unstained.exprs, spectra )
  pos.thresholds <- apply( unstained.unmixed[ , fluorophores ], 2, function( col )
    quantile( col, pos.quantile ) )

  # lapply over fluors in flow.control$clean.expr
  spectral.variants <- lapply( fluorophores, function( fluor ) {
    if ( asp$verbose )
      message( paste( "Getting spectral variants for", fluor ) )

    peak.channel <- flow.control$channel[[ fluor ]]

    # get positive events
    pos.event <- flow.control$clean.expr[
      which( flow.control$clean.event.sample == fluor ), spectral.channel ]

    pos.event <- pos.event[ which( pos.event[ , peak.channel ] > pos.thresholds[ fluor ] ), ]

    # som cluster
    set.seed( asp$variant.seed )
    map <- EmbedSOM::SOM( pos.event,
                          xdim = som.dim, ydim = som.dim )

    # subtract background
    variant.spectra <- sweep( map$codes, 2, background, FUN = "-" )
    # normalize
    variant.spectra <- t( apply( variant.spectra, 1, function( x ) x / max( x ) ) )
    variant.spectra <- as.matrix( na.omit( variant.spectra ) )

    rownames( variant.spectra ) <- paste0( fluor, "_", 1:nrow( variant.spectra ) )

    # qc to remove dissimilar spectral variants (usually AF contamination)
    original.spectrum <- spectra[ ( rownames( spectra ) == fluor ), ]

    similar <- sapply( seq_len( nrow( variant.spectra ) ), function( sp ) {
      sim <- cosine.similarity( rbind( original.spectrum, variant.spectra[ sp, ] ) )
      sim <- sim[ lower.tri( sim ) ]
      sim > sim.threshold
    } )

    variant.spectra <- variant.spectra[ similar, ]

    if ( asp$figures ) {
      spectral.variant.plot( variant.spectra, original.spectrum,
                             title = paste0( fluor, "_variants" ),
                             save = TRUE,
                             plot.dir = plot.dir )
    }

    return( variant.spectra )

  } )

  names( spectral.variants ) <- fluorophores

  variants <- list(
    thresholds = pos.thresholds,
    variants = spectral.variants
  )

  saveRDS( variants, file = file.path( plot.dir, asp$variant.filename ) )

  return( variants )

}
