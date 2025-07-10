# get_af_spectra.r

#' @title Get Autofluorescence Spectra
#'
#' @description
#' Extracts autofluorescence spectra from an unstained samples. Intended for use
#' with `unmix.autospectral`. Uses FlowSOM (EmbedSOM) clustering for rapid
#' identification of cells with similar AF profiles.
#'
#' @importFrom EmbedSOM SOM
#' @importFrom stats cutree dist hclust
#' @importFrom flowCore read.FCS exprs
#'
#' @param unstained.sample Path and file name for a unstained sample FCS file.
#' The sample type and processing (protocol) method should match the fully
#' stained samples to which the AF will be applied, ideally.
#' @param spectral.channels The fluorescence detector channels that will be used
#' to identify the autofluorescence. These should match the channels used for
#' subsequent unmixing and should normally be all fluorescent area measurements.
#' Default is `flow.control$spectral.channel`
#' @param threads Numeric. Number of threads to use for parallel processing in
#' the creation of the SOM. Default is `asp$worker.process.n`.
#' @param n.spectra Number of AF signatures to look for. Default is `20`.
#' @param figures Logical, whether to plot the spectral traces and heatmap for
#' the AF signatures. Default is `TRUE`.
#' @param plot.dir Directory (folder) where the plots will be saved. Default is
#' `asp$figure.af.dir`.
#' @param plot.title Title for the output spectral plots. Default is
#' `Autofluorescence spectra`.
#'
#' @return A matrix of autofluorescence spectra.
#'
#' @export

get.af.spectra <- function( unstained.sample,
                            spectral.channels = flow.control$spectral.channel,
                            threads = asp$worker.process.n,
                            n.spectra = 20, figures = TRUE,
                            plot.dir = asp$figure.af.dir,
                            plot.title = "Autofluorescence spectra" ) {

  # import unstained sample
  unstained.ff <- flowCore::read.FCS( unstained.sample, transformation = FALSE,
                                      truncate_max_range = FALSE, emptyValue = FALSE )

  unstained.exprs <- flowCore::exprs( unstained.ff )[ , spectral.channels ]

  # get cluster of AF from unstained unmixed
  map <- EmbedSOM::SOM( unstained.exprs, xdim = 24, ydim = 24,
                        batch = TRUE, parallel = TRUE,
                        threads = threads )

  # use hierarchical clustering to define metaclusters
  clusters <- cutree( k = n.spectra, hclust( method = 'average',
                                      dist( map$codes ) ) )[ map$mapping[ , 1 ] ]

  # get median of each cluster to obtain spectra
  af.spectra <- sapply( unique( clusters ), function( cl ) {
    clust.idx <- which( clusters == cl )
    cells <- unstained.exprs[ clust.idx, , drop = FALSE ]
    apply( cells, 2, median )
  } )

  af.spectra <- t( apply( af.spectra, 2, function( x ) x / max( x ) ) )
  rownames( af.spectra ) <- paste0( "AF", 1:nrow( af.spectra ) )

  if ( figures ) {
    message( "Plotting spectra" )
    if ( !dir.exists( plot.dir ) )
      dir.create( plot.dir )

    spectral.trace( af.spectra, plot.title = plot.title,
                    plot.dir = plot.dir, split.lasers = FALSE )

    spectral.heatmap( af.spectra, plot.prefix = plot.title,
                      output.dir = plot.dir )
  }

  return( af.spectra )
}
