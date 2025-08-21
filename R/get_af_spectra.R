# get_af_spectra.r

#' @title Get Autofluorescence Spectra
#'
#' @description
#' Extracts autofluorescence spectra from an unstained samples. Intended for use
#' with `unmix.autospectral`. Uses FlowSOM (EmbedSOM) clustering for rapid
#' identification of cells with similar AF profiles.
#'
#' @importFrom FlowSOM SOM
#' @importFrom stats cutree dist hclust as.dist
#' @importFrom flowCore read.FCS exprs
#' @importFrom flowWorkspace flowjo_biexp
#'
#' @param unstained.sample Path and file name for a unstained sample FCS file.
#' The sample type and processing (protocol) method should match the fully
#' stained samples to which the AF will be applied, ideally.
#' @param asp The AutoSpectral parameter list.
#' @param spectral.channels The fluorescence detector channels that will be used
#' to identify the autofluorescence. These should match the channels used for
#' subsequent unmixing and should normally be all fluorescent area measurements.
#' Default is `flow.control$spectral.channel`.
#' @param similarity Numeric between 0 and 1, default is `0.99`. The cosine
#' similarity value that will trigger merging of FlowSOM nodes to create
#' metaclusters that define the autofluorescence signatures.
#' @param max.spectra Maximum number of AF signatures to return. Default is `20`.
#' @param min.spectra Minimum number of AF signatures to return. Default is `5`.
#' @param biexp Logical, whether to use a biexponential transform on the data
#' prior to clustering. Default is `TRUE`.
#' @param figures Logical, whether to plot the spectral traces and heatmap for
#' the AF signatures. Default is `TRUE`.
#' @param plot.dir Directory (folder) where the plots will be saved. Default is
#' `NULL`, which inherits from `asp$figure.af.dir`.
#' @param title Title for the output spectral plots. Default is
#' `Autofluorescence spectra`.
#'
#' @return A matrix of autofluorescence spectra.
#'
#' @export

get.af.spectra <- function( unstained.sample, asp,
                            spectral.channels = flow.control$spectral.channel,
                            similarity = 0.99,
                            max.spectra = 20, min.spectra = 5,
                            biexp = TRUE,
                            figures = TRUE, plot.dir = NULL,
                            title = "Autofluorescence spectra" ) {

  # set defaults
  if ( is.null( plot.dir ) )
    plot.dir <- asp$figure.af.dir

  stopifnot( similarity > 0, similarity < 1 )

  # import unstained sample
  if ( asp$verbose )
    message( paste( "Reading FCS file", unstained.sample ) )

  unstained.ff <- suppressWarnings( flowCore::read.FCS( unstained.sample,
                                                        transformation = FALSE,
                                      truncate_max_range = FALSE, emptyValue = FALSE ) )

  unstained.exprs <- flowCore::exprs( unstained.ff )[ , spectral.channels ]

  if ( biexp ) {
    biexp.transform <- flowWorkspace::flowjo_biexp(
      channelRange = asp$default.transformation.param$length,
      maxValue = asp$default.transformation.param$max.range,
      pos = asp$default.transformation.param$pos,
      neg = asp$default.transformation.param$neg,
      widthBasis = asp$default.transformation.param$width,
      inverse = FALSE
    )

    unstained.trans <- apply( unstained.exprs, 2, biexp.transform )
  } else {
    unstained.trans <- unstained.exprs
  }

  # get cluster of AF from unstained unmixed
  if ( asp$verbose )
    message( "Creating a self-organizing map of the autofluorescence" )

  set.seed( 42 )
  fsom <- FlowSOM::SOM( unstained.trans )

  # get median of each FlowSOM node to obtain spectra
  cluster.spectra <- sapply( unique( fsom$mapping[ , 1 ] ), function( cl ) {
    clust.idx <- which( fsom$mapping[ , 1 ] == cl )
    cells <- unstained.exprs[ clust.idx, , drop = FALSE ]
    apply( cells, 2, median )
  } )
  cluster.spectra <- t( apply( cluster.spectra, 2, function( x ) x / max( x ) ) )

  # use cosine similarity of AF profiles to define metaclusters
  cluster.sim <- cosine.similarity( cluster.spectra )
  cluster.dist <- as.dist( 1 - cluster.sim )

  n.clusters <- 0
  max.iter <- 10
  iter <- 0
  increment <- 0.005

  while( ( n.clusters < min.spectra || n.clusters > max.spectra ) && iter < max.iter ) {
    metaclusters <- cutree( hclust( cluster.dist, method = "complete" ),
                            h = ( 1 - similarity ) )
    n.clusters <- length( unique( metaclusters ) )

    if ( n.clusters > max.spectra ) {
      message( paste0( "More than ", max.spectra,
                      " clusters found using similarity = ", similarity,
                      ". Trying ", similarity - increment ) )
      similarity <- similarity - increment
    } else if ( n.clusters < min.spectra ) {
      message( paste0( "Fewer than ", min.spectra,
                      " clusters found using similarity = ", similarity,
                      ". Trying ", similarity + increment ) )
      similarity <- similarity + increment
    }
    iter <- iter + 1
  }

  names( metaclusters ) <- NULL

  # get median of each cluster to obtain AF spectra
  af.spectra <- sapply( unique( metaclusters ), function( cl ) {
    clust.idx <- which( fsom$mapping[ , 1 ] == cl )
    cells <- unstained.exprs[ clust.idx, , drop = FALSE ]
    apply( cells, 2, median )
  } )

  af.spectra <- t( apply( af.spectra, 2, function( x ) x / max( x ) ) )
  rownames( af.spectra ) <- paste0( "AF", 1:nrow( af.spectra ) )

  # save as CSV
  af.file.name <- paste0( asp$af.file.name, ".csv" )
  write.csv( af.spectra, file = file.path( asp$table.spectra.dir, af.file.name ) )

  if ( figures ) {
    if ( asp$verbose )
      message( "Plotting spectra" )

    if ( !dir.exists( plot.dir ) )
      dir.create( plot.dir )

    spectral.trace( af.spectra, title = title,
                    plot.dir = plot.dir, split.lasers = FALSE )

    spectral.heatmap( af.spectra, title = title,
                      output.dir = plot.dir )
  }

  return( af.spectra )
}
