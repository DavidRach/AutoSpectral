# unmix_fcs.r

#' @title Unmix FCS Data
#'
#' @description
#' This function performs spectral unmixing on FCS data using various methods.
#'
#' @importFrom flowCore read.FCS keyword exprs flowFrame parameters pData
#' @importFrom flowCore write.FCS parameters<-
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom utils packageVersion
#'
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectra A matrix containing the spectral data.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method A character string specifying the unmixing method to use.
#' Options are `OLS`, `WLS`, `Poisson` and `FastPoisson`. Default is `OLS`.
#' `FastPoisson` requires installation of `AutoSpectralRcpp`
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `NULL`, in which case weighting will be done by
#' channel means (Poisson variance). Only used for `WLS`.
#' @param output.dir A character string specifying the directory to save the
#' unmixed FCS file. Default is `NULL`.
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`.
#' @param include.raw A logical value indicating whether to include raw
#' expression data in the written FCS file. Default is `FALSE`.
#' @param include.imaging A logical value indicating whether to include imaging
#' parameters in the written FCS file. Default is `FALSE`.
#' @param divergence.threshold Numeric. Used for `FastPoisson` only.
#' Threshold to trigger reversion towards WLS unmixing when Poisson result
#' diverges for a given point.
#' @param divergence.handling String. How to handle divergent cells from Poisson
#' IRLS. Options are `NonNeg` (non-negativity will be enforced), `WLS` (revert
#' to WLS intial unmix) or `Balance` (WLS and NonNeg will be averaged).
#' Default is `Balance`
#' @param balance.weight Numeric. Weighting to average non-convergent cells.
#' Used for `Balance` option under `divergence.handling`. Default is `0.5`.
#'
#' @return None. The function writes the unmixed FCS data to a file.
#'
#' @export

unmix.fcs <- function( fcs.file, spectra, asp, flow.control,
                       method = "OLS",
                       weights = NULL,
                       output.dir = NULL,
                       file.suffix = NULL,
                       include.raw = FALSE,
                       include.imaging = FALSE,
                       divergence.threshold = 1e4,
                       divergence.handling = "Balance",
                       balance.weight = 0.5 ){

  if ( is.null( output.dir ) ){
    output.dir <- asp$unmixed.fcs.dir
  }

  # import fcs, without warnings for fcs 3.2
  fcs.data <- suppressWarnings(
    read.FCS( fcs.file, transformation = FALSE,
              truncate_max_range = FALSE, emptyValue = FALSE )
  )

  fcs.keywords <- keyword( fcs.data )
  file.name <- keyword( fcs.data, "$FIL" )

  # deal with manufacturer peculiarities in writing fcs files
  if ( asp$cytometer == "ID7000" ) {
    file.name <- sub( "Raw", method, file.name )

  } else if ( asp$cytometer == "DiscoverS8" | asp$cytometer == "DiscoverA8" ) {
    file.name <- fcs.keywords$FILENAME
    file.name <- sub( ".*\\/", "", file.name )
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )

  } else {
    file.name <- sub( ".fcs", paste0( " ", method, ".fcs" ), file.name )
  }

  if ( !is.null( file.suffix ) )
    file.name <- sub( ".fcs", paste0( " ", file.suffix, ".fcs" ), file.name )

  # extract exprs
  fcs.exprs <- flowCore::exprs( fcs.data )
  rm( fcs.data )

  spectral.exprs <- fcs.exprs[ , flow.control$spectral.channel, drop = FALSE ]

  other.channels <- setdiff( colnames( fcs.exprs ), flow.control$spectral.channel )
  other.exprs <- fcs.exprs[ , other.channels, drop = FALSE ]

  # apply unmixing using selected method
  unmixed.data <- switch( method,
                         "OLS" = unmix.ols( spectral.exprs, spectra ),
                         "WLS" = unmix.wls( spectral.exprs, spectra, weights ),
                         "Poisson" = unmix.poisson( spectral.exprs, spectra, asp ),
                         "FastPoisson" = {
                           if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE ) &&
                               "unmix.poisson.fast" %in% ls( getNamespace( "AutoSpectralRcpp" ) ) ) {
                             tryCatch(
                               AutoSpectralRcpp::unmix.poisson.fast( spectral.exprs,
                                                                     spectra,
                                                                     maxit = asp$rlm.iter.max,
                                                                     tol = 1e-6,
                                                                     n_threads = asp$worker.process.n,
                                                                     divergence.threshold = divergence.threshold,
                                                                     divergence.handling = divergence.handling,
                                                                     balance.weight = balance.weight ),
                               error = function( e ) {
                                 warning( "FastPoisson failed, falling back to standard Poisson: ", e$message )
                                 unmix.poisson( spectral.exprs, spectra, asp )
                               }
                             )
                           } else {
                             warning( "AutoSpectralRcpp not available, falling back to standard Poisson." )
                             unmix.poisson( spectral.exprs, spectra, asp )
                           }
                         },
                         stop( "Unknown method" )
  )

  # remove imaging parameters, which are pretty useless in big panels
  if ( asp$cytometer == "DiscoverS8" | asp$cytometer == "DiscoverA8" &
       !include.imaging ) {
    other.exprs <- other.exprs[ , asp$time.and.scatter ]
  }

  if ( include.raw ){
    # add back raw exprs and others
    unmixed.data <- cbind( fcs.exprs, unmixed.data )
  } else {
    # add back other columns
    unmixed.data <- cbind( other.exprs, unmixed.data )
  }

  rm( spectral.exprs, other.exprs )

  # define new fcs file
  flow.frame <- flowFrame( unmixed.data )
  params <- pData( parameters( flow.frame) )
  params$maxRange[ params$name != "Time" ] <- asp$expr.data.max

  # add antigen labels based on match between param name and fluorophore
  fluor.match.idx <- match( params$name, flow.control$fluorophore )

  params$desc[ !is.na( fluor.match.idx ) ] <-
    flow.control$antigen[ fluor.match.idx[ !is.na( fluor.match.idx ) ] ]

  other.match.idx <- match( params$name, other.channels )

  params$desc[ !is.na( other.match.idx ) ] <- NA

  parameters( flow.frame ) <- AnnotatedDataFrame( params )

  # update keywords
  #fcs.keywords[[ "$FIL" ]] <- file.name
  #fcs.keywords[[ "$UNMIXINGMETHOD" ]] <- method
  #fcs.keywords[[ "$SPECTRA" ]] <- table( spectra )
  #fcs.keywords[[ "$AUTOSPECTRAL" ]] <- packageVersion( "AutoSpectral" )
  keyword( flow.frame ) <- fcs.keywords
  keyword( flow.frame )[[ "$FIL" ]] <- file.name
  keyword( flow.frame )[[ "$UNMIXINGMETHOD" ]] <- method
  keyword( flow.frame )[[ "$AUTOSPECTRAL" ]] <- packageVersion( "AutoSpectral" )

  for (i in seq_along(colnames(unmixed.data))) {
    keyword(flow.frame)[[paste0("$P", i, "N")]] <- colnames(unmixed.data)[i]
    keyword(flow.frame)[[paste0("$P", i, "S")]] <- colnames(unmixed.data)[i]
  }

  suppressWarnings( write.FCS( flow.frame, filename = file.path( output.dir, file.name ) ) )

}
