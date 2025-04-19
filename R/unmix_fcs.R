# unmix_fcs.r

#' @title Unmix FCS Data
#' @description This function performs spectral unmixing on FCS data using
#'     various methods.
#' @importFrom flowCore read.FCS keyword exprs flowFrame parameters pData
#' @importFrom flowCore write.FCS parameters<-
#' @importFrom Biobase AnnotatedDataFrame
#' @param fcs.file A character string specifying the path to the FCS file.
#' @param spectra A matrix containing the spectral data.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method A character string specifying the unmixing method to use.
#'     Options are "ols", "wls", and "poisson". Default is "ols".
#' @param output.dir A character string specifying the directory to save the
#'     unmixed FCS file. Default is NULL.
#' @param file.suffix A character string to append to the output file name.
#'     Default is NULL.
#' @param include.raw A logical value indicating whether to include raw expression
#'     data in the output. Default is FALSE.
#' @param include.imaging A logical value indicating whether to include imaging
#'     parameters in the output. Default is FALSE.
#' @param allow.negative A logical value indicating whether to allow negative
#'     values in the unmixing process. Default is TRUE.
#' @return None. The function writes the unmixed FCS data to a file.
#' @export


unmix.fcs <- function( fcs.file, spectra, asp, flow.control, method = "ols",
                       output.dir = NULL, file.suffix = NULL,
                       include.raw = FALSE,
                       include.imaging = FALSE,
                       allow.negative = TRUE ){

  if( is.null( output.dir ) ){
    output.dir <- asp$unmixed.fcs.dir
  }

  # import fcs, without warnings for fcs 3.2
  fcs.data <- suppressWarnings(
    read.FCS( fcs.file, transformation = FALSE,
              truncate_max_range = FALSE )
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
  fcs.exprs <- exprs( fcs.data )
  rm( fcs.data )
  spectral.exprs <- fcs.exprs[ , flow.control$spectral.channel, drop = FALSE ]

  other.channels <- setdiff( colnames( fcs.exprs ), flow.control$spectral.channel )
  other.exprs <- fcs.exprs[ , other.channels, drop = FALSE ]

  # apply unmixing using selected method
  unmixed.data <- switch( method,
                         "ols" = unmix.ols( spectral.exprs, spectra ),
                         "wls" = unmix.wls( spectral.exprs, spectra ),
                         "poisson" = unmix.poisson( spectral.exprs, spectra,
                                                    asp, allow.negative ),
                         stop( "Unknown method" )
  )

  # remove imaging parameters, which are pretty useless in big panels
  if ( asp$cytometer == "DiscoverS8" | asp$cytometer == "DiscoverA8" &
       !include.imaging ) {
    other.exprs <- other.exprs[ , asp$time.and.scatter ]
  }

  if( include.raw ){
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

  write.FCS( flow.frame, filename = file.path( output.dir, file.name ) )

}
