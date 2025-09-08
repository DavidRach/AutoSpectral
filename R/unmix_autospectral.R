# unmix_autospectral.r

#' @title Unmix AutoSpectral
#'
#' @description
#' Unmix using the AutoSpectral method to extract autofluorescence at the
#' single cell level.
#'
#' @param raw.data Expression data from raw fcs files. Cells in rows and
#' detectors in columns. Columns should be fluorescent data only and must
#' match the columns in spectra.
#' @param spectra Spectral signatures of fluorophores, normalized between 0
#' and 1, with fluorophores in rows and detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param weighted Logical, whether to use ordinary or weighted least squares
#' unmixing as the base algorithm. Default is `FALSE` and will use OLS.
#' @param weights Optional numeric vector of weights (one per fluorescent
#' detector). Default is `NULL`, in which case weighting will be done by
#' channel means (Poisson variance). Only used if `weighted`.
#' @param calculate.error Logical, whether to calculate the RMSE unmixing model
#' accuracy and include it as an output.
#' @param use.dist0 Logical, controls whether the selection of the optimal AF
#' signature for each cell is determined by which unmixing brings the cell
#' closest to 0 (`use.dist0` = `TRUE`) or by which unmixing minimizes the
#' per-cell residual (`use.dist0` = `FALSE`). Default is `FALSE`.
#'
#' @return Unmixed data with cells in rows and fluorophores in columns.
#'
#' @export

unmix.autospectral <- function( raw.data, spectra, af.spectra,
                                weighted = FALSE, weights = NULL,
                                calculate.error = FALSE,
                                use.dist0 = FALSE ) {

  # check for AF in spectra, remove if present
  if ( "AF" %in% rownames( spectra ) ) {
    message( "Removing default AF channel" )
    spectra <- spectra[ rownames( spectra ) != "AF", , drop = FALSE ]
  }

  if ( is.null( af.spectra ) )
    stop( "Multiple AF spectra must be provided." )

  fluorophores <- rownames( spectra )
  af.n <- nrow( af.spectra )
  fluorophore.n <- nrow( spectra )
  detector.n <- ncol( spectra )

  # initial no AF unmixing
  if ( weighted )
    no.af.unmixed <- unmix.wls( raw.data, spectra, weights )
  else
    no.af.unmixed <- unmix.ols( raw.data, spectra )

  no.af.residual <- rowSums( ( raw.data - ( no.af.unmixed %*% spectra ) )^2 )
  no.af.unmixed <- cbind( no.af.unmixed, no.af.residual/1e3 )

  # unmix for each af.spectrum
  message( "Extracting AF cell-by-cell" )

  model.unmixings <- vector( "list", length = af.n + 1 )
  model.residuals <- vector( "list", length = af.n + 1 )

  combined.spectra <- matrix( NA_real_, nrow = fluorophore.n + 1, ncol = detector.n )
  colnames( combined.spectra ) <- colnames( spectra )
  rownames( combined.spectra ) <- c( fluorophores, "AF" )

  for ( af in seq_len( af.n ) ) {
    combined.spectra[ 1:fluorophore.n, ] <- spectra
    combined.spectra[ fluorophore.n + 1, ] <- af.spectra[ af, , drop = FALSE ]

    if ( weighted )
      unmixed <- unmix.wls( raw.data, combined.spectra, weights )
    else
      unmixed <- unmix.ols( raw.data, combined.spectra )

    residual <- rowSums( ( raw.data - ( unmixed %*% combined.spectra ) )^2 )

    model.unmixings[[ af ]] <- unmixed
    model.residuals[[ af ]] <- residual
  }

  model.unmixings[[ af.n + 1 ]] <- no.af.unmixed
  model.residuals[[ af.n + 1 ]] <- no.af.residual

  model.residuals <- do.call( cbind, model.residuals )

  # determine best model for each cell
  if ( use.dist0 ) {
    model.dist0 <- vapply( seq_along( model.unmixings ), function( i ) {
      non.af <- model.unmixings[[ i ]][ , fluorophores, drop = FALSE ]
      rowSums( abs( non.af ) )
    }, FUN.VALUE = numeric( nrow( model.unmixings[[ 1 ]] ) ) )

    af.idx <- apply( model.dist0, 1, which.min )

  } else {
    af.idx <- apply( model.residuals, 1, which.min )
  }

  # merge
  unmixed.data <- t( vapply( seq_along( af.idx ), function( i ) {
    model.unmixings[[ af.idx[ i ] ]][ i, ]
  }, numeric( ncol( model.unmixings[[ 1 ]] ) ) ) )

  rm( model.unmixings )

  # add per-cell index of which AF has been used
  unmixed.data <- cbind( unmixed.data, af.idx )

  # calculate average per cell error
  if ( calculate.error ) {
    residual <- model.residuals[ cbind( 1:nrow( model.residuals ), af.idx ) ]
    RMSE <- sqrt( mean( residual ) )
  }

  rm( model.residuals )

  colnames( unmixed.data ) <- c( fluorophores, "AF", "AF Index" )

  # set value of AF Index to 0 for cells best fit by spectra without AF
  #unmixed.data[ , "AF Index" == ( af.n + 1 ) ] <- 0

  if ( calculate.error )
    return( list(
      RMSE = RMSE,
      unmixed.data = unmixed.data
    ) )
  else
    return( unmixed.data )

}
