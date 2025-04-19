# unmix_folder.r

#' @title Unmix All FCS Files in a Directory
#'
#' @description This function unmixes all FCS files in a specified directory
#'     using the provided spectra and method, and saves the unmixed data to an
#'     output directory.
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @param fcs.dir Directory containing FCS files to be unmixed.
#' @param raw.channels Vector of raw channel names.
#' @param spectra Matrix containing spectra information.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param method Method to be used for unmixing (default is "ols").
#' @param output.dir Directory to save the unmixed FCS files (default is asp$unmixed.fcs.dir).
#' @param include.raw Logical indicating whether to include raw data in the output.
#' @param include.imaging Logical indicating whether to include imaging data in the output.
#' @param allow.negative Logical indicating whether to allow negative coefficients.
#'
#' @return None. Saves the unmixed FCS files to the specified output directory.
#' @export


unmix.folder <- function( fcs.dir, raw.channels, spectra, asp,
                          method = "ols",
                          output.dir = NULL,
                          include.raw = FALSE,
                          include.imaging = FALSE,
                          allow.negative = TRUE ){

  if( is.null( output.dir ) ){
    output.dir <- asp$unmixed.fcs.dir
  }

  files.to.unmix <- list.files( fcs.dir, pattern = ".fcs", full.names = TRUE )

  # set up parallel processing
  if( asp$parallel ){
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  lapply.function( files.to.unmix, FUN = unmix.fcs, raw.channels, spectra, asp,
                   method, output.dir, include.raw, include.imaging, allow.negative )

}
