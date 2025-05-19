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
#' @param spectra Matrix containing spectra information.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method Method to be used for unmixing (default is "OLS").
#' @param output.dir Directory to save the unmixed FCS files (default is asp$unmixed.fcs.dir).
#' @param file.suffix A character string to append to the output file name.
#'     Default is NULL.
#' @param include.raw Logical indicating whether to include raw data in the output.
#'     Default is FALSE.
#' @param include.imaging Logical indicating whether to include imaging data in the output.
#'     Default is FALSE.
#'
#' @return None. Saves the unmixed FCS files to the specified output directory.
#' @export


unmix.folder <- function( fcs.dir, spectra, asp, flow.control,
                          method = "OLS",
                          output.dir = NULL,
                          file.suffix = NULL,
                          include.raw = FALSE,
                          include.imaging = FALSE ){

  if ( is.null( output.dir ) ){
    output.dir <- asp$unmixed.fcs.dir
  }

  files.to.unmix <- list.files( fcs.dir, pattern = ".fcs", full.names = TRUE )

  # set up parallel processing
  if ( asp$parallel & method != "Poisson" ){
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  lapply.function( files.to.unmix, FUN = unmix.fcs, spectra, asp, flow.control,
                   method, output.dir, file.suffix, include.raw,
                   include.imaging, allow.negative )

}
