# run_af_removal.r

#' @title Run Autofluorescence Removal
#'
#' @description
#' This function runs the autofluorescence removal process on a list of samples,
#' using the specified parameters and settings.
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @param clean.expr List containing cleaned expression data.
#' @param af.removal.sample Vector of sample names for which autofluorescence
#' removal is to be performed.
#' @param spectral.channel Vector of spectral channel names.
#' @param peak.channel Vector of peak detection channels for fluorophores.
#' @param universal.negative Name of the universal negative control.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param scatter.param Vector of scatter parameters.
#' @param negative.n Integer. Number of events to include in the downsampled
#' negative population. Default is `500`.
#' @param positive.n Integer. Number of events to include in the downsampled
#' positive population. Default is `1000`.
#' @param scatter.match Logical, default is `TRUE`. Whether to select negative
#' events based on scatter profiles matching the positive events. Defines a
#' region of FSC and SSC based on the distribution of selected positive events.
#' @param intermediate.figures Logical, if `TRUE` returns additional figures to
#' show the inner workings of the cleaning, including definition of low-AF cell
#' gates on the PCA-unmixed unstained and spectral ribbon plots of the AF
#' exclusion from the unstained.
#' @param main.figures Logical, if `TRUE` creates the main figures to show the
#' impact of intrusive autofluorescent event removal and scatter-matching for
#' the negatives.
#' @param parallel Logical, default is `FALSE`, in which case parallel processing
#' will not be used. Parallel processing will likely be faster when many small
#' files are read in. If the data is larger, parallel processing may not
#' accelerate the process much.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A list containing the expression data with autofluorescent events
#' removed for each sample.

run.af.removal <- function( clean.expr, af.removal.sample, spectral.channel,
                            peak.channel, universal.negative, asp, scatter.param,
                            negative.n = 500, positive.n = 1000,
                            scatter.match = TRUE,
                            intermediate.figures = FALSE, main.figures = TRUE,
                            parallel = FALSE, verbose = TRUE ) {

  # set up parallel processing
  if ( parallel ){
    future::plan( future::multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future.apply::future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  af.remove.expr <- lapply.function( af.removal.sample, function( sample.name ){

    remove.af( clean.expr, sample.name, spectral.channel, peak.channel,
               universal.negative, asp, scatter.param,
               negative.n, positive.n, scatter.match,
               main.figures, intermediate.figures, verbose )

  } )

  names( af.remove.expr ) <- af.removal.sample

  return( af.remove.expr )
}
