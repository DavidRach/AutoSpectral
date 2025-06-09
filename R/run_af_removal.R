# run_af_removal.r

#' @title Run Autofluorescence Removal
#'
#' @description
#' This function runs the autofluorescence removal process on a list of samples,
#' using the specified parameters and settings.
#'
#' @param clean.expr List containing cleaned expression data.
#' @param af.removal.sample Vector of sample names for which autofluorescence
#' removal is to be performed.
#' @param neg.artefacts.list List of negative artefacts.
#' @param spectral.channel Vector of spectral channel names.
#' @param universal.negative Name of the universal negative control.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A list containing the expression data with autofluorescent events
#' removed for each sample.
#'
#' @export

run.af.removal <- function( clean.expr, af.removal.sample, neg.artefacts.list,
                            spectral.channel, universal.negative, asp ) {

  af.remove.expr <- lapply( af.removal.sample, function( sample.name ){

    remove.af( clean.expr, sample.name, neg.artefacts.list, spectral.channel,
               universal.negative, asp )

  } )

  names( af.remove.expr ) <- af.removal.sample

  return( af.remove.expr )
}
