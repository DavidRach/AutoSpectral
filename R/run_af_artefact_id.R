# run_af_artefact_id.r

#' @title Run Autofluorescence Artefact Identification
#'
#' @description
#' This function runs the identification of autofluorescence artefacts on a
#' list of negative samples, using the specified parameters and settings.
#'
#' @param clean.expr List containing cleaned expression data.
#' @param negative.sample Vector of negative sample names for which
#' autofluorescence artefact identification is to be performed.
#' @param spectral.channel Vector of spectral channel names.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A list containing the identified autofluorescence artefacts for
#' each negative sample.
#'
#' @export

run.af.artefact.id <- function( clean.expr, negative.sample,
                                spectral.channel, asp ){

  af.artefacts <- lapply( negative.sample, function( sample.name ){

    id.af.artefacts( clean.expr, sample.name, spectral.channel, asp )

  } )

  names( af.artefacts ) <- negative.sample

  return( af.artefacts )

}
