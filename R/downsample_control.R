# downsample_control.r


#' Downsample Control Data
#'
#' This function downsamples control data by selecting a specified number of positive and negative events based on peak channel values.
#'
#' @title Downsample Control Data
#' @description Downsamples control data by selecting a specified number of
#'     positive and negative events based on peak channel values.
#'
#' @param clean.expr.data A list containing cleaned expression data for each sample.
#' @param samp The sample identifier.
#' @param peak.channels A vector of peak channels for the samples.
#' @param negative.n The number of negative events to select.
#' @param positive.n The number of positive events to select.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#'     A list of essential parameters.
#' @return A matrix with the selected expression data.
#' @export



downsample.control <- function( clean.expr.data, samp, peak.channels,
                                 negative.n, positive.n, asp ){

  if ( asp$verbose )
    message( paste( "\033[34m", "Downsampling", samp, "\033[0m" ) )

  # should just pass control's data
  pos.control.expr <- clean.expr.data[[ samp ]]
  # ditto
  peak.channel <- peak.channels[ samp ]

  pos.peak.channel <- pos.control.expr[ , peak.channel ]

  pos.selected <- sort( pos.peak.channel, decreasing = TRUE )[ 1:positive.n ]
  neg.selected <- sort( pos.peak.channel, decreasing = FALSE )[ 1:negative.n ]

  # recover full data
  selected.expr <- pos.control.expr[ names( c( pos.selected, neg.selected ) ), ]

  selected.expr

}
