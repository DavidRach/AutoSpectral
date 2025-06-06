# resolution loss function


#' @title Calculate Resolution Loss
#'
#' @description This function calculates the resolution loss for each channel
#' in the unmixed expression data, based on the standard deviation above zero.
#'
#' @importFrom stats sd
#'
#' @param unmixed.expr Matrix containing unmixed expression data.
#'
#' @return A numeric vector containing the resolution loss for each channel.
#'
#' @export



get.resolution.loss <- function( unmixed.expr ) {

  resolution.loss <- apply( unmixed.expr, 2, function( ch ){
    resolution.error <- sd( ch[ ch > 0 ] )
  } )

  resolution.loss
}
