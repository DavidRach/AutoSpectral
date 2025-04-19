# resolution loss function


#' @title Calculate Resolution Loss
#'
#' @description This function calculates the resolution loss for each channel
#'     in the unmixed expression data, based on the maximum expression data.
#'
#' @importFrom stats median mad
#'
#' @param unmixed.expr Matrix containing unmixed expression data.
#' @param expr.data.max Maximum expression data value.
#'
#' @return A numeric vector containing the resolution loss for each channel.
#' @export



get.resolution.loss <- function( unmixed.expr, expr.data.max  ) {

  resolution.loss <- apply( unmixed.expr, 2, function( ch ) {
    resolution.error <- abs( median( ch ) ) / expr.data.max +
      3 * mad( ch ) / expr.data.max

    resolution.error * 10000
  } )

  resolution.loss
}
