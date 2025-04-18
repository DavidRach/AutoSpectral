# resolution loss function

get.resolution.loss <- function( unmixed.expr, expr.data.max  ){
  
  resolution.loss <- apply( unmixed.expr, 2, function( ch ){
    resolution.error <- abs( median( ch ) ) / expr.data.max +
      3 * mad( ch ) / expr.data.max
    
    resolution.error * 10000
  })
 
  resolution.loss
}
  
  