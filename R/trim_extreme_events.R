# trim_extreme_events.r

trim.extreme.events <- function( expr.data, peak.channel, trim.factor ){
  
  # get expr.data for peak channel
  peak.channel.expr <- expr.data[ , peak.channel ]
  
  peak.channel.expr.n <- length( peak.channel.expr )
  
  expr.trim.n <- round( peak.channel.expr.n * trim.factor ) + 1
  
  peak.channel.expr.low <- sort( peak.channel.expr )[ expr.trim.n ]
  peak.channel.expr.high <- sort( peak.channel.expr,
                                  decreasing = TRUE )[ expr.trim.n ]
  
  expr.trim.idx <- which (
    peak.channel.expr > peak.channel.expr.low &
      peak.channel.expr < peak.channel.expr.high )
  
  expr.data[ expr.trim.idx, ]
  
}