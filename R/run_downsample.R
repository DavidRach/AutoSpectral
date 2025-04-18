# run_downsample.r

run.downsample <- function( clean.expr.data, downsample.sample, peak.channels,
                            negative.n, positive.n, asp ){
  
  downsample.expr <- lapply( downsample.sample, function( sample.name ){
    
    downsample.control( clean.expr.data, sample.name, peak.channels,
                            negative.n, positive.n, asp )
    
  } )
  
  names( downsample.expr ) <- downsample.sample
  
  return( downsample.expr )
  
}