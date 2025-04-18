# run_universal_negative.r

run.universal.negative <- function( clean.expr, univ.sample,
                                    universal.negatives, 
                                    scatter.param,
                                    peak.channels, downsample,
                                    negative.n, positive.n,
                                    spectral.channel, asp,
                                    control.type,
                                    scatter.match ){
  
  univ.expr <- lapply( univ.sample, function( sample.name ){
    
    get.universal.negative( clean.expr, sample.name, universal.negatives,
                                 scatter.param,
                                 peak.channels, downsample,
                                 negative.n, positive.n,
                                 spectral.channel, asp,
                                 control.type,
                                 scatter.match )
    
  } )
  
  names( univ.expr ) <- univ.sample
  
  rm( clean.expr )
  
  return( univ.expr )
  
}