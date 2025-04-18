# run_af_removal.r

run.af.removal <- function( clean.expr, af.removal.sample, neg.artefacts.list, 
                            spectral.channel, universal.negative, asp ) {
  
  af.remove.expr <- lapply( af.removal.sample, function( sample.name ){
    
    remove.af( clean.expr, sample.name, neg.artefacts.list, spectral.channel,
               universal.negative, asp )
    
  } )
  
  names( af.remove.expr ) <- af.removal.sample
  
  return( af.remove.expr )
}