# run_af_artefact_id.r

run.af.artefact.id <- function( clean.expr, negative.sample, 
                                spectral.channel, asp ){
  
  af.artefacts <- lapply( negative.sample, function( sample.name ){
    
    identify.af.artefacts( clean.expr, sample.name,
                           spectral.channel, asp )
    
  } )
  
  names( af.artefacts ) <- negative.sample
  
  return( af.artefacts )
  
}