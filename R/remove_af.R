# remove_af.r

remove.af <- function( clean.expr.data, samp, af.artefact, spectral.channel,
                       universal.negative, asp ) {
  
  if( asp$verbose )
    print( paste( "Removing autofluorescence contamination in", samp ) )
  
  # match universal negative
  matching.negative <- universal.negative[[ samp ]]
  
  # extract matching af.artefact
  matching.artefact <- af.artefact[[ matching.negative ]]
  
  af.components <- matching.artefact[[ 1 ]]
  af.boundaries <- matching.artefact[[ 2 ]]
  
  # get expr data for sample
  expr.data <- clean.expr.data[[ samp ]][ , spectral.channel ]
  
  # unmix with the autofluorescence signatures only
  gate.data <- unmix.ols( expr.data, af.components )
  
  # apply af boundary gates
  if( length( af.boundaries ) == 2 ){
    
    gate.population.pip.lower <- point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      af.boundaries$lower$x, af.boundaries$lower$y )
    
    gate.population.pip.upper <- point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      af.boundaries$upper$x, af.boundaries$upper$y )
    
    gate.population.idx <- which( gate.population.pip.lower == 0 &
                                    gate.population.pip.upper == 0 )
    
  } else {
    
    gate.population.pip <- point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      af.boundaries$upper$x, af.boundaries$upper$y )
    
    gate.population.idx <- which( gate.population.pip == 0 )
    
  }
  
  
  if( !is.null( asp$figure.clean.control.dir ) ) { 
    
    if( length( af.boundaries ) == 2 ){
      
      plot.gate.af.sample( samp, af.data = gate.data,
                           af.boundaries$lower, af.boundaries$upper, 
                           asp )
      
    } else {
      
      plot.gate.af.sample( samp, af.data = gate.data,
                           af.boundary.lower = NULL, 
                           af.boundary.upper = af.boundaries$upper,
                           asp )
    }
  
  }
  
  if( !is.null( asp$figure.spectral.ribbon.dir ) ) {
    
    spectral.ribbon.plot( expr.data, expr.data[ gate.population.idx, ],
                          spectral.channel, asp, samp, 
                          plot.prefix = asp$af.plot.filename,
                          af = TRUE,
                          removed.data = expr.data[ -gate.population.idx, ] )
  }
  
  return( clean.expr.data[[ samp ]][ gate.population.idx, ] )
  
}