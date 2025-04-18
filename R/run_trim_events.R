# run_trim_events.r

run.trim.events <- function( trim.sample.data, trim.sample, 
                            trim.peak.channels, trim.factor, asp ){
  
  trimmed.expr <- lapply( names( trim.sample.data ), function( sample.name ){
      trim.extreme.events( trim.sample.data[[ sample.name ]], 
                           trim.peak.channels[[ sample.name ]], trim.factor )
    } )

  names( trimmed.expr ) <- names( trim.sample.data )
  
  
  # issue warning if fewer than 500 events for any sample
  trimmed.expr.n <- sapply( trimmed.expr, nrow )
  
  low.sample.n <- which( trimmed.expr.n < asp$min.cell.warning.n )
  
  if( any( trimmed.expr.n < asp$min.cell.warning.n ) ){
    cat( paste( "Warning! Fewer than", asp$min.cell.warning.n,
                "gated events in", names( low.sample.n ), "\n"  ) )
  }
  
  low.sample.n <- which( trimmed.expr.n < asp$min.cell.stop.n )
  
  # for any samples that have been trimmed to fewer than 50 events,
  # return untrimmed data
  if( any( trimmed.expr.n < asp$min.cell.stop.n ) ){
    for ( low.n in names( low.sample.n ) ) {
      trimmed.expr[[ low.n ]] <- trim.sample.data[[ low.n ]]
    }
  }
  
  rm( trim.sample.data )

  if( any( trimmed.expr.n < asp$min.cell.stop.n ) ){
    cat( paste( "Warning! Fewer than", asp$min.cell.stop.n,
                "gated events in", names( low.sample.n ), "\n"  ) )
  }
  
  trimmed.expr
}






