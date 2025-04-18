# run_peacoqc_test.r

run.peacoQC <- function( expr.data, spectral.channel, 
                              all.channels, asp ){
  
  # set up parallel processing
  if( asp$parallel ){
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }
  
  # define parameters for peacoQC
  biexp.transform <- flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = FALSE
  )
  
  transform.inv <- flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue = asp$default.transformation.param$max.range,
    pos = asp$default.transformation.param$pos,
    neg = asp$default.transformation.param$neg,
    widthBasis = asp$default.transformation.param$width,
    inverse = TRUE
  )
  
  time.param <- asp$default.time.parameter

  # run peacoQC to remove flow fluctuation errors
  clean.expr <- lapply.function( names( expr.data ), function( sample.name ) {
    do.peacoQC( expr.data[[ sample.name ]], sample.name, 
                spectral.channel, biexp.transform, transform.inv,
                asp$figure.peacoqc.dir, time.param, all.channels )
    
  }, future.seed = asp$gate.downsample.seed )
  
  # note that PeacoQC will only run on MAD for low n samples like controls
  # do.peacoQC is therefore set to run with MAD only
  names( clean.expr ) <- names( expr.data )
  
  rm( expr.data )
  
  # issue warning if fewer than 500 events for any sample
  clean.expr.n <- sapply( clean.expr, nrow )
  
  low.sample.n <- which( clean.expr.n < asp$min.cell.warning.n )
  
  if( any( clean.expr.n < asp$min.cell.warning.n ) ){
    cat( paste( "Warning! Fewer than", asp$min.cell.warning.n, 
                "gated events in", names( low.sample.n ), "\n"  ) )
  }
  
  clean.expr
  
}