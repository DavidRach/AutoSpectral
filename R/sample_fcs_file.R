# sample_fcs_file.r

sample.fcs.file <- function( file.name, control.dir, downsample.n, asp ) {
  
  ff <- suppressWarnings( read.FCS( file.path( control.dir, file.name ),
                  transformation = FALSE,
                  truncate_max_range = FALSE,
                  emptyValue = FALSE ) )
  
  ff <- exprs( ff )[ , asp$default.scatter.parameter ]
  
  event.n <- nrow( ff )
  
  check.critical( event.n > asp$min.cell.warning.n, 
                  paste( "Fewer than", asp$min.cell.warning.n, "events in", file.name ) )
  
  ifelse( event.n < downsample.n, downsample.n <- event.n, downsample.n )
  
  fcs.idx <- sample( 1:event.n, downsample.n )
  
  ff <- ff[ fcs.idx, ]
  
  return( ff )
}