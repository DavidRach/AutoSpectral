# unmix_folder.r

# unmix all fcs files in a directory

unmix.folder <- function( fcs.dir, raw.channels, spectra, asp, 
                          method = "ols",
                          output.dir = NULL,
                          include.raw = FALSE,
                          include.imaging = FALSE,
                          allow.negative = TRUE ){
  
  if( is.null( output.dir ) ){
    output.dir <- asp$unmixed.fcs.dir
  }
  
  files.to.unmix <- list.files( fcs.dir, pattern = ".fcs", full.names = TRUE )
  
  # set up parallel processing
  if( asp$parallel ){
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }
  
  lapply.function( files.to.unmix, FUN = unmix.fcs, raw.channels, spectra, asp,
                   method, output.dir, include.raw, include.imaging, allow.negative )

}