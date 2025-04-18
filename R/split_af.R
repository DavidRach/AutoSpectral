# split AF into multi-AF based on PCA spectra

split.af <- function( af.spectra, flow.control, asp ){
  
  if( ! is.null( af.spectra ) ){
    if( nrow( af.spectra ) > 1 ){
      # split
      # transform/unmix AF sample using af.spectra
      # may need to reconstitute af.spectra so not PCA orthogonal?
      # pick top 200 events
      # add clean AF-removed unstained
      # set as event.sample factor
    }
      
  }
  
}
