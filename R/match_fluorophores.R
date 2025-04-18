# match_fluorophores.r

# match fluorophores in single stained fcs filenames to fluorophore database
# used to create fcs_control_file.csv

match.fluorophores <- function( control.filenames, fluorophore.database ){
  
  delim.start <- "(?:^|[[:space:]_\\.])"
  delim.end <- "(?=$|[[:space:]_\\.])"
  
  fluorophore.matches <- list()
  
  for( filename in control.filenames ) {
    
    fluorophore <- fluorophore.database$fluorophore[ sapply( fluorophore.database$fluorophore, function( fluor ) {
      mod.fluor <- gsub( " ", "\\\\s*", fluor )
      pattern <- paste0( delim.start, mod.fluor, delim.end )
      grepl( pattern, filename, perl = TRUE )
    }) ]
    
    if( length( fluorophore ) == 0 ){
      fluorophore <- fluorophore.database$fluorophore[ sapply( fluorophore.database$synonym1, function( fluor ) {
        mod.fluor <- gsub( " ", "\\\\s*", fluor )
        pattern <- paste0( delim.start, mod.fluor, delim.end )
        grepl( pattern, filename, perl = TRUE )
      }) ]
      
      if( length( fluorophore ) == 0 ){
        fluorophore <- fluorophore.database$fluorophore[ sapply( fluorophore.database$synonym2, function( fluor ) {
          mod.fluor <- gsub( " ", "\\\\s*", fluor )
          pattern <- paste0( delim.start, mod.fluor, delim.end )
          grepl( pattern, filename, perl = TRUE )
        }) ]
        
        if( length( fluorophore ) == 0 ){
          fluorophore <- fluorophore.database$fluorophore[ sapply( fluorophore.database$synonym3, function( fluor ) {
            mod.fluor <- gsub( " ", "\\\\s*", fluor )
            pattern <- paste0( delim.start, mod.fluor, delim.end )
            grepl( pattern, filename, perl = TRUE )
          }) ]
          
        }
      }
    }
    
    if( length( fluorophore ) == 0 ){ 
      
      fluorophore <- "No match"
      
    }
    
    fluorophore.matches[ filename ] <- fluorophore
    
  }
  
  fluorophore.matches <- unlist(fluorophore.matches)
  
  unlist( fluorophore.matches )
  
}