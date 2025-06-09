# match_fluorophores.r

#' @title Match Fluorophores
#'
#' @description
#' This function matches control filenames to fluorophores in the fluorophore
#' database, including synonyms, and returns the matched fluorophores.
#'
#' @param control.filenames Vector of control filenames.
#' @param fluorophore.database Data frame containing fluorophore information,
#' including synonyms.
#'
#' @return A named vector of matched fluorophores for each control filename.
#'
#' @export

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

    if ( length( fluorophore ) == 0 ){
      fluorophore <- fluorophore.database$fluorophore[ sapply( fluorophore.database$synonym1, function( fluor ) {
        mod.fluor <- gsub( " ", "\\\\s*", fluor )
        pattern <- paste0( delim.start, mod.fluor, delim.end )
        grepl( pattern, filename, perl = TRUE )
      }) ]

      if ( length( fluorophore ) == 0 ){
        fluorophore <- fluorophore.database$fluorophore[ sapply( fluorophore.database$synonym2, function( fluor ) {
          mod.fluor <- gsub( " ", "\\\\s*", fluor )
          pattern <- paste0( delim.start, mod.fluor, delim.end )
          grepl( pattern, filename, perl = TRUE )
        }) ]

        if ( length( fluorophore ) == 0 ){
          fluorophore <- fluorophore.database$fluorophore[ sapply( fluorophore.database$synonym3, function( fluor ) {
            mod.fluor <- gsub( " ", "\\\\s*", fluor )
            pattern <- paste0( delim.start, mod.fluor, delim.end )
            grepl( pattern, filename, perl = TRUE )
          }) ]

        }
      }
    }

    if ( length( fluorophore ) == 0 ){

      fluorophore <- "No match"

    }

    fluorophore.matches[ filename ] <- fluorophore

  }

  fluorophore.matches <- unlist(fluorophore.matches)

  return( unlist( fluorophore.matches ) )
}
