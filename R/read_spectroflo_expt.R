# read_spectroflo_expt.r

#' @title read.spectroflo.expt
#' 
#' @description
#' Reads an Experiment (.Expt) file from SpectroFlo and extracts the spillover
#' matrix that was used for the unmixing.
#' 
#' @importFrom xml2 read_xml xml_find_all xml_text
#' @importFrom utils write.csv
#' 
#' @param expt.file File name and path to the .Expt file to be read.
#' @param output.dir Directory where the spillover .csv file will be written.
#' @param output.filename Name for the output spillover file. Default is 
#' `"SpectroFlo_spillover_matrix.csv"`.#' 
#' 
#' @return Spillover matrix (fluorophores x detectors). Fluorophore names will
#' be automatically extracted from the control .fcs files linked to the .Expt 
#' file, with an attempt to clean them up. Detector names (columns) are not 
#' extracted.
#' 
#' @export

read.spectroflo.expt <- function( expt.file, output.dir, 
                                  output.filename = "SpectroFlo_spillover_matrix.csv" ) {
  
  # read SpectroFlo .Expt file (XML format)
  expt <- xml2::read_xml( expt.file )
  
  # find the spillover value nodes
  vector.nodes <- xml_find_all( expt, ".//*[local-name() = '_SpilloverVectorArea']" )
  
  spillover.list <- lapply( vector.nodes, function( node ) {
    float.nodes <- xml_find_all( node, ".//*[local-name()='float']" )
    as.numeric( xml_text( float.nodes ) )
  } )
  
  # convert to a matrix
  spillover.matrix <- do.call( rbind, spillover.list )
 
  # spillover.matrix <- t( apply( spillover.matrix, 1, function( x ) x / max( x ) ) )
  
  # find the associated fluorophore names
  url.nodes <- xml_find_all( expt, ".//*[local-name() = '_Url']" )
  url.paths <- xml_text( url.nodes )
  
  # extract fluorophore names from file names (before "_Controls.fcs")
  fluor.names <- gsub( ".*\\\\|_Controls\\.fcs.*", "", url.paths )
  
  # remove leading plate/rack code (e.g., "B3 ") if present
  fluor.names <- gsub( "^[A-Z]\\d+\\s+", "", fluor.names )
  
  # remove trailing " (Cells)" or similar
  fluor.names <- gsub(" \\(.*\\)", "", fluor.names )
  
  rownames( spillover.matrix ) <- fluor.names
  
  write.csv( spillover.matrix,
             file.path( output.dir, output.filename ) )
  
  return( spillover.matrix )
}
