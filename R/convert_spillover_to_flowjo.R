# convert_spillover_to_flowjo.r

#' @title Convert Spillover To FlowJo
#'
#' @description
#' A quick and dirty method for writing a FlowJo-compatible compensation matrix
#' file from any spillover csv file.
#'
#' @param spillover.file.name The file name of the spillover csv file to be imported.
#' @param flowjo.file.name The name given to the generated FlowJo compensation file.
#' @param input.dir Path to the spillover csv file.
#' @param output.dir Path to the folder where the compensation mtx file will be written.
#'
#' @return No return. Write the xml .mtx file.
#'
#' @export


convert.spillover.to.flowjo <- function( spillover.file.name, flowjo.file.name,
                                         input.dir, output.dir ) {

  # clean up any bad names
  sanitize.name <- function( name ) {
    gsub( "/", "_", name )
  }

  spill <- read.csv( file.path( input.dir, spillover.file.name ),
                     check.names = FALSE, stringsAsFactors = FALSE,
                     row.names = 1 )

  rownames( spill ) <- sanitize.name( rownames( spill ) )
  colnames( spill ) <- sanitize.name( colnames( spill ) )

  # build flowjo-style xml doc
  xml.string <- '<?xml version="1.0" encoding="UTF-8"?>\n<gating:gatingML>\n'
  xml.string <- paste0( xml.string, '  <transforms:spilloverMatrix prefix="Comp-" name="AutoSpectral-spillover" editable="1" status="FINALIZED" transforms:id="">\n' )
  xml.string <- paste0( xml.string, '    <data-type:parameters>\n' )

  for ( fluor in colnames( spill ) ) {
    xml.string <- paste0( xml.string, '      <data-type:parameter data-type:name="', fluor, '"/>\n' )
  }

  xml.string <- paste0(xml.string, '    </data-type:parameters>\n' )

  for ( i in 1:nrow( spill ) ) {
    xml.string <- paste0( xml.string, '    <transforms:spillover data-type:parameter="', rownames( spill )[ i ], '">\n' )
    for ( j in 1:ncol( spill ) ) {
      xml.string <- paste0( xml.string, '      <transforms:coefficient data-type:parameter="', colnames( spill )[ j ], '" transforms:value="', spill[ i, j ], '"/>\n' )
    }
    xml.string <- paste0( xml.string, '    </transforms:spillover>\n' )
  }

  xml.string <- paste0( xml.string, '  </transforms:spilloverMatrix>\n</gating:gatingML>' )

  # save
  writeLines( xml.string, file.path( output.dir, flowjo.file.name ) )
}
