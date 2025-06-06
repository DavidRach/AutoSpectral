# write_flowjo_mtx.r

#' @title Write FlowJo Matrix
#'
#' @description
#' A quick and dirty method for writing a FlowJo-compatible compensation matrix
#' file from the spillover matrix.
#'
#' @param spillover.matrix The spillover matrix, fluorophores x fluorophores.
#' @param flowjo.file.name The path and filename given to the generated FlowJo
#' compensation file.
#'
#' @return No return. Writes the XML .mtx file.
#'
#' @export


write.flowjo.mtx <- function( spillover.matrix, flowjo.file.name ) {

  # clean up any bad names
  sanitize.name <- function( name ) {
    gsub( "/", "_", name )
  }

  rownames( spillover.matrix ) <- sanitize.name( rownames( spillover.matrix ) )
  colnames( spillover.matrix ) <- sanitize.name( colnames( spillover.matrix ) )

  # build flowjo-style xml doc
  xml.string <- '<?xml version="1.0" encoding="UTF-8"?>\n<gating:gatingML>\n'
  xml.string <- paste0( xml.string, '  <transforms:spilloverMatrix prefix="Comp-" name="AutoSpectral-spillover" editable="1" status="FINALIZED" transforms:id="">\n' )
  xml.string <- paste0( xml.string, '    <data-type:parameters>\n' )

  for ( fluor in colnames( spillover.matrix ) ) {
    xml.string <- paste0( xml.string, '      <data-type:parameter data-type:name="', fluor, '"/>\n' )
  }

  xml.string <- paste0(xml.string, '    </data-type:parameters>\n' )

  for ( i in 1:nrow( spillover.matrix ) ) {
    xml.string <- paste0( xml.string, '    <transforms:spillover data-type:parameter="', rownames( spillover.matrix )[ i ], '">\n' )
    for ( j in 1:ncol( spillover.matrix ) ) {
      xml.string <- paste0( xml.string, '      <transforms:coefficient data-type:parameter="', colnames( spillover.matrix )[ j ], '" transforms:value="', spillover.matrix[ i, j ], '"/>\n' )
    }
    xml.string <- paste0( xml.string, '    </transforms:spillover>\n' )
  }

  xml.string <- paste0( xml.string, '  </transforms:spilloverMatrix>\n</gating:gatingML>' )

  # save
  writeLines( xml.string, flowjo.file.name )

}
