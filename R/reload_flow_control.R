# reload_flow_control.r

#' @title Reload Flow Control Information
#'
#' @description
#' This function reloads essential information from control files to permit
#' rapid unmixing at a later date, without recalculating spectra or gates.
#'
#' @importFrom utils read.csv
#' @importFrom dplyr filter
#'
#' @param control.dir file path to the single stained control fcs files
#' @param control.def.file csv file defining the single color control file
#' names, fluorophores they represent, marker names, peak channels and
#' gating requirements.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#'
#' @return A list containing the reloaded flow control information.
#'
#' @export

reload.flow.control <- function( control.dir, control.def.file, asp ) {

  # read channels from controls
  flow.set.channel.table <- read.channel( control.dir, control.def.file, asp )

  flow.set.channel <- flow.set.channel.table[[ 1 ]]
  flow.set.channel.corrected <- flow.set.channel.table[[ 2 ]]

  # read definition of controls
  non.spectral.channel <- asp$non.spectral.channel
  non.spectral.channel <- paste0( non.spectral.channel, collapse = "|" )

  flow.spectral.channel <- flow.set.channel[ !grepl( non.spectral.channel,
                                                     flow.set.channel ) ]

  flow.spectral.channel.n <- length( flow.spectral.channel )

  control.table <- read.csv( control.def.file, na.strings = "",
                             stringsAsFactors = FALSE )

  control.table <- dplyr::filter( control.table, filename != "" )

  flow.fluorophore <- control.table$fluorophore
  flow.fluorophore[ is.na( flow.fluorophore ) ] <- "Negative"

  flow.control.type <- control.table$control.type
  names( flow.control.type ) <- flow.fluorophore

  flow.antigen <- control.table$marker
  flow.channel <- control.table$channel

  check.critical( anyDuplicated( control.table$filename ) == 0,
                  "duplicated filenames in fcs data" )

  # read scatter parameters
  flow.scatter.parameter <- read.scatter.parameter( asp )

  # set labels for time, scatter parameters and channels
  flow.scatter.and.channel <- c( asp$default.time.parameter,
                                 flow.scatter.parameter, flow.channel )
  flow.scatter.and.channel.spectral <- c( asp$default.time.parameter,
                                          flow.scatter.parameter,
                                          flow.spectral.channel )

  flow.scatter.and.channel.label <- c( "Time", flow.scatter.parameter,
                                       ifelse( ! is.na( flow.antigen ),
                                               paste0( flow.antigen, " - ", flow.fluorophore ),
                                               flow.channel ) )
  names( flow.scatter.and.channel.label ) <- flow.scatter.and.channel

  # make control info
  flow.control <- list(
    filename = control.table$filename,
    fluorophore = flow.fluorophore,
    control.type = flow.control.type,
    antigen = flow.antigen,
    channel = flow.channel,
    expr.data.max = asp$expr.data.max,
    expr.data.min = asp$expr.data.min,
    spectral.channel = flow.spectral.channel,
    sample = control.table$sample,
    scatter.and.channel.label = flow.scatter.and.channel.label,
    scatter.and.channel.spectral = flow.scatter.and.channel.spectral,
    scatter.parameter = flow.scatter.parameter
  )

  return( flow.control )

}
