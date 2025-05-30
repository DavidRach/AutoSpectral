# define_flow_control.r

#' @title Define Flow Control
#' @description
#' A complex function designed to convert the single-stained control fcs files
#'     and the metadata control.def.file into a data structure for AutoSpectral.
#'     Reads the metadata file and locates the corresponding fcs files.
#'     Determines the gates that will be needed based on varibles such as beads/cells,
#'     large.gate and viability.gate, and then creates gates for each combination.
#'     Imports and gates the data in the fcs files.
#'     Assigns various factors to track the data.
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @importFrom dplyr filter
#'
#' @param control.dir file path to the single stained control fcs files
#' @param control.def.file csv file defining the single color control file names,
#'     fluorophores they represent, marker names, peak channels and gating requirements.
#' @param asp The AutoSpectral parameter list defined using get.autospectral.param.
#' @param gate Logical, default is TRUE. Not currently in use. To be implemented
#'     to allow pre-gated data to be read in without applying additional gating here.
#'
#' @return flow.control
#'     Returns the flow.control list with the following components:
#' \describe{
#'   \item{filename}{Names of the single color control files.}
#'   \item{fluorophore}{Corresponding fluorophores used in the experiment.}
#'   \item{antigen}{Corresponding markers used in the experiment.}
#'   \item{control.type}{Type of control used (beads or cells).}
#'   \item{universal.negative}{Corresponding universal negative for each control.}
#'   \item{viability}{Logical factor; whether a control represents a viability marker.}
#'   \item{large.gate}{Logical factor; large gate setting.}
#'   \item{autof.channel.idx}{Index of the autofluorescence marker channel.}
#'   \item{event.number.width}{Width of the event number.}
#'   \item{expr.data.max}{Maximum expression data value.}
#'   \item{expr.data.max.ceil}{Ceiling of the maximum expression data value.}
#'   \item{expr.data.min}{Minimum expression data value.}
#'   \item{channel}{Preliminary peak channels for the fluorophores.}
#'   \item{channel.n}{Number of channels.}
#'   \item{spectral.channel}{Spectral channel information.}
#'   \item{spectral.channel.n}{Number of spectral channels.}
#'   \item{sample}{Sample names (fluorophores).}
#'   \item{scatter.and.channel}{FSC, SSC and peak channel information.}
#'   \item{scatter.and.channel.label}{Labels for scatter and channel.}
#'   \item{scatter.and.channel.spectral}{FSC, SSC and spectral channels.}
#'   \item{scatter.parameter}{Scatter parameters used for gating.}
#'   \item{event}{Event factor.}
#'   \item{event.n}{Number of events.}
#'   \item{event.sample}{Sample information for events. Links events to samples.}
#'   \item{event.type}{Type of events.}
#'   \item{expr.data}{Expression data. The data to be used for extracting spectra.}
#' }
#'
#' @export


define.flow.control <- function( control.dir, control.def.file, asp,
                                    gate = TRUE )
{
  # set up parallel processing
  if ( asp$parallel ){
      future::plan( future::multisession, workers = asp$worker.process.n )
      options( future.globals.maxSize = asp$max.memory.n )
      lapply.function <- future.apply::future_lapply
    } else {
      lapply.function <- lapply.sequential
    }


    # read channels from controls
    if (asp$verbose )
      message( "\033[34m Reading control information \033[0m" )

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

    check.critical( anyDuplicated( control.table$filename ) == 0,
        "duplicated filenames in fcs data" )

    # define samples (replicate variously gated negatives)
    if (asp$verbose )
      message( "\033[34m Matching negatives for the controls \033[0m" )

    control.table$sample <- control.table$fluorophore

    control.table$universal.negative[ is.na( control.table$universal.negative )] <- FALSE
    control.table$is.viability[ is.na( control.table$is.viability )] <- FALSE
    control.table$large.gate[ is.na( control.table$large.gate )] <- FALSE

    negative.types <- data.frame( negative = control.table$universal.negative,
                                  large.gate = control.table$large.gate,
                                  viability = control.table$is.viability )

    unique.neg <- unique( negative.types )
    #unique.neg <- unique.neg[ unique.neg$negative == TRUE, , drop = FALSE ]

    if ( nrow( unique.neg != 0 ) ) {
      # match samples with corresponding negative by matching gates
      for ( i in 1 : nrow( unique.neg ) ) {
        match.found <- any( control.table$filename == unique.neg$negative[ i ] &
                              control.table$large.gate == unique.neg$large.gate[ i ] &
                              control.table$is.viability == unique.neg$viability[ i ] )

        if ( !match.found ) {
          matching.row <- control.table[ control.table$filename == unique.neg$negative[ i ], ]
          new.row <- matching.row[ 1, ]
          new.row$large.gate <- unique.neg$large.gate[ i ]
          new.row$is.viability <- unique.neg$viability[ i ]
          new.row$fluorophore <- "Negative"
          new.row$universal.negative <- FALSE
          new.row$sample <- paste( new.row$control.type, "Negative", i )
          control.table <- rbind( control.table, new.row )
        }
      }
    }

    # define gates needed
    if ( asp$verbose )
      message( "\033[34m Determining the gates that will be needed \033[0m" )

    gate.types <- data.frame( type = control.table$control.type,
                              viability = control.table$is.viability,
                              large.gate = control.table$large.gate )
    unique.gates <- unique( gate.types )
    sample.matches <- apply( gate.types, 1, function( row ) {
      match( TRUE, apply( unique.gates, 1,
                          function( unique.row ) all( row == unique.row ) ) )
    } )

    control.table$gate <- sample.matches
    flow.gate <- sample.matches
    gate.type <- unique( sample.matches )

    # matching universal negatives by desired sample and gating
    # with an exception for viability gating mismatch
    control.table$universal.negative <- sapply( 1:nrow( control.table ), function( i ) {
      if ( control.table$universal.negative[ i ] != "FALSE" ) {
        match.row <- control.table[ control.table$filename ==
                                      control.table$universal.negative[ i ] &
                                      control.table$gate == control.table$gate[ i ], ]
        if ( nrow( match.row ) > 0 ) {
          return( match.row$sample )
        } else {
          match.row <- control.table[ control.table$filename ==
                                        control.table$universal.negative[ i ], ]
          if ( nrow( match.row ) > 0 ) {
            return( match.row$sample )
          }
        }
      }
      return( control.table$universal.negative[ i ] )
    })

    # set samples and gate combos
    flow.sample <- control.table$sample
    flow.sample.n <- length( flow.sample )
    names( flow.gate ) <- flow.sample

    flow.file.name <- control.table$filename
    names( flow.file.name ) <- flow.sample

    flow.fluorophore <- control.table$fluorophore
    flow.fluorophore[ is.na( flow.fluorophore ) ] <- "Negative"

    flow.antigen <- control.table$marker
    flow.channel <- control.table$channel
    flow.control.type <- control.table$control.type

    flow.universal.negative <- control.table$universal.negative
    flow.viability <- control.table$is.viability
    flow.large.gate <- control.table$large.gate

    flow.antigen[ flow.fluorophore == "AF" ] <- "AF"
    flow.antigen[ is.na( flow.antigen ) ] <- "other"
    if ( is.na( flow.channel[ flow.fluorophore == "AF" ] ) )
        flow.channel[ flow.fluorophore == "AF" ] <- asp$af.channel

    flow.channel[ is.na( flow.channel ) ] <- "other"

    flow.universal.negative[ is.na( flow.universal.negative ) ] <- FALSE
    names( flow.universal.negative ) <- flow.sample
    flow.viability[ is.na( flow.viability ) ] <- FALSE
    flow.large.gate[ is.na( flow.large.gate ) ] <- FALSE

    flow.channel.n <- length( flow.channel )

    flow.autof.marker.idx <- which( flow.antigen == asp$antigen.autof )
    if ( length( flow.autof.marker.idx ) != 1 )
        flow.autof.marker.idx <- NULL

    names( flow.viability ) <- flow.sample
    names( flow.large.gate ) <- flow.sample

    # read scatter parameters
    if (asp$verbose )
      message( "\033[34m Determining channels to be used \033[0m" )

    flow.scatter.parameter <- read.scatter.parameter( asp )

    # set scatter parameters and channels
    flow.scatter.and.channel <- c( asp$default.time.parameter,
                                   flow.scatter.parameter, flow.channel )

    flow.scatter.and.channel.spectral <- c( asp$default.time.parameter,
                                            flow.scatter.parameter,
                                            flow.spectral.channel )

    flow.scatter.and.channel.matched.bool <-
      flow.scatter.and.channel.spectral %in% flow.set.channel

    if ( ! all( flow.scatter.and.channel.matched.bool ) )
    {
        channel.matched <-
          flow.scatter.and.channel.spectral[ flow.scatter.and.channel.matched.bool ]
        flow.scatter.and.channel.unmatched <- paste0(
            sort( setdiff( flow.scatter.and.channel.spectral, channel.matched ) ),
            collapse = ", "
        )
        flow.set.unmatched <- paste0(
            sort( setdiff( flow.set.channel, channel.matched ) ),
            collapse = ", "
        )
        error.msg <- sprintf(
            "wrong channel name, not found in fcs data\n\texpected: %s\n\tfound: %s",
            flow.scatter.and.channel.unmatched, flow.set.unmatched
        )
        check.critical( FALSE, error.msg )
    }

    names( flow.channel ) <- flow.sample

    check.critical( ! anyDuplicated( flow.scatter.and.channel.spectral ),
        "internal error: names for channels overlap" )

    # set labels for time, scatter parameters and channels
    flow.scatter.and.channel.label <- c( "Time", flow.scatter.parameter,
        ifelse( ! is.na( flow.antigen ),
            paste0( flow.antigen, " - ", flow.fluorophore ),
            flow.channel ) )
    names( flow.scatter.and.channel.label ) <- flow.scatter.and.channel

    # get range of fcs data
    flow.set.resolution <- asp$expr.data.max

    flow.expr.data.min <- asp$expr.data.min
    flow.expr.data.max <- asp$expr.data.max
    flow.expr.data.max.ceil <- ceiling( flow.expr.data.max / asp$data.step ) *
      asp$data.step

    # create figure and table directories
    if (asp$verbose )
      message( "\033[34m Creating output folders \033[0m" )

    create.directory( asp )

    # define gates on downsampled pooled fcs by type
    if (asp$verbose )
      message( "\033[34m Defining the gates \033[0m" )

    gate.list <- list()

    for( gate in gate.type ){

      files.to.gate <- unique( control.table[ control.table$gate == gate, ]$filename )

      files.n <- length( files.to.gate )

      downsample.n <- asp$gate.downsample.n.cells / files.n

      scatter.coords <- lapply.function( files.to.gate, sample.fcs.file,
                                         control.dir, downsample.n, asp,
                                         future.seed = asp$gate.downsample.seed )

      scatter.coords <- do.call( rbind, scatter.coords )

      viability.gate <- unique( control.table[ control.table$gate == gate, ]$is.viability )

      is.viability <- ifelse( viability.gate, "viabilityMarker", "nonViability" )

      large.gate <- unique( control.table[ control.table$gate == gate, ]$large.gate )

      is.large.gate <- ifelse( large.gate, "largeGate", "smallGate" )

      control.type <- unique( control.table[ control.table$gate == gate, ]$control.type )

      samp <- paste( control.type,is.viability, is.large.gate, sep = "_" )

      gate.boundary <- do.gate( scatter.coords, viability.gate, large.gate,
                                samp, flow.scatter.and.channel.label, control.type, asp )

      gate.list[[ gate ]] <- gate.boundary

    }

    # read in fcs files, selecting data within pre-defined gates
    if (asp$verbose )
      message( "\033[34m Reading FCS files \033[0m" )

    flow.expr.data <- lapply.function( flow.sample, FUN = get.gated.flow.expression.data,
                              flow.file.name, control.dir,
                              scatter.and.spectral.channel = flow.scatter.and.channel.spectral,
                              spectral.channel = flow.spectral.channel,
                              set.resolution = flow.set.resolution,
                              flow.gate = flow.gate, gate.list = gate.list,
                              scatter.param = flow.scatter.parameter,
                              scatter.and.channel.label = flow.scatter.and.channel.label,
                              asp = asp )

    names( flow.expr.data ) <- flow.sample

    # organize data
    if (asp$verbose )
      message( "\033[34m Organizing control info \033[0m" )

    flow.sample.event.number.max <- 0

    for ( fs.idx in 1 : flow.sample.n )
    {
        flow.sample.event.number <- nrow( flow.expr.data[[ fs.idx ]]  )

        rownames( flow.expr.data[[ fs.idx ]] ) <- paste( flow.sample[ fs.idx ],
                                seq_len( flow.sample.event.number ), sep = "_")

        if ( flow.sample.event.number < 500 )
          warning( paste( "\033[31m",  "Warning! Fewer than 500 gated events in",
                      flow.file.name[ fs.idx ],
                      "\033[0m", "\n" ) )

        if ( flow.sample.event.number > flow.sample.event.number.max )
            flow.sample.event.number.max <- flow.sample.event.number
    }

    flow.event.number.width <-
        floor( log10( flow.sample.event.number.max ) ) + 1
    flow.event.regexp <- sprintf( "\\.[0-9]{%d}$", flow.event.number.width )

    # set rownames
    for ( fs.idx in 1 : flow.sample.n )
    {
      flow.sample.event.number <- nrow( flow.expr.data[[ fs.idx ]]  )
      flow.the.sample <- flow.sample[ fs.idx ]

      flow.the.event <- sprintf( "%s.%0*d", flow.the.sample,
                                 flow.event.number.width, 1 : flow.sample.event.number )
      rownames( flow.expr.data[[ fs.idx ]] ) <- flow.the.event
    }

    flow.expr.data <- do.call( rbind, flow.expr.data )

    # set events
    flow.event <- rownames( flow.expr.data )

    flow.event.n <- length( flow.event )

    flow.event.sample <- sub( flow.event.regexp, "", flow.event )

    flow.event.sample <- factor( flow.event.sample, levels = flow.sample )

    event.type.factor <- flow.sample
    names( event.type.factor ) <- flow.control.type
    flow.event.type <- factor( flow.event.sample,
                               levels = event.type.factor,
                               labels = names( event.type.factor ) )

    names( flow.control.type ) <- flow.fluorophore

    # make control info
    flow.control <- list(
        filename = flow.file.name,
        fluorophore = flow.fluorophore,
        antigen = flow.antigen,
        control.type = flow.control.type,
        universal.negative = flow.universal.negative,
        viability = flow.viability,
        large.gate = flow.large.gate,
        autof.channel.idx = flow.autof.marker.idx,
        event.number.width = flow.event.number.width,
        expr.data.max = flow.expr.data.max,
        expr.data.max.ceil = flow.expr.data.max.ceil,
        expr.data.min = flow.expr.data.min,
        channel = flow.channel,
        channel.n = flow.channel.n,
        spectral.channel = flow.spectral.channel,
        spectral.channel.n = flow.spectral.channel.n,
        sample = flow.sample,
        scatter.and.channel = flow.scatter.and.channel,
        scatter.and.channel.label = flow.scatter.and.channel.label,
        scatter.and.channel.spectral = flow.scatter.and.channel.spectral,
        scatter.parameter = flow.scatter.parameter,
        event = flow.event,
        event.n = flow.event.n,
        event.sample = flow.event.sample,
        event.type = flow.event.type,
        expr.data = flow.expr.data
    )

    if (asp$verbose )
      message( "\033[32m Control setup complete! \n Review gates in figure_gate. \033[0m" )

    flow.control
}

