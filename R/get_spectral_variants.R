# get_spectral_variants.r

#' @title Get Spectral Variations for Fluorophores
#'
#' @description
#' This function cycles through all the fluorophores defined in `control.def.file`,
#' identifying variations in spectral profiles. It does this by performing SOM
#' clustering on the positive events in the cleaned control data. The output is
#' saved as an .rds file, and figures summarizing the variation are saved, if
#' desired. Note that the .rds file contains all the needed information for
#' downstream processing (per-cell unmixing), so you can just load that using
#' the `readRDS` function) rather than re-running this process.
#'
#' @importFrom EmbedSOM SOM
#'
#' @param control.dir File path to the single-stained control FCS files.
#' @param control.def.file CSV file defining the single-color control file names,
#' fluorophores they represent, marker names, peak channels, and gating requirements.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param spectra A matrix containing the spectral data. Fluorophores in rows,
#' detectors in columns.
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`.
#' @param n.cells Numeric, default `2000`. Number of cells to use for defining
#' the variation in spectra. Up to `n.cells` cells will be selected as positive
#' events in the peak channel for each fluorophore, above the `pos.quantile` in
#' the unstained sample.
#' @param pos.quantile Numeric, default `0.995`. Threshold for positivity. This
#' quantile will be used to define the maximum extent of the unstained sample in
#' the unmixed space. Anything above that will be considered positive in the
#' single-stained controls.
#' @param som.dim Numeric, default `10`. Number of x and y dimensions to use in
#' the SOM for clustering the spectral variation.
#' @param sim.threshold Numeric, default `0.98`. Threshold for cosine similarity-
#' based exclusion of spectral variants. Any variant less than `sim.threshold`
#' by `cosine.similarity` from the optimized spectrum for that fluorophore (from
#' `spectra`) will be excluded from output. This helps to exclude autofluorescence
#' contamination.
#' @param figures Logical, controls whether the variation in spectra for each
#' fluorophore is plotted in `output.dir`. Default is `TRUE`.
#' @param output.dir File path to whether the figures and .rds data file will be
#' saved. Default is `NULL`, in which case `asp$variant.dir` will be used.
#' @param parallel Logical, default is `FALSE`, in which case sequential processing
#' will be used. Parallel processing will likely be faster when many small
#' files are read in. If the data is larger, parallel processing may not
#' accelerate the process much or may fail outright.
#' @param verbose Logical, default is `TRUE`. Set to `FALSE` to suppress messages.
#'
#' @return A vector with the indexes of events inside the initial gate.
#'
#' @export

get.spectral.variants <- function( control.dir, control.def.file,
                                   asp, spectra, af.spectra,
                                   n.cells = 2000, pos.quantile = 0.995,
                                   som.dim = 10, sim.threshold = 0.98,
                                   figures = TRUE, output.dir = NULL,
                                   parallel = FALSE, verbose = TRUE ) {

  if ( is.null( af.spectra ) )
    stop( "Multiple AF spectra must be provided." )
  if ( nrow( af.spectra ) < 2 )
    stop( "Multiple AF spectra must be provided." )

  if ( is.null( output.dir ) )
    output.dir <- asp$variant.dir
  if ( !dir.exists( output.dir ) )
    dir.create( output.dir )

  fluorophores <- rownames( spectra )[ rownames( spectra ) != "AF" ]
  spectral.channel <- colnames( spectra )

  # read control file
  if ( !file.exists( control.def.file ) )
    stop( paste( "Unable to locate control.def.file:", control.def.file ) )

  if ( verbose ) message( "\033[34m Checking control file for errors \033[0m" )
  check.control.file( control.dir, control.def.file, asp, strict = TRUE )

  control.table <- read.csv( control.def.file, na.strings = "",
                             stringsAsFactors = FALSE )
  control.table <- dplyr::filter( control.table, filename != "" )

  # set channels to be used
  flow.set.resolution <- asp$expr.data.max
  flow.set.channel.table <- read.channel( control.dir, control.def.file, asp )
  flow.set.channel <- flow.set.channel.table[[ 1 ]]
  non.spectral.channel <- asp$non.spectral.channel
  non.spectral.channel <- paste0( non.spectral.channel, collapse = "|" )

  flow.spectral.channel <- flow.set.channel[ !grepl( non.spectral.channel,
                                                     flow.set.channel ) ]
  flow.scatter.parameter <- read.scatter.parameter( asp )
  flow.scatter.and.channel.spectral <- c( asp$default.time.parameter,
                                          flow.scatter.parameter,
                                          flow.spectral.channel )

  if ( grepl( "Discover", asp$cytometer ) )
    flow.spectral.channel <- flow.spectral.channel[ grep( asp$spectral.channel,
                                                          flow.spectral.channel ) ]

  # define list of samples
  flow.sample <- control.table$sample
  table.fluors <- control.table$fluorophore
  table.fluors <- table.fluors[ !is.na( table.fluors ) ]
  #table.fluor.n <- length( table.fluors )
  universal.negative <- control.table$universal.negative
  universal.negative[ is.na( universal.negative ) ] <- FALSE
  names( universal.negative ) <- table.fluors
  flow.channel <- control.table$channel
  names( flow.channel ) <- table.fluors
  flow.file.name <- control.table$filename
  names( flow.file.name ) <- table.fluors
  control.type <- control.table$control.type
  names( control.type ) <- table.fluors

  # stop if "AF" sample is not present, fluorophore mismatch
  if ( !( "AF" %in% table.fluors ) )
    stop( "Unable to locate `AF` control in control file. An unstained cell control is required" )
  #if ( length( rownames( spectra ) ) != table.fluor.n )
  #  stop( "The number of fluorophores in your control file doesn't match that in your spectra." )

  if ( ! all( table.fluors %in% fluorophores ) ) {
    # check for 'Negative', 'AF', check for match again
    fluor.to.match <- table.fluors[ !grepl( "Negative", table.fluors ) ]
    fluor.to.match <- fluor.to.match[ !fluor.to.match == "AF" ]
    matching.fluors <- fluor.to.match %in% fluorophores

    if ( !all( matching.fluors ) ) {
      warning( "The fluorophores in your control file don't match those in your spectra." )
      if ( !any( matching.fluors ) )
        stop()
    }
    table.fluors <- fluor.to.match[ matching.fluors ]
  }

  if ( parallel ) {
    future::plan( future::multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future.apply::future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  # get thresholds for positivity
  if ( verbose ) message( paste( "\033[33m", "Calculating positivity thresholds", "\033[0m" ) )
  unstained <- suppressWarnings(
    flowCore::read.FCS( file.path( control.dir, flow.file.name[ "AF" ] ),
                        transformation = NULL,
                        truncate_max_range = FALSE,
                        emptyValue = FALSE ) )

  # read exprs for spectral channels only
  if ( nrow( unstained ) > asp$gate.downsample.n.cells ) {
    set.seed( asp$gate.downsample.seed )
    unstained.idx <- sample( nrow( unstained ), asp$gate.downsample.n.beads )
    unstained <- flowCore::exprs( unstained )[ unstained.idx, spectral.channel ]
  } else {
    unstained <- flowCore::exprs( unstained )[ , spectral.channel ]
  }

  raw.thresholds <- apply( unstained, 2, function( col ) quantile( col, pos.quantile ) )

  unstained.unmixed <- unmix.autospectral( unstained, spectra, af.spectra, verbose = FALSE )
  unmixed.thresholds <- apply( unstained.unmixed[ , fluorophores ], 2, function( col )
    quantile( col, pos.quantile ) )

  # main loop
  if ( verbose ) message( paste( "\033[33m", "Getting spectral variants", "\033[0m" ) )
  spectral.variants <- lapply.function( table.fluors, get.fluor.variants, flow.file.name,
                               control.dir, spectra, af.spectra, flow.spectral.channel,
                               universal.negative, control.type, raw.thresholds,
                               unmixed.thresholds, flow.channel, som.dim, n.cells,
                               asp, verbose, output.dir, sim.threshold, figures,
                               future.seed = asp$variant.seed )

  names( spectral.variants ) <- fluorophores

  if ( verbose ) message( paste( "\033[33m", "Spectral variation computed!", "\033[0m" ) )

  variants <- list(
    thresholds = unmixed.thresholds,
    variants = spectral.variants
  )

  saveRDS( variants, file = file.path( output.dir, asp$variant.filename ) )

  return( variants )

}
