# clean_controls.r

#' @title Clean Controls
#'
#' @description
#' A four-part function to clean single-color controls in order to extract
#' fluorophore signatures. Any part can be run independently:
#'
#' - **Stage 1**: PeacoQC to eliminate flow artefacts.
#' - **Stage 2**: Trimming to eliminate extreme events.
#' - **Stage 3**: Autofluorescence noise removal using PCA unmixing on matching
#'   unstained (cells only).
#' - **Stage 4**: Brightest event selection from positive, universal negative
#'   from matching negative, and downsampling to speed up RLM spectra
#'   optimization.
#'
#' @param flow.control A list prepared using `define.flow.control`, containing
#' the data and essential information about the cytometer and data structure.
#' @param asp The AutoSpectral parameter list, prepared using
#' `get.autospectral.param`.
#' @param time.clean Logical, default is `FALSE`. Whether to run PeacoQC to
#' remove time-based inconsistencies in the controls.
#' @param trim Logical, default is `FALSE`. Whether to remove extreme events
#' (positive and negative) from controls.
#' @param trim.factor Numeric. Default is `asp$rlm.trim.factor`. Required if
#' `trim = TRUE`.
#' @param af.remove Logical, default is `FALSE`. Whether to remove intrusive
#' autofluorescence contamination from cell controls using PCA-based
#' identification and gating. Requires universal negatives to be defined in the
#' control file and in `flow.control`.
#' @param universal.negative Logical, default is `TRUE`. Whether to use a
#' universal negative sample as the negative for spectral extraction. Requires
#' universal negatives to be defined in the control file and in `flow.control`.
#' @param downsample Logical, default is `TRUE`. Whether to reduce cell and bead
#' control events to speed up processing.
#' @param negative.n Integer. Number of events to include in the downsampled
#' negative population. Default is `asp$negative.n`.
#' @param positive.n Integer. Number of events to include in the downsampled
#' positive population. Default is `asp$positive.n`.
#' @param scatter.match Logical, default is `TRUE`. Whether to select negative
#' events based on scatter profiles matching the positive events. Defines a
#' region of FSC and SSC based on the distribution of selected positive events.
#' @param scrub Logical, if `TRUE` allows for re-cleaning of already cleaned
#' data, provided there are clean data in `flow.control`.
#'
#' @return
#' Returns a modified `flow.control` with the original data intact. New, cleaned
#' data and corresponding factors are stored in new slots.
#'
#' @export

clean.controls <- function( flow.control, asp,
                            time.clean = FALSE,
                            trim = FALSE, trim.factor = NULL,
                            af.remove = FALSE,
                            universal.negative = TRUE, downsample = TRUE,
                            negative.n = asp$negative.n, positive.n = asp$positive.n,
                            scatter.match = TRUE, scrub = TRUE ) {

  flow.sample <- flow.control$sample
  flow.sample.n <- length( flow.sample )
  flow.control.type <- flow.control$control.type
  clean.event.type <- flow.control$event.type

  if ( scrub & !is.null( flow.control$clean.expr ) ) {

    clean.expr <- flow.control$clean.expr
    clean.event.sample <- flow.control$clean.event.sample
    flow.negative <- flow.control$clean.universal.negative
    clean.universal.negative <- flow.control$clean.universal.negative

  } else {

    clean.expr <- flow.control$expr.data
    clean.event.sample <- flow.control$event.sample
    flow.negative <- flow.control$universal.negative
    clean.universal.negative <- NULL
  }

  # split expression data by sample into a list
  clean.expr <- lapply( flow.sample, function( fs ) {
    clean.expr[ clean.event.sample == fs, ]
  } )

  names( clean.expr ) <- flow.sample

  spectral.channel <- flow.control$spectral.channel

  all.channels <- flow.control$scatter.and.channel.original

  ### Stage 1: Use PeacoQC to clean up flow  -----------------
  # clean based on irregularities in flow by TIME parameter
  # this can be slow and may not show much effect

  if ( time.clean ) {
    if ( !dir.exists( asp$figure.peacoqc.dir ) )
      dir.create( asp$figure.peacoqc.dir )

    clean.expr <- run.peacoQC( clean.expr, spectral.channel,
                               all.channels, asp )
  }

  ### Stage 2: Trimming -----------------
  # clean by trimming extreme events
  # to be used in case of aggregates
  # not recommended for regular use due to loss of brightest events

  if ( trim ) {

    if ( is.null( trim.factor ) ) {
      trim.factor <- asp$rlm.trim.factor
    }

    # trim fluorophore controls only
    trim.sample <- flow.control$fluorophore[ ! grepl( "AF|negative",
                                                      flow.control$fluorophore,
                                                      ignore.case = TRUE ) ]

    trim.peak.channels <- flow.control$channel[ flow.control$fluorophore
                                                %in% trim.sample ]

    trim.sample.data <- clean.expr[ trim.sample ]

    trimmed.expr <- run.trim.events( trim.sample.data, trim.sample,
                                     trim.peak.channels, trim.factor, asp )

    rm( trim.sample.data )

    # merge in trimmed data
    clean.expr[ names( trimmed.expr ) ] <- trimmed.expr

  }

  ### Stage 3: Remove Autofluorescence intrusions -----------------
  # clean by removing autofluorescence contamination
  # recommended for controls from tissues (e.g., mouse splenocytes)
  # not needed for PBMCs
  # can be done on a per control basis by specifying clean = TRUE in the fcs_control_file

  if ( af.remove ) {
    if ( !dir.exists( asp$figure.clean.control.dir ) )
      dir.create( asp$figure.clean.control.dir )
    if ( !dir.exists( asp$figure.spectral.ribbon.dir ) )
      dir.create( asp$figure.spectral.ribbon.dir )
    if ( !dir.exists( asp$figure.af.dir ) )
      dir.create( asp$figure.af.dir )

    # identify universal negative cell samples
    univ.neg <- unique( flow.negative )
    univ.neg <- univ.neg[ !is.null( univ.neg ) ]
    univ.neg <- univ.neg[ !is.na( univ.neg ) ]
    univ.neg <- univ.neg[ univ.neg != FALSE ]

    check.critical( length( univ.neg ) > 0,
                    "No cell-based universal negative samples could be identified.
                    To perform autofluorescence removal, you must specify a universal negative in the fcs_control_file." )

    univ.neg.sample <- flow.control$sample[ flow.control.type == "cells" &
                                              flow.sample %in% univ.neg ]

    check.critical( length( univ.neg.sample ) > 0,
                    "No cell-based universal negative samples could be identified.
                    To perform autofluorescence removal, you must specify a
                    universal negative in the fcs_control_file,
                    and you must have single-stained cell controls." )


    # select cell-based single-stained samples to be used
    # must be cells and must have a corresponding universal negative
    af.removal.sample <- flow.sample[ flow.control.type == "cells" &
                                        flow.negative %in% univ.neg ]

    # check that is not length 0 and stop if is
    check.critical( length( af.removal.sample ) > 0,
                    "Error: no cell-based samples with a corresponding universal negative could be identified.
                    To perform autofluorescence removal, you must specify a universal negative in the fcs_control_file." )

    # get af artefacts for each negative control
    univ.neg.af.artefacts <- run.af.artefact.id( clean.expr, univ.neg.sample,
                                                 spectral.channel, asp )

    # remove identified AF from single-color controls
    af.removed.expr <- run.af.removal( clean.expr, af.removal.sample,
                                       univ.neg.af.artefacts,
                                       spectral.channel,
                                       flow.negative,
                                       asp )

    # if AF is among cleaned controls, rename AF in cleaned controls to AF cleaned
    # names( af.removed.expr )[ names( af.removed.expr ) == "AF" ] <- "AF cleaned"

    # replace AF in universal negative with AF cleaned
    clean.universal.negative <- flow.negative
    # clean.universal.negative[ clean.universal.negative == "AF" ] <- "AF cleaned"

    flow.control$clean.universal.negative <- clean.universal.negative

    # store in all except AF and AF cleaned
    # get names in case of only cleaning certain samples
    # cleaned.samples <- names( af.removed.expr )[ names( af.removed.expr ) != "AF cleaned" ]

    cleaned.samples <- names( af.removed.expr )
    clean.expr[ cleaned.samples ] <- af.removed.expr[ cleaned.samples ]

    # store AF cleaned in a new slot
    #clean.expr[ "AF cleaned" ] <- af.removed.expr[ names( af.removed.expr )
    #                                                  == "AF cleaned" ]
  }

  ### Stage 4: Universal Negatives and Downsampling -----------------
  # clean by selecting universal negative (AF-removed if af.remove performed)
  # also selects brightest positive events
  # and downsamples to speed up subsequent calculations
  # recommended whenever you have an appropriate universal negative

  if ( universal.negative ) {
    if ( !dir.exists( asp$figure.scatter.dir.base ) )
      dir.create( asp$figure.scatter.dir.base )
    if ( !dir.exists( asp$figure.spectral.ribbon.dir ) )
      dir.create( asp$figure.spectral.ribbon.dir )

    # use AF-removed universal negative samples if available
    if ( !is.null( clean.universal.negative ) ) {
      flow.negative <- clean.universal.negative
    } else {
      flow.negative
    }

    # select fluorophore samples to be used
    univ.sample <- flow.control$fluorophore[ ! grepl( "negative",
                                                      flow.control$fluorophore,
                                                      ignore.case = TRUE ) &
                                               flow.negative != FALSE ]
    # probably exclude AF
    univ.sample <- univ.sample[ univ.sample != "AF" ]

    univ.peak.channels <- flow.control$channel[ flow.control$fluorophore
                                                %in% univ.sample ]

    univ.neg.expr <- run.universal.negative( clean.expr, univ.sample,
                                             flow.negative,
                                             flow.control$scatter.parameter,
                                             univ.peak.channels, downsample,
                                             negative.n, positive.n,
                                             spectral.channel, asp,
                                             flow.control.type,
                                             scatter.match )

    # merge in cleaned data
    clean.expr[ names( univ.neg.expr ) ] <- univ.neg.expr

  }

  # downsample only
  if ( !universal.negative & downsample ) {
    if ( !dir.exists( asp$figure.scatter.dir.base ) )
      dir.create( asp$figure.scatter.dir.base )
    if ( !dir.exists( asp$figure.spectral.ribbon.dir ) )
      dir.create( asp$figure.spectral.ribbon.dir )

    # select fluorophore samples to be used
    downsample.sample <- flow.control$fluorophore[ ! grepl( "AF|negative",
                                                            flow.control$fluorophore,
                                                            ignore.case = TRUE ) ]

    downsample.peak.channels <- flow.control$channel[ flow.control$fluorophore
                                                      %in% downsample.sample ]

    downsample.expr <- run.downsample( clean.expr, downsample.sample,
                                       downsample.peak.channels,
                                       negative.n, positive.n, asp$verbose )

    # merge in cleaned data
    clean.expr[ names( downsample.expr ) ] <- downsample.expr
  }

  # merge data and re-establish corresponding factors
  names( clean.expr ) <- flow.sample

  # get maximum number of events per sample to adjust event numbering
  flow.sample.event.number.max <- 0

  for ( fs.idx in 1 : flow.sample.n )
  {
    flow.sample.event.number <- nrow( clean.expr[[ fs.idx ]]  )

    rownames( clean.expr[[ fs.idx ]] ) <- paste( flow.sample[ fs.idx ],
                                                 seq_len( flow.sample.event.number ),
                                                 sep = "_")

    if ( flow.sample.event.number < asp$min.cell.warning.n )
      warning( paste( "\033[31m", "Fewer than", asp$min.cell.warning.n,
                      "gated events in",
                      names( clean.expr )[ fs.idx ], "\033[0m", "\n" ) )

    if ( flow.sample.event.number > flow.sample.event.number.max )
      flow.sample.event.number.max <- flow.sample.event.number
  }

  flow.event.number.width <-
    floor( log10( flow.sample.event.number.max ) ) + 1
  flow.event.regexp <- sprintf( "\\.[0-9]{%d}$", flow.event.number.width )

  # set rownames
  for ( fs.idx in 1 : flow.sample.n )
  {
    flow.sample.event.number <- nrow( clean.expr[[ fs.idx ]]  )
    flow.the.sample <- flow.sample[ fs.idx ]

    flow.the.event <- sprintf( "%s.%0*d", flow.the.sample,
                               flow.event.number.width, 1 : flow.sample.event.number )
    rownames( clean.expr[[ fs.idx ]] ) <- flow.the.event
  }

  clean.expr <- do.call( rbind, clean.expr )

  # set events
  flow.event <- rownames( clean.expr )
  flow.event.n <- length( flow.event )

  flow.event.sample <- sub( flow.event.regexp, "", flow.event )
  flow.event.sample <- factor( flow.event.sample, levels = flow.sample )

  event.type.factor <- flow.sample
  names( event.type.factor ) <- flow.control$control.type
  flow.event.type <- factor( flow.event.sample,
                             levels = event.type.factor,
                             labels = names( event.type.factor ) )

  # store in flow.control
  flow.control$clean.expr <- clean.expr
  flow.control$clean.event.sample <- flow.event.sample
  flow.control$clean.event.type <- flow.event.type

  return( flow.control )
}
