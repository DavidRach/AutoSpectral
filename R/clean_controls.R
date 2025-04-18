# clean_controls.r

#' @title clean.controls
#'
#' @description
#' A four-part function to clean single color controls in order to extract
#'     fluorophore signatures. Any part can be run independently.
#'     Stage 1: PeacoQC to eliminate flow artefacts.
#'     Stage 2: Trimming to eliminate extreme events.
#'     Stage 3: AF noise removal using PCA unmixing on matching unstained (cells only).
#'     Stage 4: brightest event selection from positive, universal negative from
#'     matching negative and downsampling to speed up rlm spectra optimization.
#'
#' @param flow.control Prepare using define.flow.control. List of the data and
#'     essential information about the cytometer and data structure.
#' @param asp The AutoSpectral parameter list. Prepare using get.autospectral.param.
#'     A list of essential parameters.
#' @param time.clean Logical, default is FALSE. Determines whether to run PeacoQC to
#'     remove time-based inconsistencies in the controls.
#' @param trim Logical, default is FALSE. Determines whether to remove extreme
#'     events (positive and negative) from controls. Extent of trimming is determined
#'     by trim.factor. Include trim.factor if setting trim to TRUE.
#' @param trim.factor Numeric. Default is asp$rlm.trim.factor. Include trim.factor
#'     if setting trim to TRUE.
#' @param af.remove Logical, default is FALSE. Determines whether to remove
#'     intrusive autofluorescence contamination from cell controls using PCA-based
#'     identification and gating. Requires universal negatives to be defined in the
#'     control file and subsequently in flow.control via define.flow.control.
#' @param universal.negative Logical, default is TRUE. Determines whether to
#'     use a universal negative sample as the negative for spectral extraction.
#'     Requires universal negatives to be defined in the
#'     control file and subsequently in flow.control via define.flow.control.
#' @param downsample Logical, default is TRUE. Determines whether to reduce cell
#'     and bead control events to speed up processing.
#' @param negative.n Number of events to include in the downsampled negative
#'     population. Default is asp$negative.n.
#' @param positive.n Number of events to include in the downsampled positive
#'     population. Default is asp$positive.n
#' @param scatter.match Logical, default is TRUE. Determines whether to select
#'     negative events based scatter profiles matching the positive events. Scatter
#'     matching defines a region of FSC and SSC determined by the distribution of
#'     selected positive events (see positive.n) and defines a gate for this region.
#' @param refactor Logical, default is TRUE. Reassigns event indices after cleaning
#'     steps, setting up fresh, consistent indices.
#' @param return.separately Logical, default is FALSE. Provides the option to
#'     return the cleaned expression data as an output as well as storing the clean
#'     data in flow.control.
#'
#' @return No returns unless return.separately is set to TRUE. Output is stored
#'     in flow.control so that flow.control can be used subsequenlty as a complete
#'     data set.
#' @export


clean.controls <- function( flow.control, asp,
                            time.clean = FALSE,
                            trim = FALSE, trim.factor = NULL,
                            af.remove = FALSE,
                            universal.negative = TRUE, downsample = TRUE,
                            negative.n = asp$negative.n, positive.n = asp$positive.n,
                            scatter.match = TRUE,
                            refactor = TRUE,
                            return.separately = FALSE ) {

  clean.expr <- flow.control$expr.data

  flow.sample <- flow.control$sample

  flow.sample.n <- length( flow.sample )

  clean.event.sample <- flow.control$event.sample

  flow.negative <- flow.control$universal.negative

  flow.control.type <- flow.control$control.type

  clean.universal.negative <<- NULL

  # split expression data by sample into a list
  clean.expr <- lapply( flow.sample, function( fs ) {
    clean.expr[ clean.event.sample == fs, ]
  })

  names( clean.expr ) <- flow.sample

  clean.event.type <- flow.control$event.type

  spectral.channel <- flow.control$spectral.channel

  all.channels <- flow.control$scatter.and.channel.original

  ### Stage 1: Use PeacoQC to clean up flow  -----------------
  # clean based on irregularities in flow by TIME parameter
  # this can be slow and may not show much effect

  if ( time.clean ) {

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

    # identify universal negative cell samples
    univ.neg <- unique( flow.negative )
    univ.neg <- univ.neg[ !is.null( univ.neg ) ]
    univ.neg <- univ.neg[ !is.na( univ.neg ) ]
    univ.neg <- univ.neg[ univ.neg != FALSE ]

    check.critical( length( univ.neg ) > 0,
                    "Error: no cell-based universal negative samples could be identified.
                    To perform autofluorescence removal, you must specify a universal negative in the fcs_control_file." )


    univ.neg.sample <- flow.control$fluorophore[ flow.control.type == "cells" &
                                                 flow.sample %in% univ.neg ]

    check.critical( length( univ.neg.sample ) > 0,
                    "Error: no cell-based universal negative samples could be identified.
                    To perform autofluorescence removal, you must specify a universal negative in the fcs_control_file." )


    # select cell-based single-stained samples to be used
    # must be cells and must have a corresponding universal negative
    # note: last exclusion step can be modified to include AF control
    af.removal.sample <- flow.sample[ flow.control.type == "cells" &
                                        flow.negative %in% univ.neg ]
    # & !( flow.control$fluorophore %in% univ.neg.sample )

    # check that is not length 0 and stop if is
    check.critical( length( af.removal.sample ) > 0,
                    "Error: no cell-based samples with a corresponding universal negative could be identified.
                    To perform autofluorescence removal, you must specify a universal negative in the fcs_control_file." )

    # get af artefacts for each negative control
    univ.neg.af.artefacts <- run.af.artefact.id( clean.expr, univ.neg,
                                                 spectral.channel,
                                                 asp )

    # remove identified AF from single-color controls
    af.removed.expr <- run.af.removal( clean.expr, af.removal.sample,
                                       univ.neg.af.artefacts,
                                       spectral.channel,
                                       flow.negative,
                                       asp )

    # if AF is among cleaned controls
    # rename AF in cleaned controls to AF cleaned
    names( af.removed.expr )[ names( af.removed.expr ) == "AF" ] <- "AF cleaned"

    # replace AF in universal negative with AF cleaned
    clean.universal.negative <- flow.negative
    clean.universal.negative[ clean.universal.negative == "AF" ] <- "AF cleaned"

    flow.control$clean.universal.negative <<- clean.universal.negative

    # store in all except AF and AF cleaned
    # get names in case of only cleaning certain samples
    cleaned.samples <- names( af.removed.expr )[ names( af.removed.expr ) != "AF cleaned" ]

    clean.expr[ cleaned.samples ] <- af.removed.expr[ cleaned.samples ]

    # store AF cleaned in a new slot
    clean.expr[ "AF cleaned" ] <- af.removed.expr[ names( af.removed.expr )
                                                      == "AF cleaned" ]
  }

  ### Stage 4: Universal Negatives and Downsampling -----------------
  # clean by selecting universal negative (AF-removed if af.remove performed)
  # also selects brightest positive events
  # and downsamples to speed up subsequent calculations
  # recommended whenever you have an appropriate universal negative

  if ( universal.negative ) {

    # use AF-removed universal negative samples if available
    if ( !is.null( clean.universal.negative ) ) {
      flow.negative <- clean.universal.negative
    } else {
      flow.negative
    }

    # select fluorophore samples to be used
    univ.sample <- flow.control$fluorophore[ ! grepl( "negative",
                                                      flow.control$fluorophore,
                                                      ignore.case = TRUE ) ]

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

  if ( !universal.negative & downsample ) {

    # select fluorophore samples to be used
    downsample.sample <- flow.control$fluorophore[ ! grepl( "AF|negative",
                                                            flow.control$fluorophore,
                                                            ignore.case = TRUE ) ]

    downsample.peak.channels <- flow.control$channel[ flow.control$fluorophore
                                                %in% downsample.sample ]

    downsample.expr <- run.downsample( clean.expr, downsample.sample,
                                       downsample.peak.channels,
                                       negative.n, positive.n, asp )

    # merge in cleaned data
    clean.expr[ names( downsample.expr ) ] <- downsample.expr
  }

  if ( refactor ) {

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

      if ( flow.sample.event.number < 500 )
        cat( paste( "Warning! Fewer than 500 gated events in", names( clean.expr )[ fs.idx ] ) )

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
    # yes, this is storing in a global variable from within a function
    # easier for user
    flow.control$clean.expr <<- clean.expr

    flow.control$clean.event.sample <<- flow.event.sample

    flow.control$clean.event.type <<- flow.event.type
  }

  if ( return.separately )
    return( clean.expr )
}











