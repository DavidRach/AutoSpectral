# unmix_folder.r

#' @title Unmix All FCS Files in a Directory
#'
#' @description
#' This function unmixes all FCS files in a specified directory using the
#' provided spectra and method, and saves the unmixed FCS files to an output
#' directory of the user's choice.
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#'
#' @param fcs.dir Directory containing FCS files to be unmixed.
#' @param spectra Matrix containing spectra information.
#' @param asp The AutoSpectral parameter list.
#' Prepare using `get.autospectral.param`
#' @param flow.control A list containing flow cytometry control parameters.
#' @param method A character string specifying the unmixing method to use.
#' The default is `Automatic`, which uses `AutoSpectral` for AF extraction if
#' af.spectra are provided and automatically selects `OLS` or `WLS` depending
#' on which is normal for the given cytometer in `asp$cytometer`. This means
#' that files from the ID7000, A8 and S8 will be unmixed using `WLS` while
#' others will be unmixed with `OLS`. Any option can be set manually.
#' Manual options are `OLS`, `WLS`, `AutoSpectral`, `Poisson` and `FastPoisson`.
#' Default is `OLS`. `FastPoisson` requires installation of `AutoSpectralRcpp`.
#' @param weights Optional numeric vector of weights: one per fluorescent
#' detector. Default is `NULL`, in which case weighting will be done by
#' channel means. Only used for `WLS`
#' @param af.spectra Spectral signatures of autofluorescences, normalized
#' between 0 and 1, with fluorophores in rows and detectors in columns. Prepare
#' using `get.af.spectra`. Required for `AutoSpectral` unmixing. Default is
#' `NULL` and will thus provoke failure if no spectra are provided and
#' `AutoSpectral` is selected.
#' @param output.dir Directory to save the unmixed FCS files
#' (default is asp$unmixed.fcs.dir).
#' @param file.suffix A character string to append to the output file name.
#' Default is `NULL`
#' @param include.raw Logical indicating whether to include raw data in the
#' written FCS file. Default is `FALSE`
#' @param include.imaging Logical indicating whether to include imaging data in
#' the written FCS file: relevant for S8 and A8. Default is `FALSE`
#' @param divergence.threshold Numeric. Used for `FastPoisson` only. Threshold
#' to trigger reversion towards WLS unmixing when Poisson result diverges.
#' Default is `1e4`
#' @param divergence.handling String. How to handle divergent cells from Poisson
#' IRLS. Options are `NonNeg`, in which case non-negativity will be enforced,
#' `WLS`, where values will revert to the WLS initial unmix or `Balance`,
#' where `WLS` and `NonNeg` results will be averaged. Default is `Balance`
#' @param balance.weight Numeric. Weighting to average non-convergent cells.
#' Used for `Balance` option under `divergence.handling`. Default is `0.5`
#'
#' @return None. Saves the unmixed FCS files to the specified output directory.
#'
#' @export

unmix.folder <- function( fcs.dir, spectra, asp, flow.control,
                          method = "Automatic",
                          weights = NULL,
                          af.spectra = NULL,
                          output.dir = NULL,
                          file.suffix = NULL,
                          include.raw = FALSE,
                          include.imaging = FALSE,
                          divergence.threshold = 1e4,
                          divergence.handling = "Balance",
                          balance.weight = 0.5 ){

  if ( is.null( output.dir ) ){
    output.dir <- asp$unmixed.fcs.dir
  }

  files.to.unmix <- list.files( fcs.dir, pattern = ".fcs", full.names = TRUE )

  # set up parallel processing
  if ( asp$parallel && method != "Poisson" && method != "FastPoisson" ){
    plan( multisession, workers = asp$worker.process.n )
    options( future.globals.maxSize = asp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  lapply.function( files.to.unmix, FUN = unmix.fcs, spectra, asp, flow.control,
                   method, weights, af.spectra, output.dir, file.suffix, include.raw,
                   include.imaging, divergence.threshold, divergence.handling,
                   balance.weight )
}
