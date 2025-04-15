# create_directory.r

#' @title create.directory
#'
#' @description
#' Creates figure and table directories.
#'
#' @param asp The AutoSpectral parameter list defined using get.autospectral.param.
#' @return No returns. Creates directories in the current working directory.
#' @export

create.directory <- function( asp )
{

    figure.dir <- c(
      asp$figure.scatter.dir.base,
      asp$figure.af.dir,
      asp$figure.gate.dir,
      asp$figure.peacoqc.dir,
      asp$figure.clean.control.dir,
      asp$figure.spectral.ribbon.dir,
      asp$figure.convergence.dir,
      asp$figure.spectra.dir,
      asp$figure.slope.error.dir,
      asp$figure.skewness.dir,
      asp$figure.similarity.heatmap.dir
    )

    table.dir <- c(
      asp$table.convergence.dir,
      asp$table.spectra.dir,
      asp$table.slope.error.dir,
      asp$table.skewness.dir
    )

    for ( ftd in c( figure.dir, table.dir, asp$unmixed.fcs.dir ) )
        if ( ! is.null( ftd ) && ! file.exists( ftd ) )
            dir.create( ftd, recursive = TRUE )

}

