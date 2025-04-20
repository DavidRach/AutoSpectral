
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AutoSpectral

<!-- badges: start -->
<!-- badges: end -->

The goal of AutoSpectral is to provide you with the best possible
spectral signatures of the fluorophores in your single-stained controls.
Whether or not these accurately model your fully stained samples will
depend on what you’ve chosen to use for the controls, how they were
prepared and other factors such as any divergence in handling between
samples and controls as well as machine variation.

More to the point, AutoSpectral is intended to make working with messy
cell-based controls as easy as compensation beads. This should give you
better accuracy and precision in your spectral definition and thus in
your unmixing.

## Installation

You can install the development version of AutoSpectral from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("DrCytometer/AutoSpectral")
```

## Example

This is a basic example of the workflow, using samples from the ID7000:

``` r
library(AutoSpectral)

# define the location of the single-stained control files
control.dir <- "./Cell controls"

# get the parameters for your cytometer
asp <- get.autospectral.param( cytometer = "id7000", figures = TRUE )

# optionally, create a control file
# after creating it, manually edit to ensure it is correct
create.control.file( control.dir, asp )

# load the edited control file
control.def.file <- "fcs_control_file.csv"

# adjustments to automated gating
# this will be covered more extensively elsewhere
# in this case, the cells on scatter are very small, so we modify the target
asp$gate.bound.density.max.target.cells <- 0

# load the control files, gate, prepare for spectral extraction
# this is usually the slowest step
flow.control <- define.flow.control( control.dir, control.def.file, asp )

# for speed and accuracy, cleaning the controls by using scatter-matched universal
# negatives are recommended for both beads and cells
clean.controls( flow.control, asp, universal.negative = TRUE, downsample = TRUE,
                scatter.match = TRUE )

# extract the clean spectra
univ.neg.spectra <- get.fluorophore.spectra( flow.control, asp, 
                                            use.clean.expr = TRUE,
                                            plot.prefix = "Universal Negative Cells" )

# by default AutoSpectral extracts autofluorescence as a spectrum
# to remove this for unmixing:
no.af.spectra <- univ.neg.spectra[ ! rownames( univ.neg.spectra ) == "AF", ]

## unmixing
# set the location of the raw files to be unmixed
# sample.dir <- "./Raw samples"

# perform unmixing, here we will unmix all fcs files in the folder using wls
unmix.folder( sample.dir, flow.control$spectral.channel, no.af.spectra, asp,
              method = "wls" )
```
