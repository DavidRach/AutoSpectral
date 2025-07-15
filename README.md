
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AutoSpectral

<!-- badges: start -->

<!-- badges: end -->

AutoSpectral is AutoSpill updated for the spectral flow era.

It is also pretty complex and newly released, so there will be bugs.
Sorry. Please read the help pages and articles before submitting bug
reports. There’s a lot of info there.

The goal of AutoSpectral is to provide you with the best possible
spectral signatures of the fluorophores in your single-stained controls.
Whether or not these accurately model your fully stained samples will
depend on what you’ve chosen to use for the controls, how they were
prepared and other factors such as machine condition and any divergence
in handling between samples and controls.

More to the point, AutoSpectral is intended to make working with messy
cell-based controls as easy as compensation beads. This should give you
better accuracy and precision in your spectral definition and thus in
your unmixing.

Plus, you can extract each cell’s individual autofluorescent background
in a manner specific to that cell, producing better unmixing with less
spread.

At the moment, the following cytometers are supported:

- Cytek Aurora/Northern Lights (“aurora”)
- Sony ID7000 (“id7000”)
- BD FACSDiscoverS8 (“s8”)
- BD FACSDiscoverA8 (“a8”)
- Agilent NovoCyte Opteon (“opteon”)

There will likely be some unresolved issues with plotting data from
certain configurations of the ID7000 and Aurora. Please note that FCS
3.2 files from the S8 and A8 cytometers are not fully supported in
flowCore. You may receive warnings, but things should still work.

If you want to use data from another cytometer and are wiling to provide
files for establishing the workflow, contact the author/maintainer.

This work has received funding from the KU Leuven C1 program, the
European Union’s Horizon 2020 research and innovation programme under
grant agreement No 874707 (EXIMIOUS), Wellcome Investigator Award,
222442/A/21/Z and UKRI Proactive Vaccinology Award, MR/Y004450/1
(IMMPROVE).

AutoSpectral is provided under a GPL3 licence for academic use. If you
wish to use AutoSpectral in a non-academic setting, please contact
Adrian Liston regarding licencing options.

## Installation

You can install the development version of AutoSpectral from
[GitHub](https://github.com/) with:

``` r

# Optional: Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowWorkspace", "flowCore", "PeacoQC", "Biobase"))

# You'll need devtools or remotes to install from GitHub.
# install.packages("devtools")
devtools::install_github("DrCytometer/AutoSpectral")
```

## Example

This is a basic example of the workflow, using samples from the ID7000:

``` r
library( AutoSpectral )

# Define the location of the single-stained control files.
control.dir <- "./Cell controls"

# Get the parameters for your cytometer.
# Supported cytometers include "aurora", "id7000", "a8", "s8" and "opteon".
asp <- get.autospectral.param( cytometer = "id7000", figures = TRUE )

# Optionally, create a control file (see article on this).
# After creating it, manually edit to ensure it is correct.
create.control.file( control.dir, asp )

# Locate the edited control file.
control.def.file <- "fcs_control_file.csv"

# check the control file for errors
control.file.errors <- check.control.file(control.dir, control.def.file, asp)

## Adjustments to automated gating
# This will be covered more extensively elsewhere.
# in this case, the cells on scatter are very small, so we modify the target.
asp$gate.bound.density.max.target.cells <- 0

# Load the control files, gate, prepare for spectral extraction.
# This is usually the slowest step.
flow.control <- define.flow.control( control.dir, control.def.file, asp )

# For speed and accuracy, the recommended approach is to clean the controls
# using separate (universal) negative samples for both beads and cells.
# In doing so, we select cells with matching autofluorescent background and
# restrict the calculations to the brightest events.
clean.controls( flow.control, asp, universal.negative = TRUE, downsample = TRUE,
                scatter.match = TRUE )

# Extract the clean spectra
univ.neg.spectra <- get.fluorophore.spectra( flow.control, asp, 
                                            use.clean.expr = TRUE,
                                            plot.prefix = "Universal Negative Cells" )

# By default AutoSpectral extracts autofluorescence as a spectrum.
# This is comparable to "Autofluorescence as a Fluorescent Tag" in SpectroFlo
# More on this later.
# To remove this AF parameter for unmixing:
no.af.spectra <- univ.neg.spectra[ ! rownames( univ.neg.spectra ) == "AF", ]

## Unmixing
# Set the location of the raw files to be unmixed.
sample.dir <- "./Raw samples"

# Perform unmixing, here we will unmix all fcs files in the folder using 
# Weighted least squares (WLS or in Sony lingo, WLSM).
unmix.folder( sample.dir, no.af.spectra, asp, flow.control, method = "WLS" )
```
