# fluorophore_database.r

#' @title Fluorophore Database
#'
#' @description Information about fluorophores and their detection on various
#'     cytometers.
#' @format A data frame with columns:
#' \describe{
#'  \item{fluorophore}{Fluorophore name}
#'  \item{synonym1}{Fluorophore synonym 1}
#'  \item{synonym2}{Fluorophore synonym 2}
#'  \item{synonym3}{Fluorophore synonym 3}
#'  \item{channel.Aurora}{Peak channels for Fluorophores on the 5L Aurora}
#'  \item{channel.ID7000}{Peak channels for Fluorophores on the 5L ID7000}
#'  \item{channel.s8}{Peak channels for Fluorophores on the DiscoverS8}
#'  \item{channel.opteon}{Peak channels for Fluorophores on the Opteon}
#'  \item{excitation.laser}{Excitation laser}
#'  \item{nominal.wavelength}{Numeric: nominal peak emission wavelength}
#'  \item{is.viability}{Logical: is the fluorophore a viability dye}
#' }
