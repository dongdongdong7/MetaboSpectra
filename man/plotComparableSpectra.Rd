% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_function.R
\name{plotComparableSpectra}
\alias{plotComparableSpectra}
\title{plotComparableSpectra}
\usage{
plotComparableSpectra(
  spMat_query,
  spMat_lib,
  num = 10,
  ms2_tolerance_in_da = -1,
  ms2_tolerance_in_ppm = 5
)
}
\arguments{
\item{spMat_query}{Query spectra matrix}

\item{spMat_lib}{Library spectra matrix}

\item{num}{max peak number}

\item{ms2_tolerance_in_da}{The MS2 tolerance in Da, set to -1 to disable}

\item{ms2_tolerance_in_ppm}{The MS2 tolerance in ppm, set to -1 to disable}
}
\value{
A ggplot object.
}
\description{
Compare the two mass spectra
}
\examples{
mz <- c(21.3303, 40.1326, 86.3402, 138.3291, 276.5717, 276.5836)
intensity <- c(100, 1300, 4030, 10000, 31600, 1000)
standardRow1 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
spMat1 <- get_spMat(standardRow1)
mz <- c(21.3002, 40.1128, 86.3201, 138.3214, 276.5312)
intensity <- c(100, 1300, 4030, 10000, 31600)
standardRow2 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
spMat2 <- get_spMat(standardRow2)
plotComparableSpectra(spMat1, spMat2, tol_ppm2 = 200, tol_da2 = -1)
plotComparableSpectra(spMat1, spMat2, tol_ppm2 = -1, tol_da2 = 0.02)
}
