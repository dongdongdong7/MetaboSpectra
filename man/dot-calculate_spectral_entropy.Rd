% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.calculate_spectral_entropy}
\alias{.calculate_spectral_entropy}
\title{Calculate spectral entropy of a spectrum}
\usage{
.calculate_spectral_entropy(peaks)
}
\arguments{
\item{peaks}{A matrix of peaks, with two columns: m/z and intensity.}
}
\value{
A double value of spectral entropy.
}
\description{
Calculate spectral entropy of a spectrum
}
\examples{
mz <- c(100.212, 300.321, 535.325)
intensity <- c(37.16, 66.83, 999.0)
peaks <- matrix(c(mz, intensity), ncol = 2, byrow = FALSE)
calculate_spectral_entropy(peaks)

}
