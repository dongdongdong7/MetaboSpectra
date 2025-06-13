# raw_sps: Spectra object from multi samples
# sps: Spectra object for ms1 or ms2 from one sample. It is reconmmended to first perform empty spectra filtering.
# ppm: ppm
# noise: noise
#' @title Convert Spectra to data points data.table
#'
#' @param sps Spectra object for ms1 or ms2 from one sample.
#' It is reconmmended to first perform empty spectra filtering.
#' @param ppm The minimum mz difference in ppm to merge peaks
#' @param noise noise value of spectra
#' @param mz_method Method for calculating the merged peak mz
#' @param intensity_method Method for calculating the merged peak intensity
#'
#' @returns data.table
#' @export
sps2dps <- function(sps, ppm = 10, noise = 100, mz_method = "max_intensity", intensity_method = "sum"){
  # data.table::setDTthreads(threads = 1L)
  # print(data.table::getDTthreads())
  scanIndex <- Spectra::scanIndex(sps)
  rtime <- Spectra::rtime(sps)
  peaksData <- Spectra::peaksData(sps)
  peaksData <- as.list(peaksData)
  peaksData <- batch_slim_peaksMatrix(peaksMatrixList = peaksData,
                                      ppm = ppm,
                                      mz_method = "max_intensity",
                                      intensity_method = "sum")
  # s_time <- Sys.time()
  dps <- data.table::rbindlist(
    lapply(seq_along(peaksData), function(i) {
      dt <- as.data.table(peaksData[[i]])
      dt[, scanIndex := scanIndex[i]]
      dt
    })
  )
  # print(Sys.time() - s_time)
  dps[intensity > noise, ]
}
