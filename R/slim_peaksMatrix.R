# Slim peaks matrix (spectrum) based on ppm
# (Merge adjacent peaks within the ppm range)
# 250605
# Barry Song

#' @rdname slim_peaksMatrix
#' @title Slim peaks matrix (spectrum)
#' @description
#' Merge adjacent peaks within the ppm range
#'
#' @param peaksMatrix  Peaks matrix (spectrum) from ms1 or ms2; also named spMat. It must be centroided
#' @param ppm The minimum mz difference in ppm to merge peaks
#' @param mz_method Method for calculating the merged peak mz
#' @param intensity_method Method for calculating the merged peak intensity
#'
#' @returns peaksMatrix
#' @export
#'
#' @examples
#' set.seed(123)
#' spMat <- cbind(
#'   mz = sort(rnorm(1000, mean = 100, sd = 0.1)),
#'   intensity = rgamma(1000, shape = 2, scale = 1)
#' )
#' system.time({
#' lapply(1:1000, function(i) {
#'   slim_peaksMatrix(peaksMatrix = spMat, ppm = 30)
#'  })
#' })
slim_peaksMatrix <- function(peaksMatrix, ppm = 5,
                                   mz_method = c("median", "mean", "max_intensity"),
                                   intensity_method = c("sum", "mean", "max")) {
  mz_method <- match.arg(mz_method)
  intensity_method <- match.arg(intensity_method)

  .slim_peaksMatrix_rcpp(peaksMatrix = peaksMatrix, ppm = ppm, mz_method = mz_method, intensity_method = intensity_method)
}

#' @rdname slim_peaksMatrix
#' @returns peaksMatrixList
#' @export
#'
#' @examples
#' spMatList <- lapply(1:1000, function(i) {spMat})
#' system.time({
#' batch_slim_peaksMatrix(peaksMatrixList = spMatList)
#' })
batch_slim_peaksMatrix <- function(peaksMatrixList, ppm = 5,
                                   mz_method = c("median", "mean", "max_intensity"),
                                   intensity_method = c("sum", "mean", "max")){
  mz_method <- match.arg(mz_method)
  intensity_method <- match.arg(intensity_method)

  if(!is.list(peaksMatrixList)) stop("peaksMatrixList should be a list")

  .batch_slim_peaksMatrix_rcpp(peaksMatrixList = peaksMatrixList, ppm = ppm, mz_method = mz_method, intensity_method = intensity_method)
}
