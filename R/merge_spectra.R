# Merge multiple spectra
# 250408
# Barry Song

#' @rdname merge_spectra
#' @title Merge spectra
#' @description
#' When we are comparing multiple spectra that have inconsistent mz precision, we often need to
#' first change the mz precision to consistent. In addition, there are often multiple spectra for
#' one metabolite in the database, and combining these spectra is useful for metabolite identification.
#'
#' @param spMat spectra matrix with mz and intensity
#' @param digits mz digits
#' @export
#' @examples
#' hmdbMsTb <- pubmsR::load_hmdbMsTb()
#' spMatList <- lapply(1:nrow(hmdbMsTb[hmdbMsTb$accession == "HMDB0000062", ]), function(i) {
#'   MetaboSpectra::get_spMat(hmdbMsTb[hmdbMsTb$accession == "HMDB0000062", ][i, ])
#' })
#' round_spMat(spMat = spMatList[[2]], digits = 0)
#' round_spMat(spMat = spMatList[[2]], digits = 1)
#' round_spMat(spMat = spMatList[[2]], digits = 2)
round_spMat <- function(spMat, digits = 0){
  spMat[, "mz"] <- round(spMat[, "mz"], digits = digits)
  spDf <- as.data.frame(spMat)
  spDf <- aggregate(intensity ~ mz, data = spDf, FUN = sum)
  spMat <- as.matrix(spDf)
  clean_spMat(spMat = spMat, noise_threshold = 0)
}

#' @rdname merge_spectra
#' @param spMatList spMat list
#' @param digits mz digits
#' @param method sum, median, mean
#'
#' @returns spMat
#' @export
#'
#' @examples
#' merge_spectra(spMatList = spMatList, method = sum)
merge_spectra <- function(spMatList, digits = 0, method = sum){
  spMatList <- lapply(spMatList, function(x) {
    round_spMat(spMat = x, digits = digits)
  })
  spDfList <- lapply(spMatList, function(x) {
    as.data.frame(x)
  })
  spDf <- purrr::list_rbind(spDfList)
  spDf <- aggregate(intensity ~ mz, data = spDf, FUN = method)
  spMat <- as.matrix(spDf)
  clean_spMat(spMat = spMat, noise_threshold = 0)
}
