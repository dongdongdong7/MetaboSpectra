#' Pipe operator
#' @name %>%
#' @export
#' @importFrom magrittr %>%
NULL

## Rcpp
## usethis::use_rcpp()
#' @useDynLib MetaboSpectra
#' @importFrom Rcpp evalCpp
NULL

## usethis namespace: start
#' @useDynLib MetaboSpectra, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @title Get a spectrum matrix
#' @description
#' Get a spectrum matrix from a tibble.
#'
#' @param standardRow
#' A tibble which need to include two columns mz and intensity.
#'
#' @return `matrix()` with two columns mz and intensity
#' @export
#'
#' @examples
#' mz <- c(40.1320, 86.3400, 138.3290, 276.5710)
#' intensity <- c(1300, 4030, 10000, 31600)
#' standardRow <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' get_spMat(standardRow)
get_spMat <- function(standardRow){
  if(nrow(standardRow) != 1) stop("standardRow must has 1 row!")
  mz <- standardRow$mz[[1]];intensity <- standardRow$intensity[[1]]
  if(length(mz) != length(intensity)) stop("length of mz and intensity don't match")
  spMat <- matrix(c(mz, intensity), ncol = 2,
                  dimnames = list(NULL, c("mz", "intensity")))
  spMat <- spMat[complete.cases(spMat), , drop = FALSE] # remove NA
  spMat <- spMat[order(spMat[, "mz"]), , drop = FALSE]
  return(spMat)
}

#' @title Clean spectrum matrix
#' @description
#' Clean a spMat includes removing noise, merging similar peaks, and intensity normaliaztion.
#'
#' @param spMat The matrix from the mass spectrum needs to contain two columns, mz and intensity.
#' @param min_ms2_difference_in_da The minimum mz difference in Da to merge peaks, set to -1 to disable, any two peaks with mz difference < min_ms2_difference_in_da will be merged
#' @param min_ms2_difference_in_ppm The minimum mz difference in ppm to merge peaks, set to -1 to disable, any two peaks with mz difference < min_ms2_difference_in_ppm will be merged
#' @param min_mz The minimum mz value to keep, set to -1 to disable
#' @param max_mz The maximum mz value to keep, set to -1 to disable
#' @param noise_threshold The noise threshold, set to -1 to disable, all peaks have intensity < noise_threshold * max_intensity will be removed
#' @param max_peak_num The maximum number of peaks to keep, set to -1 to disable
#' @param scale_int Normalized target value, set to -1 to disable
#' @param normalize_intensity Whether to normalize the intensity to sum to 1
#'
#' @return spMat.
#'
#' @export
#'
#' @examples
#' mz <- c(21.3300, 40.1320, 86.3400, 138.3290, 276.5710, 276.5830)
#' intensity <- c(100, 1300, 4030, 10000, 31600, 1000)
#' standardRow <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat <- get_spMat(standardRow)
#' spMat <- clean_spMat(spMat)
#'
clean_spMat <- function(spMat,
                        min_ms2_difference_in_da = -1,
                        min_ms2_difference_in_ppm = 5,
                        min_mz = -1, max_mz = -1,
                        noise_threshold = 0.01,
                        max_peak_num = -1,
                        scale_int = 100,
                        normalize_intensity = FALSE){
  if(min_ms2_difference_in_da > 0 & min_ms2_difference_in_ppm > 0) min_ms2_difference_in_ppm <- -1
  spMat <- clean_spectrum(spMat = spMat,
                          min_mz = min_mz, max_mz = max_mz,
                          min_ms2_difference_in_da = min_ms2_difference_in_da,
                          min_ms2_difference_in_ppm = min_ms2_difference_in_ppm,
                          noise_threshold = noise_threshold,
                          max_peak_num = max_peak_num,
                          normalize_intensity = normalize_intensity)
  if(!normalize_intensity){
    if(scale_int != -1) spMat[, "intensity"] <- scale_int *spMat[, "intensity"] / max(spMat[, "intensity"])
  }
  return(spMat)
}
#' @title Compare spMat entropy
#' @description Compare two spMat using entropy. Two spMat must be cleaned before calculate entropy.
#' The intensity of spMat must sum to 1.
#'
#' @param x Query spMat
#' @param y Library spMat
#' @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 to disable
#' @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 to disable
#' @param unweighted Whether to calculate the unweighted entropy similarity
#'
#' @return The entropy similarity
#' @export
#'
#' @examples
#' mz <- c(21.3300, 40.1320, 86.3400, 138.3290, 276.5710, 276.5830)
#' intensity <- c(100, 1300, 4030, 10000, 31600, 1000)
#' standardRow1 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat1 <- get_spMat(standardRow1)
#' mz <- c(21.3000, 40.1120, 86.3200, 138.3210, 276.5310)
#' intensity <- c(100, 1300, 4030, 10000, 31600)
#' standardRow2 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat2 <- get_spMat(standardRow2)
#' # Normalize intensity (sum to 1)
#' spMat1 <- clean_spMat(spMat1, normalize_intensity = TRUE)
#' spMat2 <- clean_spMat(spMat2, normalize_intensity = TRUE)
#' compare_spMat_entropy(x = spMat1, y = spMat2, ms2_tolerance_in_da = 0.02)
compare_spMat_entropy <- function(x, y,
                                  ms2_tolerance_in_da = -1,
                                  ms2_tolerance_in_ppm = 5,
                                  unweighted = FALSE){
  if(ms2_tolerance_in_da > 0 & ms2_tolerance_in_ppm > 0) ms2_tolerance_in_ppm <- -1
  x <- x[order(x[, 1]), , drop = FALSE]
  y <- y[order(y[, 1]), , drop = FALSE]
  if(unweighted){
    score <- calculate_unweighted_entropy_similarity(spMat1 = x, spMat2 = y,
                                                     ms2_tolerance_in_da = ms2_tolerance_in_da,
                                                     ms2_tolerance_in_ppm = ms2_tolerance_in_ppm)
  }else{
    score <- calculate_entropy_similarity(spMat1 = x, spMat2 = y,
                                          ms2_tolerance_in_da = ms2_tolerance_in_da,
                                          ms2_tolerance_in_ppm = ms2_tolerance_in_ppm)
  }
  return(score)
}

#' @title compare_spMat_ndotproduct
#' @description
#' The dot product method is used to compare the similarity of two mass spectrum matrices.
#'
#' @param x
#' spMat_query.
#' @param y
#' spMat_lib.
#' @param joinpeak inner, left, right, outer; More details please see [Spectra::joinPeaks]
#' @param m `numeric`, weighting for the first column of x and y ("m/z")
#' @param n `numeric`, weighting for the second column of x and y ("intensity")
#' @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 to disable
#' @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 to disable
#'
#' @return
#' A value.
#' @export
#'
#' @examples
#' mz <- c(21.3300, 40.1320, 86.3400, 138.3290, 276.5710, 276.5830)
#' intensity <- c(100, 1300, 4030, 10000, 31600, 1000)
#' standardRow1 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat1 <- get_spMat(standardRow1)
#' mz <- c(21.3000, 40.1120, 86.3200, 138.3210, 276.5310)
#' intensity <- c(100, 1300, 4030, 10000, 31600)
#' standardRow2 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat2 <- get_spMat(standardRow2)
#' # Spectra::joinPeaks(x = spMat1, y = spMat2, type = "inner", tolerance = 0, ppm = 200)
#' compare_spMat_ndotproduct(x = spMat1, y = spMat2, joinpeak = "outer", ms2_tolerance_in_da = 0.03)
compare_spMat_ndotproduct <- function(x, y, joinpeak = "inner",
                                      m = 0L, n = 0.5,
                                      ms2_tolerance_in_da = -1,
                                      ms2_tolerance_in_ppm = 5){
  if(ms2_tolerance_in_da >0 & ms2_tolerance_in_ppm > 0) ms2_tolerance_in_ppm <- -1
  if(ms2_tolerance_in_da < 0) ms2_tolerance_in_da <- 0
  if(ms2_tolerance_in_ppm < 0) ms2_tolerance_in_ppm <- 0
  x <- x[order(x[, 1]), , drop = FALSE]
  y <- y[order(y[, 1]), , drop = FALSE]

  tmpList <- Spectra::joinPeaks(x = x, y = y,type = joinpeak, tolerance = ms2_tolerance_in_da, ppm = ms2_tolerance_in_ppm)
  x <- tmpList$x;y <- tmpList$y
  score <- MsCoreUtils::ndotproduct(x = x, y = y, m = m, n = n)
  return(score)
}

# TODO: whether to keep searchLib function?

#' @title searchLib_entropy
#' @description
#' The mass spectrometry database was screened using the spectral entropy algorithm.
#'
#' @param standardInput
#' A tibble, including six columns of information,
#' feature: name of the feature,
#' precursorMz: parent ion mz,
#' rt: retention time,
#' adduct: adduct form,
#' mz: mz of secondary mass spectrum,
#' intensity: intensity of secondary mass spectrum;
#' Both mz and intensity are stored in tibble as a list.
#' @param lib
#' From any of the publicMs2List databases.
#' @param st
#' Threshold for scoring.
#' @param tol_da1
#' Under the Da unit, the tolerance of the same parent ion mz is searched.
#' @param tol_ppm1
#' Under the ppm unit, the tolerance of the same parent ion mz is searched.
#' @param tol_da2
#' Under the Da unit, two peaks are considered the tolerance of one peak.
#' @param tol_ppm2
#' Under the ppm unit, two peaks are considered the tolerance of one peak.
#' @param predicted
#' Whether to search only predictive mass spectra;All, Yes, No.
#' @param thread
#' The number of threads running the function.
#' @param min_mz
#' Minimum value of mz.
#' @param max_mz
#' Maximum value of mz.
#' @param noise_threshold
#' Peak intensity below noise_threshold * max(intensity) will be considered noise.
#' @param max_peak_num
#' The number of peaks in the mass spectrum.
#' @details
#' min_mz, max_mz, noise_threshold, max_peak_num, scale_int are parameter of function clean_spectra.
#'
#' @return
#' searchRes_entropy, a list.
#'
#' @examples
#' data("standardInput", package = "MetaboSpectra")
#' load("publicMs2List.RData") # please see https://github.com/dongdongdong7/MetaboLib.ms2 and download publicMs2List.RData
#' searchRes_entropy <- searchLib_entropy(standardInput = standardInput, lib = publicMs2List$hmdb, thread = 8, st = 0.5)
# searchLib_entropy <- function(standardInput, lib, st = 0.8,
#                               tol_da1 = 0.01, tol_ppm1 = -1, tol_da2 = 0.02, tol_ppm2 = -1,
#                               predicted = "All", thread = 1,
#                               min_mz = -1, max_mz = -1, noise_threshold = 0.01, max_peak_num = -1){
#   t1 <- Sys.time()
#   if(tol_da1 != -1 & tol_ppm1 != -1) tol_da1 = -1
#   else if(tol_da1 == -1 & tol_ppm1 == -1) stop("tol is wrong!")
#
#   if(tol_da2 != -1 & tol_ppm2 != -1) tol_da2 = -1
#   else if(tol_da2 == -1 & tol_ppm2 == -1) stop("tol is wrong!")
#
#   query_number <- nrow(standardInput)
#   cmps_number <- nrow(lib$cmps)
#   ms2_number <- nrow(lib$ms2)
#   meta <- lib$meta
#   message(paste0("You are using ", "Spectral entropy", " to search library!"))
#   print(meta)
#   message(paste0("query number: ", query_number))
#   message(paste0("cmps number: ", cmps_number))
#   message(paste0("ms2 number: ", ms2_number))
#
#   loop <- function(i){
#     standardInput_tmp <- standardInput[i, ]
#     if(tol_ppm1 != -1)  tol_da1 <- standardInput_tmp$precursorMz * tol_ppm1 * 10^-6
#     lib_ms2_tmp <- lib$ms2 %>%
#       dplyr::filter(dplyr::near(precursorMz, standardInput_tmp$precursorMz, tol = tol_da1)) %>%
#       dplyr::filter(adduct == standardInput_tmp$adduct)
#     if(predicted == "All"){}
#     else if(predicted == "Yes") lib_ms2_tmp <- lib_ms2_tmp %>% dplyr::filter(predicted)
#     else if(predicted == "No") lib_ms2_tmp <- lib_ms2_tmp %>% dplyr::filter(!predicted)
#     else stop("predicted wrong!")
#
#     if(nrow(lib_ms2_tmp) == 0) return(NULL)
#
#     spMat_query <- get_spMat(standardInput_tmp)
#     spMat_query <- clean_spMat(spMat_query, tol_da2 = tol_da2, tol_ppm2 = tol_ppm2,
#                                min_mz = min_mz, max_mz = max_mz, noise_threshold = noise_threshold, max_peak_num = max_peak_num, normalize_intensity = TRUE)
#     score_vec <- sapply(1:nrow(lib_ms2_tmp), function(j) {
#       spMat_lib <- get_spMat(lib_ms2_tmp[j, ])
#       spMat_lib <- clean_spMat(spMat_lib, tol_da2 = tol_da2, tol_ppm2 = tol_ppm2,
#                                min_mz = min_mz, max_mz = max_mz, noise_threshold = noise_threshold, max_peak_num = max_peak_num, normalize_intensity = TRUE)
#       score <- msentropy::calculate_entropy_similarity(spMat_query, spMat_lib,
#                                                        ms2_tolerance_in_da = tol_da2, ms2_tolerance_in_ppm = tol_ppm2,
#                                                        clean_spectra = FALSE, min_mz = min_mz, max_mz = max_mz,
#                                                        noise_threshold = noise_threshold,
#                                                        max_peak_num = max_peak_num)
#       return(score)
#     })
#     lib_ms2_tmp$score <- score_vec
#     lib_ms2_tmp <- lib_ms2_tmp %>%
#       dplyr::filter(score >= st) %>%
#       dplyr::arrange(dplyr::desc(score))
#     return(lib_ms2_tmp)
#   }
#
#   pb <- utils::txtProgressBar(max = query_number, style = 3)
#   if(thread == 1){
#     searchRes <- lapply(1:query_number, function(i) {
#       utils::setTxtProgressBar(pb, i)
#       print(i)
#       loop(i)
#     })
#   }else if(thread > 1){
#     cl <- snow::makeCluster(thread)
#     doSNOW::registerDoSNOW(cl)
#     opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
#                                                                  n))
#     searchRes <- foreach::`%dopar%`(foreach::foreach(i = 1:query_number,
#                                                      .options.snow = opts,
#                                                      .packages = c("dplyr", "msentropy"),
#                                                      .export = c("get_spMat", "clean_spMat", "tol_da1")),
#                                     loop(i))
#     snow::stopCluster(cl)
#     gc()
#   }else stop("thread is wrong!")
#   t2 <- Sys.time()
#   stime <- sprintf("%.3f", t2 - t1)
#   message(paste0("\nTime used: ", stime, " ", attr(t2 - t1, "units"), "\n"))
#   return(searchRes)
# }

#' @title searchLib_ndotproduct
#' @description
#' The mass spectrometry database was screened using the ndotproduct algorithm.
#'
#' @param standardInput
#' A tibble, including six columns of information,
#' feature: name of the feature,
#' precursorMz: parent ion mz,
#' rt: retention time,
#' adduct: adduct form,
#' mz: mz of secondary mass spectrum,
#' intensity: intensity of secondary mass spectrum;
#' Both mz and intensity are stored in tibble as a list.
#' @param lib
#' From any of the publicMs2List databases.
#' @param joinpeak
#' See Spectra::joinPeaks type.
#' inner, left, right, outer
#' @param st
#' Threshold for scoring.
#' @param tol_da1
#' Under the Da unit, the tolerance of the same parent ion mz is searched.
#' @param tol_ppm1
#' Under the ppm unit, the tolerance of the same parent ion mz is searched.
#' @param tol_da2
#' Under the Da unit, two peaks are considered the tolerance of one peak.
#' @param tol_ppm2
#' Under the ppm unit, two peaks are considered the tolerance of one peak.
#' @param predicted
#' Whether to search only predictive mass spectra;All, Yes, No.
#' @param thread
#' The number of threads running the function.
#' @param min_mz
#' Minimum value of mz.
#' @param max_mz
#' Maximum value of mz.
#' @param noise_threshold
#' Peak intensity below noise_threshold * max(intensity) will be considered noise.
#' @param max_peak_num
#' The number of peaks in the mass spectrum.
#' @param scale_int
#' Normalized target value.
#'
#' @return
#' searchRes_ndotproduct, a list.
#'
#' @examples
#' data("standardInput", package = "MetaboSpectra")
#' load("publicMs2List.RData") # please see https://github.com/dongdongdong7/MetaboLib.ms2 and download publicMs2List.RData
#' searchRes_ndotproduct <- searchLib_ndotproduct(standardInput, lib = publicMs2List$hmdb)
# searchLib_ndotproduct <- function(standardInput, lib, joinpeak = "inner", st = 0.8,
#                                   tol_da1 = 0.01, tol_ppm1 = -1, tol_da2 = 0.02, tol_ppm2 = -1,
#                                   predicted = "All", thread = 1,
#                                   min_mz = -1, max_mz = -1, noise_threshold = 0.01, max_peak_num = -1, scale_int = 100){
#   t1 <- Sys.time()
#   query_number <- nrow(standardInput)
#   cmps_number <- nrow(lib$cmps)
#   ms2_number <- nrow(lib$ms2)
#   meta <- lib$meta
#   message(paste0("You are using ", "ndotproduct with ", joinpeak, " to search library!"))
#   print(meta)
#   message(paste0("query number: ", query_number))
#   message(paste0("cmps number: ", cmps_number))
#   message(paste0("ms2 number: ", ms2_number))
#
#   loop <- function(i){
#     standardInput_tmp <- standardInput[i, ]
#     if(tol_ppm1 != -1)  tol_da1 <- standardInput_tmp$precursorMz * tol_ppm1 * 10^-6
#     lib_ms2_tmp <- lib$ms2 %>%
#       dplyr::filter(dplyr::near(precursorMz, standardInput_tmp$precursorMz, tol = tol_da1)) %>%
#       dplyr::filter(adduct == standardInput_tmp$adduct)
#     if(predicted == "All"){}
#     else if(predicted == "Yes") lib_ms2_tmp <- lib_ms2_tmp %>% dplyr::filter(predicted)
#     else if(predicted == "No") lib_ms2_tmp <- lib_ms2_tmp %>% dplyr::filter(!predicted)
#     else stop("predicted wrong!")
#
#     if(nrow(lib_ms2_tmp) == 0) return(NULL)
#
#     spMat_query <- get_spMat(standardInput_tmp)
#     spMat_query <- clean_spMat(spMat_query, tol_da2 = tol_da2, tol_ppm2 = tol_ppm2,
#                                min_mz = min_mz, max_mz = max_mz, noise_threshold = noise_threshold, max_peak_num = max_peak_num, scale_int = scale_int)
#     score_vec <- sapply(1:nrow(lib_ms2_tmp), function(j) {
#       spMat_lib <- get_spMat(lib_ms2_tmp[j, ])
#       spMat_lib <- clean_spMat(spMat_lib, tol_da2 = tol_da2, tol_ppm2 = tol_ppm2,
#                                min_mz = min_mz, max_mz = max_mz, noise_threshold = noise_threshold, max_peak_num = max_peak_num, scale_int = scale_int)
#       score <- compare_spMat_ndotproduct(x = spMat_query, y = spMat_lib, joinpeak = joinpeak, tol_da2 = tol_da2, tol_ppm2 = tol_ppm2)
#       return(score)
#     })
#     lib_ms2_tmp$score <- score_vec
#     lib_ms2_tmp <- lib_ms2_tmp %>%
#       dplyr::filter(score >= st) %>%
#       dplyr::arrange(dplyr::desc(score))
#     return(lib_ms2_tmp)
#   }
#
#   pb <- utils::txtProgressBar(max = query_number, style = 3)
#   if(thread == 1){
#     searchRes <- lapply(1:query_number, function(i) {
#       utils::setTxtProgressBar(pb, i)
#       print(i)
#       loop(i)
#     })
#   }else if(thread > 1){
#     cl <- snow::makeCluster(thread)
#     doSNOW::registerDoSNOW(cl)
#     opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
#                                                                  n))
#     searchRes <- foreach::`%dopar%`(foreach::foreach(i = 1:query_number,
#                                                      .options.snow = opts,
#                                                      .packages = c("dplyr", "Spectra", "MsCoreUtils"),
#                                                      .export = c("get_spMat", "clean_spMat", "compare_spMat_ndotproduct", "tol_da1")),
#                                     loop(i))
#     snow::stopCluster(cl)
#     gc()
#   }else stop("thread is wrong!")
#   t2 <- Sys.time()
#   stime <- sprintf("%.3f", t2 - t1)
#   message(paste0("\nTime used: ", stime, " ", attr(t2 - t1, "units"), "\n"))
#   return(searchRes)
# }
