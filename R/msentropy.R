#' @title Calculate spectral entropy
#'
#' @param spMat A matrix of spectrum, with two columns: m/z and intensity
#'
#' @returns A double value of spectral entropy
#' @export
#'
#' @examples
#' mz <- c(100.212, 300.321, 535.325)
#' intensity <- c(37.16, 66.83, 999.0)
#' spMat <- matrix(c(mz, intensity), ncol = 2, byrow = FALSE)
#' calculate_spectral_entropy(spMat)
calculate_spectral_entropy <- function(spMat){
  .calculate_spectral_entropy(peaks = spMat)
}

#' @title Clean a spectrum
#' @description
#' Clean a spectrum
#'
#' @param spMat A matrix of spectrum, with two columns: m/z and intensity
#' @param min_mz The minimum mz value to keep, set to -1 to disable
#' @param max_mz The maximum mz value to keep, set to -1 to disable
#' @param noise_threshold The noise threshold, set to -1 to disable, all peaks have intensity < noise_threshold * max_intensity will be removed
#' @param ppm The minimum mz difference in ppm to merge peaks, any two peaks with mz difference < ppm will be merged
#' @param max_peak_num The maximum number of peaks to keep, set to -1 to disable
#' @param normalize_intensity Whether to normalize the intensity to sum to 1
#'
#' @returns A matrix of spectral peaks, with two columns: mz and intensity
#' @export
#'
#' @examples
#' mz <- c(100.212, 169.071, 169.078, 300.321)
#' intensity <- c(0.3716, 7.917962, 100., 66.83)
#' spMat <- matrix(c(mz, intensity), ncol = 2, byrow = FALSE)
#' clean_spectrum(spMat, min_mz = 0, max_mz = 1000, noise_threshold = 0.01,
#'                ppm = -1,
#'                max_peak_num = 100, normalize_intensity = TRUE)
clean_spectrum <- function(spMat,
                           min_mz = -1, max_mz = -1,
                           noise_threshold = 0.01,
                           ppm = 5,
                           max_peak_num = -1, normalize_intensity = FALSE){
  .clean_spectrum(peaks = spMat,
                  min_mz = min_mz, max_mz = max_mz,
                  noise_threshold = noise_threshold,
                  min_ms2_difference_in_da = -1,
                  min_ms2_difference_in_ppm = ppm,
                  max_peak_num = max_peak_num,
                  normalize_intensity = normalize_intensity)
}

#' @title Unweighted entropy similarity between two spectra
#' @description Calculate the unweighted entropy similarity between two spectra.
#' Two spectra sholud be cleaned before calculating.
#'
#'
#' @param spMat1 A matrix of spectral peaks, with two columns: mz and intensity
#' @param spMat2 A matrix of spectral peaks, with two columns: mz and intensity
#' @param ppm The MS2 tolerance in ppm, set to -1 to disable
#'
#' @return The unweighted entropy similarity
#' @export
#'
#' @examples
#' mz_a <- c(169.071, 186.066, 186.0769)
#' intensity_a <- c(7.917962, 1.021589, 100.0)
#' mz_b <- c(120.212, 169.071, 186.066)
#' intensity_b <- c(37.16, 66.83, 999.0)
#' spMat1 <- matrix(c(mz_a, intensity_a), ncol = 2, byrow = FALSE)
#' spMat2 <- matrix(c(mz_b, intensity_b), ncol = 2, byrow = FALSE)
#' calculate_unweighted_entropy_similarity(spMat1, spMat2, ppm = 5)
calculate_unweighted_entropy_similarity <- function(spMat1, spMat2, ppm = 5){
  .calculate_unweighted_entropy_similarity(peaks_a = spMat1, peaks_b = spMat2,
                                           ms2_tolerance_in_da = -1, ms2_tolerance_in_ppm = ppm,
                                           clean_spectra = FALSE,
                                           min_mz = -1, max_mz = -1,
                                           noise_threshold = -1,
                                           max_peak_num = -1)
}

#' @title Entropy similarity between two spectra
#' @description Calculate the entropy similarity between two spectra
#' Two spectra sholud be cleaned before calculating.
#'
#' @param spMat1 A matrix of spectral peaks, with two columns: mz and intensity
#' @param spMat2 A matrix of spectral peaks, with two columns: mz and intensity
#' @param ppm The MS2 tolerance in ppm, set to -1 to disable
#' @return The entropy similarity
#' @export
#' @examples
#' mz_a <- c(169.071, 186.066, 186.0769)
#' intensity_a <- c(7.917962, 1.021589, 100.0)
#' mz_b <- c(120.212, 169.071, 186.066)
#' intensity_b <- c(37.16, 66.83, 999.0)
#' spMat1 <- matrix(c(mz_a, intensity_a), ncol = 2, byrow = FALSE)
#' spMat2 <- matrix(c(mz_b, intensity_b), ncol = 2, byrow = FALSE)
#' calculate_entropy_similarity(spMat1, spMat2, ppm = 5)
calculate_entropy_similarity <- function(spMat1, spMat2, ppm = 5){
  .calculate_entropy_similarity(peaks_a = spMat1, peaks_b = spMat2,
                                ms2_tolerance_in_da = -1, ms2_tolerance_in_ppm = ppm,
                                clean_spectra = FALSE,
                                min_mz = -1, max_mz = -1,
                                noise_threshold = -1,
                                max_peak_num = -1)
}
