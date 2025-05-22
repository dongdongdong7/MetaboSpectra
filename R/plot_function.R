#' @title plotSpectra
#' @description
#' Plot a ggplot-style mass spectrum.
#'
#' @param spMat spectra matrix
#' @param num max peak number
#' @param min_mz The minimum value of the horizontal coordinate
#' @param max_mz The maximum value of the horizontal coordinate
#'
#' @return
#' A ggplot object
#' @export
#'
#' @examples
#' mz <- c(40.1320, 86.3400, 138.3290, 276.5710)
#' intensity <- c(1300, 4030, 10000, 31600)
#' standardRow <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat <- get_spMat(standardRow)
#' plotSpectra(spMat)
plotSpectra <- function(spMat, num = 10, min_mz = NA, max_mz = NA){
  mz <- spMat[, "mz"];intensity <- spMat[,"intensity"]
  df <- tibble::tibble(mz = mz, intensity = intensity) %>%
    dplyr::arrange(dplyr::desc(intensity))
  n <- nrow(df)
  if(num > n){
    num <- n
  }
  label_point <- c(df$intensity[1:num], rep(NA, n - num))
  label_text <- as.character(c(sprintf("%.4f", df$mz[1:num]), rep(NA, n-num)))
  df$label_point <- label_point
  df$label_text <- label_text
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mz, y = intensity)) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz, y = 0, yend = intensity)) +
    ggplot2::geom_point(ggplot2::aes(x = mz, y = label_point), color = "orange", size = 4) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label_text), size = 4, vjust = -1, min.segment.length = Inf) +
    #ggplot2::labs(title = sp2$peak_id, subtitle = paste0("Retention time: ", sprintf("%.4f", Spectra::rtime(sp2)))) +
    ggplot2::theme_light() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )
  if(!is.na(min_mz) | !is.na(max_mz)) p <- p+ggplot2::xlim(min_mz, max_mz)
  return(p)
}

#' @title plotComparableSpectra
#' @description
#' Compare the two mass spectra
#'
#' @param spMat_query Query spectra matrix
#' @param spMat_lib Library spectra matrix
#' @param num max peak number
#' @param ms2_tolerance_in_da The MS2 tolerance in Da, set to -1 to disable
#' @param ms2_tolerance_in_ppm The MS2 tolerance in ppm, set to -1 to disable
#'
#' @return
#' A ggplot object.
#' @export
#'
#' @examples
#' mz <- c(21.3303, 40.1326, 86.3402, 138.3291, 276.5717, 276.5836)
#' intensity <- c(100, 1300, 4030, 10000, 31600, 1000)
#' standardRow1 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat1 <- get_spMat(standardRow1)
#' mz <- c(21.3002, 40.1128, 86.3201, 138.3214, 276.5312)
#' intensity <- c(100, 1300, 4030, 10000, 31600)
#' standardRow2 <- tibble::tibble(mz = list(mz), intensity = list(intensity))
#' spMat2 <- get_spMat(standardRow2)
#' plotComparableSpectra(spMat1, spMat2, tol_ppm2 = 200, tol_da2 = -1)
#' plotComparableSpectra(spMat1, spMat2, tol_ppm2 = -1, tol_da2 = 0.02)
plotComparableSpectra <- function(spMat_query, spMat_lib, num = 10,
                                  ms2_tolerance_in_da = -1, ms2_tolerance_in_ppm = 5){
  if(ms2_tolerance_in_da >0 & ms2_tolerance_in_ppm > 0) ms2_tolerance_in_ppm <- -1
  if(ms2_tolerance_in_da < 0) ms2_tolerance_in_da <- 0
  if(ms2_tolerance_in_ppm < 0) ms2_tolerance_in_ppm <- 0

  mz_query <- spMat_query[, "mz"];intensity_query <- spMat_query[,"intensity"]
  df_query <- tibble::tibble(mz = mz_query, intensity = intensity_query) %>%
    dplyr::arrange(dplyr::desc(intensity))
  n_query <- nrow(df_query)
  mz_lib <- spMat_lib[, "mz"];intensity_lib <- spMat_lib[,"intensity"]
  df_lib <- tibble::tibble(mz = mz_lib, intensity = intensity_lib) %>%
    dplyr::arrange(dplyr::desc(intensity))
  n_lib <- nrow(df_lib)
  if(n_lib > num) n_lib <- num
  if(n_query > num) n_query <- num

  tmpList <- Spectra::joinPeaks(x = spMat_query, y = spMat_lib,type = "inner", tolerance = tol_da2, ppm = tol_ppm2)
  logical_x <- sapply(df_query$mz, function(x) {
    return(any(dplyr::near(x, tmpList$x[, "mz"])))
  })
  logical_y <- sapply(df_lib$mz, function(x) {
    return(any(dplyr::near(x, tmpList$y[, "mz"])))
  })
  df_query$type <- sapply(logical_x, function(x) {
    if(x) return("match")
    else return("mismatch")
  })
  df_lib$type <- sapply(logical_y, function(x) {
    if(x) return("match")
    else return("mismatch")
  })
  df_lib$intensity <- -df_lib$intensity

  df_query <- df_query[1:n_query, ]
  df_lib <- df_lib[1:n_lib, ]

  df <- rbind(df_query, df_lib)
  df$label_point <- sapply(1:nrow(df), function(i) {
    if(df[i, ]$type == "match") return(df[i, ]$intensity)
    else return(NA)
  })
  df$label_text1 <- sapply(1:nrow(df), function(i) {
    if(df[i, ]$type == "match" & df[i, ]$intensity > 0) return(round(df[i, ]$mz, 4))
    else return(NA)
  })
  df$label_text2 <- sapply(1:nrow(df), function(i) {
    if(df[i, ]$type == "match" & df[i, ]$intensity < 0) return(round(df[i, ]$mz, 4))
    else return(NA)
  })

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mz, y = intensity)) +
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz, y = 0, yend = intensity)) +
    # ggplot2::geom_point(ggplot2::aes(x = mz, y = label_point), color = "orange", size = 4) +
    # ggrepel::geom_text_repel(ggplot2::aes(label = label_text1), size = 4, vjust = 1, min.segment.length = Inf) +
    # ggrepel::geom_text_repel(ggplot2::aes(label = label_text2), size = 4, vjust = 0, min.segment.length = Inf) +
    ggplot2::theme_light() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )
  if(any(!is.na(df$label_point))){
    p <- p + ggplot2::geom_point(ggplot2::aes(x = mz, y = label_point), color = "orange", size = 4) +
      ggrepel::geom_text_repel(ggplot2::aes(label = label_text1), size = 4, vjust = 1, min.segment.length = Inf) +
      ggrepel::geom_text_repel(ggplot2::aes(label = label_text2), size = 4, vjust = 0, min.segment.length = Inf)
  }
  return(p)
}
