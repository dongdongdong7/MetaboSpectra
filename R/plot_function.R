#' @title plotSpectra
#' @description
#' Plot a ggplot-style mass spectrum.
#'
#' @param spMat spectra matrix
#' @param num max peak number
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
plotSpectra <- function(spMat, num = 10){
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
    ggplot2::geom_segment(ggplot2::aes(x = mz, xend = mz, y = 0, yend = intensity), color = "grey") +
    ggplot2::geom_point(ggplot2::aes(x = mz, y = label_point), color = "orange", size = 4) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label_text), size = 4, vjust = -1, min.segment.length = Inf) +
    #ggplot2::labs(title = sp2$peak_id, subtitle = paste0("Retention time: ", sprintf("%.4f", Spectra::rtime(sp2)))) +
    ggplot2::theme_classic()
  return(p)
}
