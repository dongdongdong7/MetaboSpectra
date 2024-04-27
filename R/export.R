#' @title searchRes2idenTibble
#' @description
#' Change searchRes to idenTibble
#'
#' @param standardInput
#' A tibble, including six columns of information,
#' id: name of the feature,
#' precursorMz: parent ion mz,
#' rt: retention time,
#' adduct: adduct form,
#' mz: mz of secondary mass spectrum,
#' intensity: intensity of secondary mass spectrum;
#' Both mz and intensity are stored in tibble as a list.
#' @param searchRes Return of the searchLib function
#' @param top Retain the results of the top scores
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' data("standardInput", package = "MetaboSpectra")
#' load("publicMs2List.RData") # please see https://github.com/dongdongdong7/MetaboLib.ms2 and download publicMs2List.RData
#' searchRes_entropy <- searchLib_entropy(standardInput = standardInput, lib = publicMs2List$hmdb, thread = 8, st = 0.5)
#' idenTibble_entropy <- searchRes2idenTibble(standardInput, searchRes_entropy, top = 5)
searchRes2idenTibble <- function(standardInput, searchRes, top = 5){
  idenTibble <- standardInput
  idenTibble$iden_inchikey <- sapply(searchRes, function(x) {
    if(is.null(x)) return("")
    else if(nrow(x) == 0) return("")
    else{
      x <- x %>% # x should be arrange by score before!
        dplyr::distinct(inchikey, .keep_all = TRUE)
      num <- ifelse(nrow(x) > top, top, nrow(x))
      x <- x[1:num, ]
      inchikey_vec <- x$inchikey
      return(paste0(inchikey_vec, collapse = ";"))
    }
  })
  idenTibble$iden_id <- sapply(searchRes, function(x) {
    if(is.null(x)) return("")
    else if(nrow(x) == 0) return("")
    else{
      x <- x %>%
        dplyr::distinct(inchikey, .keep_all = TRUE)
      num <- ifelse(nrow(x) > top, top, nrow(x))
      x <- x[1:num, ]
      accession_vec <- x$accession
      return(paste0(accession_vec, collapse = ";"))
    }
  })
  idenTibble$iden_score <- sapply(searchRes, function(x) {
    if(is.null(x)) return("")
    else if(nrow(x) == 0) return("")
    else{
      x <- x %>%
        dplyr::distinct(inchikey, .keep_all = TRUE)
      num <- ifelse(nrow(x) > top, top, nrow(x))
      x <- x[1:num, ]
      score_vec <- x$score
      return(paste0(score_vec, collapse = ";"))
    }
  })
  return(idenTibble)
}
