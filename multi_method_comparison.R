load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/publicMs2List.RData")
devtools::document()
data("standardInput")

searchRes_entropy <- searchLib_entropy(standardInput = standardInput, lib = publicMs2List$hmdb, thread = 8,
                                       st = 0.2, tol_da2 = 0.08, predicted = "All")
idenTibble_entropy <- searchRes2idenTibble(standardInput, searchRes_entropy, top = 5)
searchRes_entropy[[5]]
spMat_lib <- get_spMat(searchRes_entropy[[5]][3, ])
spMat_lib <- clean_spMat(spMat_lib)
plotSpectra(spMat_lib)
plotSpectra(spMat_lib, min_mz = 50, max_mz = 80)
spMat_query <- get_spMat(standardInput[5, ])
spMat_query <- clean_spMat(spMat_query)
plotSpectra(spMat_query)
plotSpectra(spMat_query, min_mz = 50, max_mz = 80)

low_vec <- sapply(1:nrow(standardInput), function(i) {
  standardInput_tmp <- standardInput[i, ]
  if(length(standardInput_tmp$mz[[1]]) == 2) return(TRUE)
  else return(FALSE)
})
standardInput_low <- standardInput[low_vec, ]
searchRes_low <- searchLib_entropy(standardInput = standardInput_low, lib = publicMs2List$hmdb, thread = 8,
                                       st = 0, tol_da1 = 0.01,tol_da2 = 0.02, predicted = "All")
idenTibble_low <- searchRes2idenTibble(standardInput_low, searchRes_low)
View(idenTibble_low)
i <- 29
ms2_lib <- searchRes_low[[i]] %>%
  dplyr::distinct(inchikey, .keep_all = TRUE)
ms2_lib <- searchRes_low[[i]] %>%
  dplyr::filter(accession == "HMDB0000783")
spMat_lib <- get_spMat(ms2_lib[2, ])
#spMat_lib <- get_spMat(searchRes_low[[i]][which(searchRes_low[[i]]$accession == "HMDB0001310")[2], ] )
spMat_lib <- clean_spMat(spMat_lib, noise_threshold = 0.01)
plotSpectra(spMat_lib)
plotSpectra(spMat_lib, min_mz = 70, max_mz = 140)
spMat_query <- get_spMat(idenTibble_low[i, ])
spMat_query <- clean_spMat(spMat_query, noise_threshold = 0.01)
plotSpectra(spMat_query)
plotSpectra(spMat_query, min_mz = 70, max_mz = 140)

plotComparableSpectra(spMat_query, spMat_lib, num = 20)

standardInput
idenTibble_entropy <- standardInput
idenTibble_entropy$iden_inchikey <- sapply(searchRes_entropy, function(x) {
  if(is.null(x)) return("")
  else if(nrow(x) == 0) return("")
  else{
    x <- x %>%
      dplyr::distinct(inchikey, .keep_all = TRUE)
    inchikey_vec <- x$inchikey
    return(paste0(inchikey_vec, collapse = ";"))
  }
})
idenTibble_entropy$iden_id <- sapply(searchRes_entropy, function(x) {
  if(is.null(x)) return("")
  else if(nrow(x) == 0) return("")
  else{
    x <- x %>%
      dplyr::distinct(inchikey, .keep_all = TRUE)
    accession_vec <- x$accession
    return(paste0(accession_vec, collapse = ";"))
  }
})
idenTibble_entropy$iden_score <- sapply(searchRes_entropy, function(x) {
  if(is.null(x)) return("")
  else if(nrow(x) == 0) return("")
  else{
    x <- x %>%
      dplyr::distinct(inchikey, .keep_all = TRUE)
    score_vec <- x$score
    return(paste0(score_vec, collapse = ";"))
  }
})

print(idenTibble_entropy[300:400, c("id", "iden_id", "iden_score")], n = 100)
delete_vec <- sapply(1:nrow(idenTibble_entropy), function(i) {
  if(length(idenTibble_entropy[i,]$mz[[1]]) == 0) return(FALSE)
  else return(TRUE)
})
idenTibble_entropy <- idenTibble_entropy[delete_vec, ]
iden_good1 <- sapply(1:nrow(idenTibble_entropy), function(i) {
  if(idenTibble_entropy[i, ]$id %in% strsplit(idenTibble_entropy[i, ]$iden_id, ";")[[1]][1:1]) return(1)
  else return(0)
})
sum(iden_good1) / nrow(idenTibble_entropy)

iden_good3 <- sapply(1:nrow(idenTibble_entropy), function(i) {
  if(idenTibble_entropy[i, ]$id %in% strsplit(idenTibble_entropy[i, ]$iden_id, ";")[[1]][1:3]) return(1)
  else return(0)
})
sum(iden_good3) / nrow(idenTibble_entropy)

iden_good5 <- sapply(1:nrow(idenTibble_entropy), function(i) {
  if(idenTibble_entropy[i, ]$id %in% strsplit(idenTibble_entropy[i, ]$iden_id, ";")[[1]][1:5]) return(1)
  else return(0)
})
sum(iden_good5) / nrow(idenTibble_entropy)

spMat_query <- get_spMat(standardInput[111, ])
searchRes_entropy[[111]]
spMat_lib <- get_spMat(searchRes_entropy[[111]][1, ])
spMat_query <- clean_spMat(spMat_query, noise_threshold = 0.01, normalize_intensity = TRUE, tol_da2 = 0.02)
spMat_lib <- clean_spMat(spMat_lib, noise_threshold = 0.01, normalize_intensity = TRUE, tol_da2 = 0.02)
plotSpectra(spMat_query)
plotSpectra(spMat_lib)
p <- plotComparableSpectra(spMat_query, spMat_lib, num = 20, tol_da2 = 0.02)
score <- msentropy::calculate_entropy_similarity(spMat_query, spMat_lib,
                                                 ms2_tolerance_in_da = 0.02, ms2_tolerance_in_ppm = -1,
                                                 clean_spectra = FALSE, min_mz = -1, max_mz = -1,
                                                 noise_threshold = 0.01,
                                                 max_peak_num = -1)
p <- p+ggplot2::labs(title = "HMDB0012275", subtitle = paste0("Score: ", sprintf("%.4f", score)))
p

searchRes_ndotproduct <- searchLib_ndotproduct(standardInput, lib = publicMs2List$hmdb, thread = 8,
                                               st = 0.8, tol_da2 = 0.02, joinpeak = "inner")
searchRes_ndotproduct[[1]]

{
  idenTibble_ndotproduct <- standardInput
  idenTibble_ndotproduct$iden_inchikey <- sapply(searchRes_ndotproduct, function(x) {
    if(is.null(x)) return("")
    else if(nrow(x) == 0) return("")
    else{
      x <- x %>%
        dplyr::distinct(inchikey, .keep_all = TRUE)
      inchikey_vec <- x$inchikey
      return(paste0(inchikey_vec, collapse = ";"))
    }
  })
  idenTibble_ndotproduct$iden_id <- sapply(searchRes_ndotproduct, function(x) {
    if(is.null(x)) return("")
    else if(nrow(x) == 0) return("")
    else{
      x <- x %>%
        dplyr::distinct(inchikey, .keep_all = TRUE)
      accession_vec <- x$accession
      return(paste0(accession_vec, collapse = ";"))
    }
  })
  idenTibble_ndotproduct$iden_score <- sapply(searchRes_ndotproduct, function(x) {
    if(is.null(x)) return("")
    else if(nrow(x) == 0) return("")
    else{
      x <- x %>%
        dplyr::distinct(inchikey, .keep_all = TRUE)
      score_vec <- x$score
      return(paste0(score_vec, collapse = ";"))
    }
  })

  print(idenTibble_ndotproduct[300:400, c("id", "iden_id", "iden_score")], n = 100)
  delete_vec <- sapply(1:nrow(idenTibble_ndotproduct), function(i) {
    if(length(idenTibble_ndotproduct[i,]$mz[[1]]) == 0) return(FALSE)
    else return(TRUE)
  })
  idenTibble_ndotproduct <- idenTibble_ndotproduct[delete_vec, ]
  iden_good1 <- sapply(1:nrow(idenTibble_ndotproduct), function(i) {
    if(idenTibble_ndotproduct[i, ]$id %in% strsplit(idenTibble_ndotproduct[i, ]$iden_id, ";")[[1]][1:1]) return(1)
    else return(0)
  })
  print(sum(iden_good1) / nrow(idenTibble_ndotproduct))

  iden_good3 <- sapply(1:nrow(idenTibble_ndotproduct), function(i) {
    if(idenTibble_ndotproduct[i, ]$id %in% strsplit(idenTibble_ndotproduct[i, ]$iden_id, ";")[[1]][1:3]) return(1)
    else return(0)
  })
  print(sum(iden_good3) / nrow(idenTibble_ndotproduct))

  iden_good5 <- sapply(1:nrow(idenTibble_ndotproduct), function(i) {
    if(idenTibble_ndotproduct[i, ]$id %in% strsplit(idenTibble_ndotproduct[i, ]$iden_id, ";")[[1]][1:5]) return(1)
    else return(0)
  })
  print(sum(iden_good5) / nrow(idenTibble_ndotproduct))

}

idenTibble_ndotproduct[204, ] # HMDB0031580 HMDB0002287 HMDB0001863 HMDB0000525 HMDB0001538
standardInput_tmp <- standardInput[standardInput$id == "HMDB0000525", ][1, ]
spMat_query <- get_spMat(standardInput_tmp)
spMat_query <- clean_spMat(spMat_query, normalize_intensity = TRUE)
lib_ms2_tmp <- publicMs2List$hmdb$ms2 %>%
  dplyr::filter(dplyr::near(precursorMz, standardInput_tmp$precursorMz, tol = 0.1)) %>%
  dplyr::filter(adduct == standardInput_tmp$adduct) %>%
  dplyr::filter(accession == "HMDB0033173")
spMat_lib <- get_spMat(lib_ms2_tmp[3, ])
spMat_lib <- clean_spMat(spMat_lib, normalize_intensity = TRUE)
plotSpectra(spMat_query)
plotSpectra(spMat_lib)
p <- plotComparableSpectra(spMat_query, spMat_lib, tol_da2 = 0.02, num = 20)
p <- p + ggplot2::labs(title = "HMDB0000525 - HMDB0033173", subtitle =  paste0("Score: ", sprintf("%.4f", 0.0735404), "(entropy)"))
p

compare_spMat_ndotproduct(spMat_query, spMat_lib, joinpeak = "outer", tol_da2 = 0.02)

msentropy::calculate_entropy_similarity(spMat_query, spMat_lib,
                                        ms2_tolerance_in_da = 0.02, ms2_tolerance_in_ppm = -1,
                                        clean_spectra = FALSE, min_mz = -1, max_mz = -1,
                                        noise_threshold = -1,
                                        max_peak_num = -1)
