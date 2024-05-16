load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/inHouseMs2List.RData")
devtools::document()
standardInput_Amide_Neg <- inHouseMs2List$Amide_Neg$ms2
standardInput_Amide_Neg <- standardInput_Amide_Neg %>%
  dplyr::filter(stringr::str_detect(accession, "HMDB")) %>%
  dplyr::arrange(accession) %>%
  dplyr::select(accession, precursorMz, adduct, mz, intensity)
standardInput_Amide_Neg <- dplyr::left_join(standardInput_Amide_Neg, inHouseMs2List$Amide_Neg$cmps, by = c("accession" = "compound_id")) %>%
  dplyr::select(accession, precursorMz, rt, adduct, mz, intensity)

standardInput_Amide_Pos <- inHouseMs2List$Amide_Pos$ms2
standardInput_Amide_Pos <- standardInput_Amide_Pos %>%
  dplyr::filter(stringr::str_detect(accession, "HMDB")) %>%
  dplyr::arrange(accession) %>%
  dplyr::select(accession, precursorMz, adduct, mz, intensity)
standardInput_Amide_Pos <- dplyr::left_join(standardInput_Amide_Pos, inHouseMs2List$Amide_Pos$cmps, by = c("accession" = "compound_id")) %>%
  dplyr::select(accession, precursorMz, rt, adduct, mz, intensity)

standardInput_T3_Neg <- inHouseMs2List$T3_Neg$ms2
standardInput_T3_Neg <- standardInput_T3_Neg %>%
  dplyr::filter(stringr::str_detect(accession, "HMDB")) %>%
  dplyr::arrange(accession) %>%
  dplyr::select(accession, precursorMz, adduct, mz, intensity)
standardInput_T3_Neg <- dplyr::left_join(standardInput_T3_Neg, inHouseMs2List$T3_Neg$cmps, by = c("accession" = "compound_id")) %>%
  dplyr::select(accession, precursorMz, rt, adduct, mz, intensity)

standardInput_T3_Pos <- inHouseMs2List$T3_Pos$ms2
standardInput_T3_Pos <- standardInput_T3_Pos %>%
  dplyr::filter(stringr::str_detect(accession, "HMDB")) %>%
  dplyr::arrange(accession) %>%
  dplyr::select(accession, precursorMz, adduct, mz, intensity)
standardInput_T3_Pos <- dplyr::left_join(standardInput_T3_Pos, inHouseMs2List$T3_Pos$cmps, by = c("accession" = "compound_id")) %>%
  dplyr::select(accession, precursorMz, rt, adduct, mz, intensity)

standardInput <- rbind(standardInput_Amide_Neg, standardInput_Amide_Pos, standardInput_T3_Neg, standardInput_T3_Pos) %>%
  dplyr::arrange(accession)
colnames(standardInput)[1] <- "feature"

standardInput_identified <- standardInput
standardInput$feature <- sapply(1:nrow(standardInput), function(i) {
  tmp <- standardInput[i, ]
  rt_tmp <- round(tmp$rt, 2)
  mz_tmp <- round(tmp$precursorMz, 4)
  feature_name <- paste0(rt_tmp, "_", mz_tmp, "m/z")
})
usethis::use_data(standardInput, overwrite = TRUE)
usethis::use_data(standardInput_identified, overwrite = TRUE)
