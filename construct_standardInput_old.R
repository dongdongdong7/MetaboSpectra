load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/publicMs2List.RData")
devtools::document()
Amide_Pos <- Spectra::Spectra("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/Amide_Pos.msp",
                              source = MsBackendMsp::MsBackendMsp())
T3_Pos <- Spectra::Spectra("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/T3_Pos.msp",
                           source = MsBackendMsp::MsBackendMsp())
# HMDB
{
  AmideInput_id_vec <- Amide_Pos$name[stringr::str_detect(Amide_Pos$name, "HMDB")] # 689
  AmideInput_precursorMz_vec <- Amide_Pos$precursorMz[stringr::str_detect(Amide_Pos$name, "HMDB")]
  AmideInput_rt_vec <- Amide_Pos$rtime[stringr::str_detect(Amide_Pos$name, "HMDB")]
  AmideInput_adduct_vec <- Amide_Pos$adduct[stringr::str_detect(Amide_Pos$name, "HMDB")]
  AmideInput_mz_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "HMDB")]), function(x) {
    as.double(x[, "mz"])
  })
  AmideInput_intensity_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "HMDB")]), function(x) {
    as.double(x[, "intensity"])
  })
  AmideInput <- tibble::tibble(id = AmideInput_id_vec, precursorMz = AmideInput_precursorMz_vec, rt = AmideInput_rt_vec,
                               adduct = AmideInput_adduct_vec, mz = AmideInput_mz_list, intensity = AmideInput_intensity_list)


  T3Input_id_vec <- T3_Pos$name[stringr::str_detect(T3_Pos$name, "HMDB")] # 959
  T3Input_precursorMz_vec <- T3_Pos$precursorMz[stringr::str_detect(T3_Pos$name, "HMDB")]
  T3Input_rt_vec <- T3_Pos$rtime[stringr::str_detect(T3_Pos$name, "HMDB")]
  T3Input_adduct_vec <- T3_Pos$adduct[stringr::str_detect(T3_Pos$name, "HMDB")]
  T3Input_mz_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "HMDB")]), function(x) {
    as.double(x[, "mz"])
  })
  T3Input_intensity_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "HMDB")]), function(x) {
    as.double(x[, "intensity"])
  })
  T3Input <- tibble::tibble(id = T3Input_id_vec, precursorMz = T3Input_precursorMz_vec, rt = T3Input_rt_vec,
                            adduct = T3Input_adduct_vec, mz = T3Input_mz_list, intensity = T3Input_intensity_list)

  AmideInput <- AmideInput %>%
    dplyr::arrange(id)
  T3Input <- T3Input %>%
    dplyr::arrange(id)

  T3Input_single <- T3Input[!(duplicated(T3Input$id, fromLast = TRUE) | duplicated(T3Input$id, fromLast = FALSE)),]
  T3Input_multi <- T3Input[(duplicated(T3Input$id, fromLast = TRUE) | duplicated(T3Input$id, fromLast = FALSE)),]
  searchRes_entropy <- searchLib_entropy(standardInput = T3Input_multi, lib = publicMs2List$hmdb, thread = 8,
                                         st = 0.2, tol_da2 = 0.08, predicted = "All")
  idenTibble_entropy <- searchRes2idenTibble(T3Input_multi, searchRes_entropy, top = 5)
  openxlsx::write.xlsx(idenTibble_entropy, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/idenTibble_T3.xlsx")

  T3Input_new <- T3Input_multi[0, ]
  for(id in unique(idenTibble_entropy$id)){
    T3Input_tmp <- T3Input_multi[idenTibble_entropy$id == id, ]
    tmp_tibble <- idenTibble_entropy[idenTibble_entropy$id == id, ]
    score_vec <- sapply(1:nrow(tmp_tibble), function(i) {
      tmp <- tmp_tibble[i, ]
      idx <- which(strsplit(tmp$iden_id, ";")[[1]] == id)
      if(length(idx) == 0) return(0)
      score <- as.numeric(strsplit(tmp$iden_score, ";")[[1]][idx])
      return(score)
    })
    max_idx <- which.max(score_vec)
    T3Input_new <- rbind(T3Input_new, T3Input_tmp[max_idx, ])
  }
  T3Input <- rbind(T3Input_single, T3Input_new)

  spMat_query <- get_spMat(T3Input_multi[2, ])
  spMat_query <- clean_spMat(spMat_query)
  plotSpectra(spMat_query)
  spMat_lib <- get_spMat(searchRes_entropy[[135]][1, ])
  spMat_lib <- clean_spMat(spMat_lib)
  plotSpectra(spMat_lib)
  plotComparableSpectra(spMat_query, spMat_lib, num = 20, tol_da2 = 0.08) + ggplot2::xlim(c(60, 180))

  AmideInput_single <- AmideInput[!(duplicated(AmideInput$id, fromLast = TRUE) | duplicated(AmideInput$id, fromLast = FALSE)), ]
  AmideInput_multi <- AmideInput[(duplicated(AmideInput$id, fromLast = TRUE) | duplicated(AmideInput$id, fromLast = FALSE)), ]
  searchRes_entropy <- searchLib_entropy(standardInput = AmideInput_multi, lib = publicMs2List$hmdb, thread = 8,
                                         st = 0.2, tol_da2 = 0.08, predicted = "All")
  idenTibble_entropy <- searchRes2idenTibble(AmideInput_multi, searchRes_entropy, top = 5)
  openxlsx::write.xlsx(idenTibble_entropy, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/idenTibble_Amide.xlsx")

  AmideInput_new <- AmideInput_multi[0, ]
  for(id in unique(idenTibble_entropy$id)){
    AmideInput_tmp <- AmideInput_multi[idenTibble_entropy$id == id, ]
    tmp_tibble <- idenTibble_entropy[idenTibble_entropy$id == id, ]
    score_vec <- sapply(1:nrow(tmp_tibble), function(i) {
      tmp <- tmp_tibble[i, ]
      idx <- which(strsplit(tmp$iden_id, ";")[[1]] == id)
      if(length(idx) == 0) return(0)
      score <- as.numeric(strsplit(tmp$iden_score, ";")[[1]][idx])
      return(score)
    })
    max_idx <- which.max(score_vec)
    AmideInput_new <- rbind(AmideInput_new, AmideInput_tmp[max_idx, ])
  }
  AmideInput <- rbind(AmideInput_single, AmideInput_new)

  # 有一些是二级hmdb_id, 手动更改
  AmideInput[AmideInput$id == "HMDB0015157", ]$id <- "HMDB0003555"
  # 这两个直接改了msp, 因为改后有重复
  # AmideInput[AmideInput$id == "HMDB0006055", ]$id <- "HMDB0000725"
  # AmideInput[AmideInput$id == "HMDB0000846", ]$id <- "HMDB0000222"
  # AmideInput[AmideInput$id == "HMDB00870", ]$id <- "HMDB0000870"

  data("metabolitesList", package = "MetaboRich")
  Amide_cmps <- AmideInput %>%
    dplyr::select(id, rt)
  Amide_cmps <- dplyr::left_join(Amide_cmps, metabolitesList$hmdb, by = c("id" = "hmdb_id")) %>%
    dplyr::select(id, name, inchikey, formula, monoisotop_mass, synonyms, rt)
  Amide_cmps[which(is.na(Amide_cmps$name)), ]
  idx <- which(is.na(Amide_cmps$name))
  Amide_cmps[idx, ]$name <- c("5,7-dihydroxy-2-(4-hydroxy-3,5-dimethoxyphenyl)-4H-chromen-4-one", "2-(3,4-dihydroxyphenyl)-3,4-dihydro-2H-1-benzopyran-3,5,7-triol", "6-hydroxy-2-phenyl-4H-chromen-4-one", "D-Glucosaminic acid")
  Amide_cmps[idx, ]$inchikey <- c(NA, NA, "GPZYYYGYCRFPBU-UHFFFAOYSA-N", "UFYKDFXCZBTLOO-TXICZTDVSA-N")
  Amide_cmps[idx, ]$formula <- c("C17H14O7", "C15H14O6", "C15H10O3", "C6H13NO6")
  Amide_cmps[idx, ]$monoisotop_mass <- c(330.074, 290.079, 238.0630, 195.0743)
  Amide_ms2 <- AmideInput %>%
    dplyr::select(id, precursorMz, adduct, mz, intensity)
  Amide_ms2 <- dplyr::left_join(Amide_ms2, Amide_cmps, by = c("id" = "id")) %>%
    dplyr::select(inchikey, id, adduct, precursorMz, mz, intensity)
  Amide_ms2$collision_energy <- NA
  Amide_ms2$instrument_type <- NA
  Amide_ms2$instrument <- NA
  Amide_ms2$predicted <- FALSE
  colnames(Amide_ms2)[2] <- "accession"
  Amide_ms2 <- Amide_ms2 %>%
    dplyr::select(inchikey, accession, adduct, collision_energy, instrument_type, instrument, precursorMz, predicted, mz, intensity)
  # 手动修改
  T3Input[T3Input$id == "HMDB0014330", ]$id <- "HMDB0001934"
  # 这两个直接改了msp, 因为改后有重复
  # T3Input[T3Input$id == "HMDB0000805", ]$id <- "HMDB0000267"
  # T3Input[T3Input$id == "HMDB0000846", ]$id <- "HMDB0000222"
  # T3Input[T3Input$id == "HMDB0006055", ]$id <- "HMDB0000725"
  # T3Input[T3Input$id == "HMDB00870", ]$id <- "HMDB0000870"
  # 没有hmdb一级记载: HMDB0124861 HMDB0127721, HMDB0141530
  T3_cmps <- T3Input %>%
    dplyr::select(id, rt)
  T3_cmps <- dplyr::left_join(T3_cmps, metabolitesList$hmdb, by = c("id" = "hmdb_id")) %>%
    dplyr::select(id, name, inchikey, formula, monoisotop_mass, synonyms, rt)
  T3_cmps[which(is.na(T3_cmps$name)), ]
  idx <- which(is.na(T3_cmps$name))
  T3_cmps[idx, ]
  T3_cmps[idx, ]$name <- c("5,7-dihydroxy-2-(4-hydroxy-3,5-dimethoxyphenyl)-4H-chromen-4-one",
                           "2-(3,4-dihydroxyphenyl)-3,4-dihydro-2H-1-benzopyran-3,5,7-triol",
                           "6-hydroxy-2-phenyl-4H-chromen-4-one")
  T3_cmps[idx, ]$inchikey <- c(NA, NA, "GPZYYYGYCRFPBU-UHFFFAOYSA-N")
  T3_cmps[idx, ]$formula <- c("C17H14O7", "C15H14O6", "C15H10O3")
  T3_cmps[idx, ]$monoisotop_mass <- c(330.074, 290.079, 238.0630)
  T3_ms2 <- T3Input %>%
    dplyr::select(id, precursorMz, adduct, mz, intensity)
  T3_ms2 <- dplyr::left_join(T3_ms2, T3_cmps, by = c("id" = "id")) %>%
    dplyr::select(inchikey, id, adduct, precursorMz, mz, intensity)
  T3_ms2$collision_energy <- NA
  T3_ms2$instrument_type <- NA
  T3_ms2$instrument <- NA
  T3_ms2$predicted <- FALSE
  colnames(T3_ms2)[2] <- "accession"
  T3_ms2 <- T3_ms2 %>%
    dplyr::select(inchikey, accession, adduct, collision_energy, instrument_type, instrument, precursorMz, predicted, mz, intensity)

  Amide_cmps_hmdb <- Amide_cmps
  Amide_ms2_hmdb <- Amide_ms2
  T3_cmps_hmdb <- T3_cmps
  T3_ms2_hmdb <- T3_ms2
}
# CID
{
  AmideInput_id_vec <- Amide_Pos$name[stringr::str_detect(Amide_Pos$name, "CID")] # 62
  AmideInput_precursorMz_vec <- Amide_Pos$precursorMz[stringr::str_detect(Amide_Pos$name, "CID")]
  AmideInput_rt_vec <- Amide_Pos$rtime[stringr::str_detect(Amide_Pos$name, "CID")]
  AmideInput_adduct_vec <- Amide_Pos$adduct[stringr::str_detect(Amide_Pos$name, "CID")]
  AmideInput_mz_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "CID")]), function(x) {
    as.double(x[, "mz"])
  })
  AmideInput_intensity_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "CID")]), function(x) {
    as.double(x[, "intensity"])
  })
  AmideInput <- tibble::tibble(id = AmideInput_id_vec, precursorMz = AmideInput_precursorMz_vec, rt = AmideInput_rt_vec,
                               adduct = AmideInput_adduct_vec, mz = AmideInput_mz_list, intensity = AmideInput_intensity_list)

  T3Input_id_vec <- T3_Pos$name[stringr::str_detect(T3_Pos$name, "CID")] # 99
  T3Input_precursorMz_vec <- T3_Pos$precursorMz[stringr::str_detect(T3_Pos$name, "CID")]
  T3Input_rt_vec <- T3_Pos$rtime[stringr::str_detect(T3_Pos$name, "CID")]
  T3Input_adduct_vec <- T3_Pos$adduct[stringr::str_detect(T3_Pos$name, "CID")]
  T3Input_mz_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "CID")]), function(x) {
    as.double(x[, "mz"])
  })
  T3Input_intensity_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "CID")]), function(x) {
    as.double(x[, "intensity"])
  })
  T3Input <- tibble::tibble(id = T3Input_id_vec, precursorMz = T3Input_precursorMz_vec, rt = T3Input_rt_vec,
                            adduct = T3Input_adduct_vec, mz = T3Input_mz_list, intensity = T3Input_intensity_list)

  AmideInput <- AmideInput %>%
    dplyr::arrange(id)
  T3Input <- T3Input %>%
    dplyr::arrange(id)

  Amide_cmps <- AmideInput %>%
    dplyr::select(id, rt)
  openxlsx::write.xlsx(Amide_cmps, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/Amide_cmps_CID.xlsx")
  T3_cmps <- T3Input %>%
    dplyr::select(id, rt)
  openxlsx::write.xlsx(T3_cmps, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/T3_cmps_CID.xlsx")

  Amide_cmps <- dplyr::as_tibble(openxlsx::read.xlsx("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/Amide_cmps_CID_records.xlsx", sheet = 1)) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  T3_cmps <- dplyr::as_tibble(openxlsx::read.xlsx("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/T3_cmps_CID_records.xlsx", sheet = 1)) %>%
    dplyr::distinct(id, .keep_all = TRUE)

  Amide_ms2 <- AmideInput %>%
    dplyr::select(id, precursorMz, adduct, mz, intensity)
  Amide_ms2 <- dplyr::left_join(Amide_ms2, Amide_cmps, by = c("id" = "id")) %>%
    dplyr::select(inchikey, id, adduct, precursorMz, mz, intensity)
  Amide_ms2$collision_energy <- NA
  Amide_ms2$instrument_type <- NA
  Amide_ms2$instrument <- NA
  Amide_ms2$predicted <- FALSE
  colnames(Amide_ms2)[2] <- "accession"
  Amide_ms2 <- Amide_ms2 %>%
    dplyr::select(inchikey, accession, adduct, collision_energy, instrument_type, instrument, precursorMz, predicted, mz, intensity)

  T3_ms2 <- T3Input %>%
    dplyr::select(id, precursorMz, adduct, mz, intensity)
  T3_ms2 <- dplyr::left_join(T3_ms2, T3_cmps, by = c("id" = "id")) %>%
    dplyr::select(inchikey, id, adduct, precursorMz, mz, intensity)
  T3_ms2$collision_energy <- NA
  T3_ms2$instrument_type <- NA
  T3_ms2$instrument <- NA
  T3_ms2$predicted <- FALSE
  colnames(T3_ms2)[2] <- "accession"
  T3_ms2 <- T3_ms2 %>%
    dplyr::select(inchikey, accession, adduct, collision_energy, instrument_type, instrument, precursorMz, predicted, mz, intensity)

  Amide_cmps_cid <- Amide_cmps
  Amide_ms2_cid <- Amide_ms2
  T3_cmps_cid <- T3_cmps
  T3_ms2_cid <- T3_ms2
}
# CAS
{
  AmideInput_id_vec <- Amide_Pos$name[stringr::str_detect(Amide_Pos$name, "CAS")] # 7
  AmideInput_precursorMz_vec <- Amide_Pos$precursorMz[stringr::str_detect(Amide_Pos$name, "CAS")]
  AmideInput_rt_vec <- Amide_Pos$rtime[stringr::str_detect(Amide_Pos$name, "CAS")]
  AmideInput_adduct_vec <- Amide_Pos$adduct[stringr::str_detect(Amide_Pos$name, "CAS")]
  AmideInput_mz_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "CAS")]), function(x) {
    as.double(x[, "mz"])
  })
  AmideInput_intensity_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "CAS")]), function(x) {
    as.double(x[, "intensity"])
  })
  AmideInput <- tibble::tibble(id = AmideInput_id_vec, precursorMz = AmideInput_precursorMz_vec, rt = AmideInput_rt_vec,
                               adduct = AmideInput_adduct_vec, mz = AmideInput_mz_list, intensity = AmideInput_intensity_list)

  T3Input_id_vec <- T3_Pos$name[stringr::str_detect(T3_Pos$name, "CAS")] # 14
  T3Input_precursorMz_vec <- T3_Pos$precursorMz[stringr::str_detect(T3_Pos$name, "CAS")]
  T3Input_rt_vec <- T3_Pos$rtime[stringr::str_detect(T3_Pos$name, "CAS")]
  T3Input_adduct_vec <- T3_Pos$adduct[stringr::str_detect(T3_Pos$name, "CAS")]
  T3Input_mz_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "CAS")]), function(x) {
    as.double(x[, "mz"])
  })
  T3Input_intensity_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "CAS")]), function(x) {
    as.double(x[, "intensity"])
  })
  T3Input <- tibble::tibble(id = T3Input_id_vec, precursorMz = T3Input_precursorMz_vec, rt = T3Input_rt_vec,
                            adduct = T3Input_adduct_vec, mz = T3Input_mz_list, intensity = T3Input_intensity_list)

  AmideInput <- AmideInput %>%
    dplyr::arrange(id)
  T3Input <- T3Input %>%
    dplyr::arrange(id)

  Amide_cmps <- AmideInput %>%
    dplyr::select(id, rt)
  openxlsx::write.xlsx(Amide_cmps, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/Amide_cmps_CAS.xlsx")
  T3_cmps <- T3Input %>%
    dplyr::select(id, rt)
  openxlsx::write.xlsx(T3_cmps, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/T3_cmps_CAS.xlsx")

  Amide_cmps <- dplyr::as_tibble(openxlsx::read.xlsx("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/Amide_cmps_CAS_records.xlsx", sheet = 1)) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  T3_cmps <- dplyr::as_tibble(openxlsx::read.xlsx("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/T3_cmps_CAS_records.xlsx", sheet = 1)) %>%
    dplyr::distinct(id, .keep_all = TRUE)

  Amide_ms2 <- AmideInput %>%
    dplyr::select(id, precursorMz, adduct, mz, intensity)
  Amide_ms2 <- dplyr::left_join(Amide_ms2, Amide_cmps, by = c("id" = "id")) %>%
    dplyr::select(inchikey, id, adduct, precursorMz, mz, intensity)
  Amide_ms2$collision_energy <- NA
  Amide_ms2$instrument_type <- NA
  Amide_ms2$instrument <- NA
  Amide_ms2$predicted <- FALSE
  colnames(Amide_ms2)[2] <- "accession"
  Amide_ms2 <- Amide_ms2 %>%
    dplyr::select(inchikey, accession, adduct, collision_energy, instrument_type, instrument, precursorMz, predicted, mz, intensity)

  T3_ms2 <- T3Input %>%
    dplyr::select(id, precursorMz, adduct, mz, intensity)
  T3_ms2 <- dplyr::left_join(T3_ms2, T3_cmps, by = c("id" = "id")) %>%
    dplyr::select(inchikey, id, adduct, precursorMz, mz, intensity)
  T3_ms2$collision_energy <- NA
  T3_ms2$instrument_type <- NA
  T3_ms2$instrument <- NA
  T3_ms2$predicted <- FALSE
  colnames(T3_ms2)[2] <- "accession"
  T3_ms2 <- T3_ms2 %>%
    dplyr::select(inchikey, accession, adduct, collision_energy, instrument_type, instrument, precursorMz, predicted, mz, intensity)

  Amide_cmps_cas <- Amide_cmps
  Amide_ms2_cas <- Amide_ms2
  T3_cmps_cas <- T3_cmps
  T3_ms2_cas <- T3_ms2
}
# YMDB Find nothing
{
  AmideInput_id_vec <- Amide_Pos$name[stringr::str_detect(Amide_Pos$name, "YMDB")] # 2
  AmideInput_precursorMz_vec <- Amide_Pos$precursorMz[stringr::str_detect(Amide_Pos$name, "YMDB")]
  AmideInput_rt_vec <- Amide_Pos$rtime[stringr::str_detect(Amide_Pos$name, "YMDB")]
  AmideInput_adduct_vec <- Amide_Pos$adduct[stringr::str_detect(Amide_Pos$name, "YMDB")]
  AmideInput_mz_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "YMDB")]), function(x) {
    as.double(x[, "mz"])
  })
  AmideInput_intensity_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "YMDB")]), function(x) {
    as.double(x[, "intensity"])
  })
  AmideInput <- tibble::tibble(id = AmideInput_id_vec, precursorMz = AmideInput_precursorMz_vec, rt = AmideInput_rt_vec,
                               adduct = AmideInput_adduct_vec, mz = AmideInput_mz_list, intensity = AmideInput_intensity_list)

  T3Input_id_vec <- T3_Pos$name[stringr::str_detect(T3_Pos$name, "YMDB")] # 4
  T3Input_precursorMz_vec <- T3_Pos$precursorMz[stringr::str_detect(T3_Pos$name, "YMDB")]
  T3Input_rt_vec <- T3_Pos$rtime[stringr::str_detect(T3_Pos$name, "YMDB")]
  T3Input_adduct_vec <- T3_Pos$adduct[stringr::str_detect(T3_Pos$name, "YMDB")]
  T3Input_mz_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "YMDB")]), function(x) {
    as.double(x[, "mz"])
  })
  T3Input_intensity_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "YMDB")]), function(x) {
    as.double(x[, "intensity"])
  })
  T3Input <- tibble::tibble(id = T3Input_id_vec, precursorMz = T3Input_precursorMz_vec, rt = T3Input_rt_vec,
                            adduct = T3Input_adduct_vec, mz = T3Input_mz_list, intensity = T3Input_intensity_list)

  AmideInput <- AmideInput %>%
    dplyr::arrange(id)
  T3Input <- T3Input %>%
    dplyr::arrange(id)

  Amide_cmps <- AmideInput %>%
    dplyr::select(id, rt)
  openxlsx::write.xlsx(Amide_cmps, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/Amide_cmps_YMDB.xlsx")
  T3_cmps <- T3Input %>%
    dplyr::select(id, rt)
  openxlsx::write.xlsx(T3_cmps, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/T3_cmps_YMDB.xlsx")
}

Amide_meta <- data.frame(name = c("source", "url", "source_version", "source_date", "organism"),
                         value = c("TangLab", NA, NA, "2024-05", NA))
Amide_cmps_cid$synonyms <- lapply(Amide_cmps_cid$id, function(x) {return(NULL)})
Amide_cmps_cid <- Amide_cmps_cid %>%
  dplyr::select(id, name, inchikey, formula, monoisotop_mass, synonyms, rt)
Amide_cmps_cas$synonyms <- lapply(Amide_cmps_cas$id, function(x) {return(NULL)})
Amide_cmps_cas <- Amide_cmps_cas %>%
  dplyr::select(id, name, inchikey, formula, monoisotop_mass, synonyms, rt)
Amide_cmps <- rbind(Amide_cmps_hmdb, Amide_cmps_cid, Amide_cmps_cas)
Amide_ms2 <- rbind(Amide_ms2_hmdb, Amide_ms2_cid, Amide_ms2_cas)

inHouseMs2List <- list()
inHouseMs2List$Amide <- list(cmps = Amide_cmps, ms2 = Amide_ms2, meta = Amide_meta)
T3_meta <- data.frame(name = c("source", "url", "source_version", "source_date", "organism"),
                         value = c("TangLab", NA, NA, "2024-05", NA))
inHouseMs2List$T3 <- list(cmps = T3_cmps, ms2 = T3_ms2, meta = T3_meta)
save(inHouseMs2List, file = "D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/240505/inHouseMs2List.RData")

standardInput <- rbind(AmideInput, T3Input) %>%
  dplyr::arrange(precursorMz)
remain_vec <- sapply(1:nrow(standardInput), function(i) {
  standardInput_tmp <- standardInput[i, ]
  if(length(standardInput_tmp$mz[[1]]) >= 2) return(TRUE)
  else return(FALSE)
})
standardInput <- standardInput[remain_vec, ]
standardInput <- standardInput %>%
  dplyr::arrange(id)

usethis::use_data(standardInput, overwrite = TRUE)

