# Amide_Pos <- Spectra::Spectra("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/Amide_Pos.msp",
#                               source = MsBackendMsp::MsBackendMsp())
# T3_Pos <- Spectra::Spectra("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/In_house/T3_Pos.msp",
#                            source = MsBackendMsp::MsBackendMsp())
# AmideInput_id_vec <- Amide_Pos$name[stringr::str_detect(Amide_Pos$name, "HMDB")] # 665
# AmideInput_precursorMz_vec <- Amide_Pos$precursorMz[stringr::str_detect(Amide_Pos$name, "HMDB")]
# AmideInput_rt_vec <- Amide_Pos$rtime[stringr::str_detect(Amide_Pos$name, "HMDB")]
# AmideInput_adduct_vec <- Amide_Pos$adduct[stringr::str_detect(Amide_Pos$name, "HMDB")]
# AmideInput_mz_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "HMDB")]), function(x) {
#   as.double(x[, "mz"])
# })
# AmideInput_intensity_list <- lapply(Spectra::peaksData(Amide_Pos[stringr::str_detect(Amide_Pos$name, "HMDB")]), function(x) {
#   as.double(x[, "intensity"])
# })
# AmideInput <- tibble::tibble(id = AmideInput_id_vec, precursorMz = AmideInput_precursorMz_vec, rt = AmideInput_rt_vec,
#                      adduct = AmideInput_adduct_vec, mz = AmideInput_mz_list, intensity = AmideInput_intensity_list)
#
#
# T3Input_id_vec <- T3_Pos$name[stringr::str_detect(T3_Pos$name, "HMDB")] # 665
# T3Input_precursorMz_vec <- T3_Pos$precursorMz[stringr::str_detect(T3_Pos$name, "HMDB")]
# T3Input_rt_vec <- T3_Pos$rtime[stringr::str_detect(T3_Pos$name, "HMDB")]
# T3Input_adduct_vec <- T3_Pos$adduct[stringr::str_detect(T3_Pos$name, "HMDB")]
# T3Input_mz_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "HMDB")]), function(x) {
#   as.double(x[, "mz"])
# })
# T3Input_intensity_list <- lapply(Spectra::peaksData(T3_Pos[stringr::str_detect(T3_Pos$name, "HMDB")]), function(x) {
#   as.double(x[, "intensity"])
# })
# T3Input <- tibble::tibble(id = T3Input_id_vec, precursorMz = T3Input_precursorMz_vec, rt = T3Input_rt_vec,
#                   adduct = T3Input_adduct_vec, mz = T3Input_mz_list, intensity = T3Input_intensity_list)
#
# standardInput <- rbind(AmideInput, T3Input) %>%
#   dplyr::arrange(precursorMz)
#
# usethis::use_data(standardInput)

