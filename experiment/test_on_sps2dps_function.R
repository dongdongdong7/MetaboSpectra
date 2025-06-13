demo_path <- "D:/fudan/Projects/2024/MetaboDeconv/Progress/build_package/demo_data/pos/plate1_series05.mzML"
sps <- Spectra::Spectra(
  object = demo_path,
  backend = Spectra::MsBackendMzR()
)
sps <- Spectra::filterEmptySpectra(object = sps)
sps_ms1 <- Spectra::filterMsLevel(object = sps, 1L)

dps_ms1 <- sps2dps(sps = sps_ms1, ppm = 10, noise = 100)
