load("D:/fudan/Projects/2024/MetaboSpectra/Progress/Database/build_package/MetaboLib.ms2/publicMs2List.RData")
devtools::document()
data("standardInput")

searchRes_entropy <- searchLib_entropy(standardInput = standardInput, lib = publicMs2List$hmdb, thread = 8, st = 0.5)
searchRes_entropy[[1]]

searchRes_ndotproduct <- searchLib_ndotproduct(standardInput, lib = publicMs2List$hmdb, thread = 8, st = 0.5)
