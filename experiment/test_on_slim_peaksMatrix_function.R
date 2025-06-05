# test on slim_peaksMatrix function
# 250605
# Barry Song

set.seed(123)
spMat <- cbind(
  mz = sort(rnorm(1000, mean = 100, sd = 0.1)),
  intensity = rgamma(1000, shape = 2, scale = 1)
)
# original function
.slim_peaksMatrix = function(peaksMatrix, ppm = 5){
  mz_diff <- diff(peaksMatrix[, "mz"])
  mz_tol <- MsCoreUtils::ppm(peaksMatrix[, "mz"], ppm = ppm)
  mz_tol <- mz_tol[-length(mz_tol)]
  mz_merger_i <- which(mz_diff < mz_tol)
  if(length(mz_merger_i) != 0){
    l <- split(mz_merger_i, cumsum(c(1, diff(mz_merger_i) != 1)))
    l <- lapply(l, function(x) {c(x, x[length(x)] + 1)})
    mz_merger_i <- unlist(l)
    peakMatrix_1 <- peaksMatrix[-mz_merger_i, ]
    peakMatrix_2_list <- lapply(l, function(i) {
      matrix(data = c(median(peaksMatrix[i, "mz"]), sum(peaksMatrix[i, "intensity"])),
             ncol = 2, dimnames = list(NULL, c("mz", "intensity")))
    })
    peakMatrix_2 <- do.call("rbind", peakMatrix_2_list)
    peakMatrix_new <- rbind(peakMatrix_1, peakMatrix_2)
    peaksMatrix <- peakMatrix_new[order(peakMatrix_new[, "mz"]), ]
  }
  return(peaksMatrix)
}
system.time({
  lapply(1:10000, function(i) {
    .slim_peaksMatrix(spMat, ppm = 5)
  })
}) # 81.42 secs

system.time({
  lapply(1:10000, function(i) {
    slim_peaksMatrix(spMat, ppm = 5)
  })
}) # 2.18 secs

system.time({
  lapply(1:10000, function(i) {
    .slim_peaksMatrix_rcpp(spMat, ppm = 5)
  })
}) # 1.98 secs

spMatList <- lapply(1:10000, function(i) {spMat})
system.time({
  batch_slim_peaksMatrix(spMatList, ppm = 5)
}) # 1.99 secs

system.time({
  .batch_slim_peaksMatrix_rcpp(spMatList, ppm = 5)
}) # 1.91 secs

