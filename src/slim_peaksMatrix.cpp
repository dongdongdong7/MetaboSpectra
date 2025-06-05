#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>

Rcpp::NumericMatrix slim_peaksMatrix_rcpp_old(Rcpp::NumericMatrix peaksMatrix, double ppm = 5,
                                          std::string mz_method = "median",
                                          std::string intensity_method = "sum") {
  // peaksMatrix is a matrix with mz and intensity, and mz should be sorted
  int n = peaksMatrix.nrow();
  if (n == 0) return peaksMatrix;

  Rcpp::NumericVector mz = peaksMatrix(Rcpp::_, 0);
  Rcpp::NumericVector intensity = peaksMatrix(Rcpp::_, 1);

  // cal mz_tol by ppm
  Rcpp::NumericVector mz_tol = mz * ppm * 1e-6;

  // idx of merged peaks
  std::vector<int> mz_merger_i;
  for (int i = 0; i < n-1; ++i) {
    if ((mz[i+1] - mz[i]) < mz_tol[i]) {
      mz_merger_i.push_back(i);
    }
  }

  if (mz_merger_i.empty()) {
    return peaksMatrix;
  }

  // split + cumsum
  std::vector<std::vector<int>> l;
  if (!mz_merger_i.empty()) {
    std::vector<int> current_group = {mz_merger_i[0]};
    for (size_t i = 1; i < mz_merger_i.size(); ++i) {
      if (mz_merger_i[i] == mz_merger_i[i-1] + 1) {
        current_group.push_back(mz_merger_i[i]);
      } else {
        l.push_back(current_group);
        current_group = {mz_merger_i[i]};
      }
    }
    l.push_back(current_group);

    // change like (1,2,3) to (1,2,3,4)
    // Use & to directly modify elements in the original container
    for (auto& group : l) {
      group.push_back(group.back() + 1);
    }
  }

  // Collect all indexes that need to be deleted.
  std::vector<int> to_remove;
  for (const auto& group : l) { // Prohibiting modification of variable l
    to_remove.insert(to_remove.end(), group.begin(), group.end());
  }
  std::sort(to_remove.begin(), to_remove.end());
  to_remove.erase(std::unique(to_remove.begin(), to_remove.end()), to_remove.end()); // erase-remove

  // Create keep parts (peakMatrix_1)
  std::vector<double> keep_mz, keep_intensity;
  for (int i = 0; i < n; ++i) {
    if (std::find(to_remove.begin(), to_remove.end(), i) == to_remove.end()) { // to_remove.end() is the memory address of the next valid element in the container
      keep_mz.push_back(mz[i]);
      keep_intensity.push_back(intensity[i]);
    }
  }

  // Create merged parts (peakMatrix_2)
  std::vector<double> merged_mz, merged_intensity;
  for (const auto& group : l) {
    // 提取当前组的mz和intensity
    std::vector<double> group_mz, group_intensity;
    for (int i : group) {
      group_mz.push_back(mz[i]);
      group_intensity.push_back(intensity[i]);
    }

    // 计算合并后的mz（三种方法）
    double new_mz;
    if (mz_method == "median") {
      std::sort(group_mz.begin(), group_mz.end());
      size_t mid = group_mz.size() / 2;
      new_mz = (group_mz.size() % 2 == 0) ?
      (group_mz[mid-1] + group_mz[mid]) / 2.0 : group_mz[mid];
    } else if (mz_method == "mean") {
      new_mz = std::accumulate(group_mz.begin(), group_mz.end(), 0.0) / group_mz.size();
    } else { // max_intensity
      auto max_it = std::max_element(group_intensity.begin(), group_intensity.end());
      new_mz = group_mz[std::distance(group_intensity.begin(), max_it)];
    }

    // 计算合并后的intensity（三种方法）
    double new_intensity;
    if (intensity_method == "sum") {
      new_intensity = std::accumulate(group_intensity.begin(), group_intensity.end(), 0.0);
    } else if (intensity_method == "mean") {
      new_intensity = std::accumulate(group_intensity.begin(), group_intensity.end(), 0.0) / group_intensity.size();
    } else { // max
      new_intensity = *std::max_element(group_intensity.begin(), group_intensity.end());
    }

    merged_mz.push_back(new_mz);
    merged_intensity.push_back(new_intensity);
  }

  // Merge two part and sort
  keep_mz.insert(keep_mz.end(), merged_mz.begin(), merged_mz.end());
  keep_intensity.insert(keep_intensity.end(), merged_intensity.begin(), merged_intensity.end());

  // sort
  std::vector<size_t> idx(keep_mz.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&keep_mz](size_t i, size_t j) {
    return keep_mz[i] < keep_mz[j];
  });

  Rcpp::NumericMatrix result(keep_mz.size(), 2);
  for (size_t i = 0; i < idx.size(); ++i) {
    result(i, 0) = keep_mz[idx[i]];
    result(i, 1) = keep_intensity[idx[i]];
  }

  // Set colnames
  result.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create(),
    Rcpp::CharacterVector::create("mz", "intensity")
  );

  return result;
}

// [[Rcpp::export(name = ".slim_peaksMatrix_rcpp")]]
Rcpp::NumericMatrix slim_peaksMatrix_rcpp(Rcpp::NumericMatrix peaksMatrix, double ppm = 5,
                                           std::string mz_method = "median",
                                           std::string intensity_method = "sum") {
  // Check for empty input
  int n = peaksMatrix.nrow();
  if (n == 0) return peaksMatrix;

  // Direct access to matrix columns
  Rcpp::NumericVector mz = peaksMatrix(Rcpp::_, 0);
  Rcpp::NumericVector intensity = peaksMatrix(Rcpp::_, 1);

  // Pre-calculate mz tolerances
  Rcpp::NumericVector mz_tol = mz * ppm * 1e-6;

  // Identify peaks to merge
  std::vector<int> merge_starts;
  std::vector<int> merge_ends;

  int start = 0;
  for (int i = 0; i < n - 1; ++i) {
    if ((mz[i+1] - mz[i]) >= mz_tol[i]) {
      if (i > start) {  // Found a merge group
        merge_starts.push_back(start);
        merge_ends.push_back(i);
      }
      start = i + 1;
    }
  }

  // Check for last group
  if (n - 1 > start) {
    merge_starts.push_back(start);
    merge_ends.push_back(n - 1);
  }

  // If no merges needed, return original matrix
  if (merge_starts.empty()) {
    return peaksMatrix;
  }

  // Prepare output vectors
  std::vector<double> result_mz;
  std::vector<double> result_intensity;
  result_mz.reserve(n);  // Reserve space to avoid reallocations
  result_intensity.reserve(n);

  // Process non-merged peaks before first merge group
  for (int i = 0; i < merge_starts[0]; ++i) {
    result_mz.push_back(mz[i]);
    result_intensity.push_back(intensity[i]);
  }

  // Process merge groups
  for (size_t g = 0; g < merge_starts.size(); ++g) {
    int start_idx = merge_starts[g];
    int end_idx = merge_ends[g];
    int group_size = end_idx - start_idx + 1;

    // Extract group data
    Rcpp::NumericVector group_mz = mz[Rcpp::Range(start_idx, end_idx)];
    Rcpp::NumericVector group_intensity = intensity[Rcpp::Range(start_idx, end_idx)];

    // Calculate merged mz
    double new_mz;
    if (mz_method == "median") {
      Rcpp::NumericVector sorted_mz = Rcpp::clone(group_mz);
      std::sort(sorted_mz.begin(), sorted_mz.end());
      new_mz = sorted_mz[group_size / 2];
      if (group_size % 2 == 0) {
        new_mz = (new_mz + sorted_mz[group_size / 2 - 1]) / 2.0;
      }
    } else if (mz_method == "mean") {
      new_mz = Rcpp::mean(group_mz);
    } else { // max_intensity
      int max_pos = Rcpp::which_max(group_intensity);
      new_mz = group_mz[max_pos];
    }

    // Calculate merged intensity
    double new_intensity;
    if (intensity_method == "sum") {
      new_intensity = Rcpp::sum(group_intensity);
    } else if (intensity_method == "mean") {
      new_intensity = Rcpp::mean(group_intensity);
    } else { // max
      new_intensity = Rcpp::max(group_intensity);
    }

    result_mz.push_back(new_mz);
    result_intensity.push_back(new_intensity);

    // Add non-merged peaks between current and next merge group
    if (g < merge_starts.size() - 1) {
      for (int i = end_idx + 1; i < merge_starts[g+1]; ++i) {
        result_mz.push_back(mz[i]);
        result_intensity.push_back(intensity[i]);
      }
    }
  }

  // Process non-merged peaks after last merge group
  for (int i = merge_ends.back() + 1; i < n; ++i) {
    result_mz.push_back(mz[i]);
    result_intensity.push_back(intensity[i]);
  }

  // Create result matrix
  Rcpp::NumericMatrix result(result_mz.size(), 2);
  std::copy(result_mz.begin(), result_mz.end(), result.column(0).begin());
  std::copy(result_intensity.begin(), result_intensity.end(), result.column(1).begin());

  // Set column names
  result.attr("dimnames") = Rcpp::List::create(
    Rcpp::CharacterVector::create(),
    Rcpp::CharacterVector::create("mz", "intensity")
  );

  return result;
}

// [[Rcpp::export(name = ".batch_slim_peaksMatrix_rcpp")]]
Rcpp::List batch_slim_peaksMatrix_rcpp(Rcpp::List peaksMatrixList, double ppm = 5,
                                        std::string mz_method = "median",
                                        std::string intensity_method = "sum"){
  int n = peaksMatrixList.size();
  Rcpp::List out(n);

  Rcpp::List::iterator it;
  int i = 0;
  for(it = peaksMatrixList.begin(); it != peaksMatrixList.end(); ++it, ++i){
    out[i] = slim_peaksMatrix_rcpp(*it, ppm, mz_method, intensity_method);
  }

  return out;
}
