#ifndef CENTRALIYHELPER_H_
#define CENTRALIYHELPER_H_

#include "TH1.h"

#include <iostream>

namespace cent
{

enum GlauberQuantity {
  kNpart = 0,
  kNcoll,
  kNanc,
  kEcc,
  kB
};

// Returns the value X such that lPercentileRequested% of entries lie above X.
inline double getBoundaryForPercentile(TH1* histo, double lPercentileRequested)
{
  if (!histo || histo->GetEntries() <= 0) {
    return 0.0;
  }

  double lPercentile = 100.0 - lPercentileRequested;
  const long lNBins = histo->GetNbinsX();
  double lTotal = histo->Integral();
  double lCountDesired = lPercentile * lTotal / 100.0;
  double lCount = 0.0;

  for (long ibin = 1; ibin <= lNBins; ibin++) {
    double lPrevCount = lCount;
    lCount += histo->GetBinContent(ibin);
    if (lCount >= lCountDesired) {
      double lWidth = histo->GetBinWidth(ibin);
      double lLeftPercentile  = 100.0 * lPrevCount / lTotal;
      double lRightPercentile = 100.0 * lCount / lTotal;
      double lProportion = (lRightPercentile > lLeftPercentile) ? (lPercentile - lLeftPercentile) / (lRightPercentile - lLeftPercentile) : 0.0;

      return histo->GetBinLowEdge(ibin) + lProportion * lWidth;
    }
  }
  return histo->GetBinLowEdge(lNBins + 1); // fallback: top edge of last bin
}

template <typename T>
inline T sumInQuadrature(const std::vector<T>& v)
{
  T sumSq{};
  for (const T& x : v) {
    sumSq += x * x;
  }
  return std::sqrt(sumSq);
}

template <typename T>
inline T weightedAverage(const std::vector<T>& vals, const std::vector<T>& weights)
{
  if (vals.size() != weights.size() || vals.empty()) {
    return T{};
  }

  T weightedSum{};
  T totalWeight{};
  for (size_t ii = 0; ii < vals.size(); ++ii) {
    weightedSum += vals[ii] * weights[ii];
    totalWeight += weights[ii];
  }

  if (totalWeight == T{}) {
    return T{};
  }

  return weightedSum / totalWeight;
}

template <typename T>
inline std::vector<T> calculateDifference(const std::vector<T>& vals, const std::vector<T>& refs)
{
  if (vals.size() != refs.size()) {
    std::cerr << "Error in calculateDifference: vals(" << vals.size()
               << ") and refs(" << refs.size()
               << ") do not have the same size!" << std::endl;
    return {};
  }

  std::vector<T> diff;
  diff.reserve(vals.size());
  for (size_t i = 0; i < vals.size(); ++i) {
    diff.push_back(std::abs(vals[i] - refs[i]));
  }
  return diff;
}

template <typename T>
inline std::vector<T> calculateRelativeDifference(const std::vector<T>& vals, const std::vector<T>& refs)
{
  if (vals.size() != refs.size()) {
    std::cerr << "Error in calculateRelativeDifference: vals(" << vals.size()
               << ") and refs(" << refs.size()
               << ") do not have the same size!" << std::endl;
    return {};
  }

  std::vector<T> relDiff;
  relDiff.reserve(vals.size());
  for (size_t i = 0; i < vals.size(); ++i) {
    if (refs[i] == T{}) {
      std::cerr << "Warning in calculateRelativeDifference: refs[" << i
                 << "] is zero, pushing 0 to avoid division by zero." << std::endl;
      relDiff.push_back(T{});
      continue;
    }
    relDiff.push_back(std::abs(vals[i] - refs[i]) / refs[i]);
  }
  return relDiff;
}

template <typename T>
inline std::vector<T> averageOverVariations(const std::vector<std::vector<T>>& variations)
{
  if (variations.empty()) {
    return {};
  }

  const size_t nBins = variations[0].size();
  std::vector<T> result(nBins, T{});

  for (const auto& variation : variations) {
    for (size_t ibin = 0; ibin < nBins; ++ibin) {
      result[ibin] += variation[ibin];
    }
  }

  for (auto& val : result) {
    val /= static_cast<T>(variations.size());
  }

  return result;
}

template <typename T>
inline std::vector<T> sumInQuadratureOverVariations(const std::vector<std::vector<T>>& variations)
{
  if (variations.empty()) {
    return {};
  }

  const size_t nBins = variations[0].size();
  std::vector<T> result(nBins, T{});

  for (const auto& variation : variations) {
    if (variation.size() != nBins) {
      std::cerr << "Error in sumInQuadratureOverVariations: inconsistent vector sizes!" << std::endl;
      return {};
    }
    for (size_t ibin = 0; ibin < nBins; ++ibin) {
      result[ibin] += variation[ibin] * variation[ibin];
    }
  }

  for (auto& val : result) {
    val = std::sqrt(val);
  }

  return result;
}

} // namespace cent

#endif //CENTRALIYHELPER_H_