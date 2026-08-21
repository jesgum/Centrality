#ifndef CENTRALIYHELPER_H_
#define CENTRALIYHELPER_H_

#include "TH1.h"

namespace cent
{

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

} // namespace cent

#endif //CENTRALIYHELPER_H_