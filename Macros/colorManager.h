#ifndef COLORMANAGER_H_
#define COLORMANAGER_H_

#include <vector>
#include "TColor.h"

class ColorManager
{
 public:
  ColorManager() = default;

  explicit ColorManager(size_t nColors)
  {
    reset(nColors);
  }

  void reset(size_t nColors)
  {
    mColors.clear();
    mColors.resize(nColors);

    mFirstIndex = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, nColors);
    for (size_t i = 0; i < nColors; ++i) {
      mColors[i] = mFirstIndex + nColors - 1 - i;
    }
  }

  int getColor(size_t i) const
  {
    if (i >= mColors.size()) {
      return mColStd;
    }

    return mColors[i];
  }

 private:
  static constexpr int NRGBs = 5;

  double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  double red[NRGBs] = { 0.00, 0.00, 0.9 * 0.87, 1.00, 0.51 };
  double green[NRGBs] = { 0.00, 0.81, 0.9 * 1.00, 0.20, 0.00 };
  double blue[NRGBs] = { 0.51, 0.9 * 1.00, 0.12, 0.00, 0.00 };

  int mColStd = kBlack;
  int mFirstIndex = 0;
  std::vector<int> mColors;
};

#endif //COLORMANAGER_H_