#ifndef COLORMANAGER_H_
#define COLORMANAGER_H_

#include <vector>
#include "TColor.h"

class ColorManager
{
 public:
  enum class Palette {
    Default, // the original rainbow-like palette
    Red,
    Blue,
    Green,
    Orange,
    Purple,
    Gray,
    CoolWarm, // blue -> white -> red
    Custom // set via setCustomPalette()
  };

  ColorManager() = default;

  explicit ColorManager(size_t nColors, Palette palette = Palette::Default, bool reversed = false)
    : mPalette(palette), mReversed(reversed)
  {
    reset(nColors);
  }

  void setPalette(Palette palette, bool reversed = false)
  {
    mPalette = palette;
    mReversed = reversed;
    if (!mColors.empty()) {
      reset(mColors.size());
    }
  }

  // Define your own gradient. Stops must be increasing in [0, 1] and all vectors
  // must have the same length; RGB values are in [0, 1].
  void setCustomPalette(const std::vector<double>& stops, const std::vector<double>& red,
                        const std::vector<double>& green, const std::vector<double>& blue,
                        bool reversed = false)
  {
    mCustom = { stops, red, green, blue };
    setPalette(Palette::Custom, reversed);
  }

  void reset(size_t nColors)
  {
    mColors.clear();
    mColors.resize(nColors);

    Gradient g = getGradient(mPalette);
    mFirstIndex = TColor::CreateGradientColorTable(g.stops.size(), g.stops.data(), g.red.data(),
                                                   g.green.data(), g.blue.data(), nColors);
    for (size_t i = 0; i < nColors; ++i) {
      // Default ordering (not reversed) matches the original: high end of the gradient first
      mColors[i] = mReversed ? mFirstIndex + i : mFirstIndex + nColors - 1 - i;
    }
  }

  int getColor(size_t i) const
  {
    if (i >= mColors.size() || i < 0) {
      return mColStd;
    }

    return mColors[i];
  }

  size_t size() const { return mColors.size(); }

 private:
  struct Gradient {
    std::vector<double> stops, red, green, blue;
  };

  // Converts 0-255 RGB anchors at equally spaced stops into a Gradient
  static Gradient fromRGB255(const std::vector<std::vector<int>>& rgb)
  {
    Gradient g;
    const size_t n = rgb.size();
    for (size_t i = 0; i < n; ++i) {
      g.stops.push_back(n > 1 ? double(i) / (n - 1) : 0.);
      g.red.push_back(rgb[i][0] / 255.);
      g.green.push_back(rgb[i][1] / 255.);
      g.blue.push_back(rgb[i][2] / 255.);
    }
    return g;
  }

  Gradient getGradient(Palette palette) const
  {
    switch (palette) {
      // Single-hue palettes: light -> dark (the light end is kept visible on white)
      case Palette::Red:
        return fromRGB255({ { 252, 146, 114 }, { 239, 59, 44 }, { 165, 15, 21 }, { 80, 0, 10 } });
      case Palette::Blue:
        return fromRGB255({ { 158, 202, 225 }, { 66, 146, 198 }, { 8, 81, 156 }, { 8, 30, 80 } });
      case Palette::Green:
        return fromRGB255({ { 161, 217, 155 }, { 65, 171, 93 }, { 0, 109, 44 }, { 0, 55, 20 } });
      case Palette::Orange:
        return fromRGB255({ { 253, 174, 107 }, { 241, 105, 19 }, { 166, 54, 3 }, { 90, 30, 0 } });
      case Palette::Purple:
        return fromRGB255({ { 188, 189, 220 }, { 128, 125, 186 }, { 84, 39, 143 }, { 45, 0, 80 } });
      case Palette::Gray:
        return fromRGB255({ { 200, 200, 200 }, { 0, 0, 0 } });
      case Palette::CoolWarm:
        return fromRGB255({ { 59, 76, 192 }, { 221, 221, 221 }, { 180, 4, 38 } });
      case Palette::Custom:
        if (!mCustom.stops.empty()) {
          return mCustom;
        }
        [[fallthrough]];
      case Palette::Default:
      default:
        return { { 0.00, 0.34, 0.61, 0.84, 1.00 },
                 { 0.00, 0.00, 0.9 * 0.87, 1.00, 0.51 },
                 { 0.00, 0.81, 0.9 * 1.00, 0.20, 0.00 },
                 { 0.51, 0.9 * 1.00, 0.12, 0.00, 0.00 } };
    }
  }

  Palette mPalette = Palette::Default;
  bool mReversed = false;
  Gradient mCustom;

  int mColStd = kBlack;
  int mFirstIndex = 0;
  std::vector<int> mColors;
};

#endif //COLORMANAGER_H_
