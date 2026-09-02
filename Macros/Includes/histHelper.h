#ifndef HISTHELPER_H_
#define HISTHELPER_H_

#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"

#include <array>

namespace hist {

template <typename THist2D = TH2F>
inline THist2D* createTH2(const char* name, const char* title, std::array<float, 3> xaxis, std::array<float, 3> yaxis)
{
  return new THist2D(name, title, xaxis[0], xaxis[1], xaxis[2], yaxis[0], yaxis[1], yaxis[2]);
}

template <typename THist = TH1F>
inline THist* createTH1(const char* name, const char* title, std::array<float, 3> axis)
{
  return new THist(name, title, axis[0], axis[1], axis[2]);
}

inline void drawTH1(TH1* h1d, const char* title)
{
  TCanvas* canv = new TCanvas(Form("canv%s", title), "", 1600, 1200);
  canv->SetLogy();
  h1d->Draw("");
  canv->SaveAs(Form("h%s.pdf", title));
}

inline void drawTH2(TH2* h2d, const char* title)
{
  TCanvas* canv = new TCanvas(Form("canv%s", title), "", 1600, 1200);
  canv->SetLogz();
  h2d->Draw("colz");
  canv->SaveAs(Form("h2d%s.pdf", title));
}

} // namespace hist


#endif // HISTHELPER_H_