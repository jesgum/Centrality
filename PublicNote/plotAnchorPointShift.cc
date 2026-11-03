#include <cstdio>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TString.h"

#include "plotHelpers.h"

// Pb-Pb FT0C amplitude near the 90% anchor point: data, the Glauber + NBD fit and the
// renormalised Glauber + NBD, with a vertical line at the anchor point of each.
// Input: Files/hAnchorPointShift_zoom.root (from Glauber/drawSummaryPlotsAnchorPoint.cc),
// canvas c2 (log y), contents:
//   hV0MzoomedFull  data (frame)
//   hFittedGlauber  Glauber + NBD fit (blue)
//   hCalibGlauber   Glauber + NBD renormalised (red)
//   TLine x2        anchor points (fit in blue, renormalised in red)
//   TLatex          "ALICE Pb-Pb 5.36 TeV"
//   TLegend         one entry per histogram and line

void plotAnchorPointShift()
{
  plot::init();

  TFile* f = TFile::Open("Files/hAnchorPointShift_zoom.root");
  if (!f || f->IsZombie()) {
    printf("cannot open Files/hAnchorPointShift_zoom.root\n");
    return;
  }
  TCanvas* c = (TCanvas*)f->Get("c2");
  c->SetName("canvAnchorPointShift");
  c->Draw();
  plot::removeStats(c);

  // Canvas size in pixels (original 900x600)
  // c->SetCanvasSize(900, 600);
  // plot::setMargins(c, 0.10, 0.03, 0.12, 0.03); // left, right, bottom, top

  printf("\n=== anchor point shift ===\n");
  plot::printCanvasContents(c);

  TH1* hData = plot::get<TH1>(c, "hV0MzoomedFull"); // also the frame
  TH1* hFit = plot::get<TH1>(c, "hFittedGlauber");
  TH1* hRenorm = plot::get<TH1>(c, "hCalibGlauber");

  hData->SetTitle("");
  hData->GetXaxis()->SetTitle("FT0C Amplitude");
  hData->GetYaxis()->SetTitle("Counts");
  // hData->GetXaxis()->SetRangeUser(0, 540);
  // hData->GetYaxis()->SetRangeUser(1.25e+3, 2.4e+7); // log y; also update the line heights below

  // Line widths of the histograms (the legend follows them, see below)
  const int lineWidth = 2;
  for (const auto& h : { hData, hFit, hRenorm }) {
    h->SetLineWidth(lineWidth);
  }

  // Anchor-point lines: both are named "TLine", fit first, renormalised second
  for (const auto& obj : *c->GetListOfPrimitives()) {
    if (obj->IsA() == TLine::Class()) {
      TLine* line = (TLine*)obj;
      line->SetLineWidth(lineWidth);
      // line->SetY2(1.2e+5); // top of the line
    }
  }

  // Legend: position in NDC (x1, y1, x2, y2), originals from the file
  TLegend* leg = plot::get<TLegend>(c, "TPave");
  leg->SetX1NDC(0.5);
  leg->SetY1NDC(0.62);
  leg->SetX2NDC(0.92);
  leg->SetY2NDC(0.94);
  // leg->SetTextSize(0.038);
  // A legend entry is drawn with the line style of the object it points to. In the
  // file the histogram entries point to separate copies of the histograms, so point
  // them to the ones on the canvas to follow their style (entries: data, fit,
  // renormalised, then the two anchor lines)
  TList* entries = leg->GetListOfPrimitives();
  ((TLegendEntry*)entries->At(0))->SetObject(hData);
  ((TLegendEntry*)entries->At(1))->SetObject(hFit);
  ((TLegendEntry*)entries->At(2))->SetObject(hRenorm);
  for (const auto& obj : *entries) {
    ((TLegendEntry*)obj)->SetLineWidth(lineWidth);
  }

  // Replace the ALICE label
  plot::removeLatex(c, "ALICE");
  c->cd();
  TLatex* lat = new TLatex();
  // lat->SetTextFont(42);
  lat->SetTextSize(0.045);
  lat->DrawLatexNDC(0.12, 0.89, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");

  c->Modified();
  c->Update();
  c->SaveAs("Figures/hAnchorPointShift_zoom_PbPb.pdf");
}
