#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <vector>

// Draws the Glauber + NBD fits of several calibration files (output of Glauber/runCalibration.cc) on the same
// pad together with the data, bottom pad: data / fit, both with an inset zooming on the peripheral region (layout
// as in drawSummaryPlots.cc). The data is taken from the first file. The anchor point of every calibration is
// marked with a dotted line in the color of its fit.

struct GlauberInput {
  TString fileName;
  TString label;
  int color;
};

struct GlauberCurve {
  TH1* hData = nullptr;
  TH1* hFit = nullptr;   // in counts per data bin
  TH1* hRatio = nullptr; // data / fit, rebinned
  TH1* hRatioFine = nullptr; // data / fit in data bins, for the zoomed inset
  double anchorPoint = -1;
};

bool loadGlauber(const GlauberInput& input, int index, int rebinRatio, GlauberCurve& curve)
{
  TFile* file = TFile::Open(input.fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: cannot open " << input.fileName << std::endl;
    return false;
  }
  TH1* hData = dynamic_cast<TH1*>(file->Get("hData"));
  TH1* hFit = dynamic_cast<TH1*>(file->Get("hGlauberFine"));
  if (!hFit) {
    hFit = dynamic_cast<TH1*>(file->Get("hGlauber"));
  }
  if (!hData || !hFit) {
    std::cerr << "Error: hData or hGlauberFine/hGlauber not found in " << input.fileName << std::endl;
    return false;
  }
  curve.hData = dynamic_cast<TH1*>(hData->Clone(Form("hData_%d", index)));
  curve.hFit = dynamic_cast<TH1*>(hFit->Clone(Form("hFit_%d", index)));
  if (TH1* hAnchor = dynamic_cast<TH1*>(file->Get("hAnchorPoint"))) {
    curve.anchorPoint = hAnchor->GetBinContent(1);
  }
  file->Close();

  // data / fit, both rebinned to rebinRatio data bins (the fit may be finer than the data)
  const int rebinFit = std::max(1, TMath::Nint(curve.hData->GetBinWidth(1) / curve.hFit->GetBinWidth(1)));
  curve.hRatio = dynamic_cast<TH1*>(curve.hData->Clone(Form("hRatio_%d", index)));
  curve.hRatio->Rebin(rebinRatio);
  TH1* hFitRebinned = dynamic_cast<TH1*>(curve.hFit->Clone(Form("hFitRebinned_%d", index)));
  hFitRebinned->Rebin(rebinRatio * rebinFit);
  curve.hRatio->Divide(hFitRebinned);
  delete hFitRebinned;
  curve.hRatioFine = dynamic_cast<TH1*>(curve.hData->Clone(Form("hRatioFine_%d", index)));
  TH1* hFitFine = dynamic_cast<TH1*>(curve.hFit->Clone(Form("hFitFine_%d", index)));
  hFitFine->Rebin(rebinFit);
  curve.hRatioFine->Divide(hFitFine);
  delete hFitFine;

  // fit in counts per data bin, to be drawn on top of the data
  curve.hFit->Scale(curve.hData->GetBinWidth(1) / curve.hFit->GetBinWidth(1));
  return true;
}

// Draws the fits of inputs (any calibration files with hData, hGlauber(Fine) and hAnchorPoint: Glauber or Trento)
void drawMultipleFits(const std::vector<GlauberInput>& inputs, const TString& outFile)
{
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  const TString xTitle = "FT0C Amplitude";
  const TString header = "ALICE Pb-Pb 5.36 TeV";
  const double xMax = 54000;
  const int rebinRatio = 20;
  const double ratioMin = 0.45, ratioMax = 2.495;  // the upper half of the ratio pad is covered by the inset
  const double xMaxInsetTop = 590 * xMax / 59000.;  // peripheral zoom, as in drawSummaryPlots.cc
  const double xMaxInsetBot = 1150 * xMax / 59000.;
  const double ratioMinInset = 0.45, ratioMaxInset = 1.55;
  const bool drawAnchorPoints = true;

  std::vector<GlauberCurve> curves;
  for (size_t i = 0; i < inputs.size(); ++i) {
    GlauberCurve curve;
    if (!loadGlauber(inputs[i], i, rebinRatio, curve)) {
      return;
    }
    curves.push_back(curve);
  }
  if (curves.empty()) {
    std::cerr << "Error: no calibration files given" << std::endl;
    return;
  }

  // canvas split in two equal halves, as in drawSummaryPlots.cc
  TCanvas* canv = new TCanvas("canvMultipleGlaubers", "", 900, 900);
  TPad* padTop = new TPad("padTop", "", 0, 0.5, 1, 1);
  padTop->SetTicks(1, 1);
  padTop->SetLogy();
  padTop->SetTopMargin(0.02);
  padTop->SetBottomMargin(0.001);
  padTop->SetLeftMargin(0.08);
  padTop->SetRightMargin(0.05);
  padTop->Draw();
  canv->cd();
  TPad* padBot = new TPad("padBot", "", 0, 0, 1, 0.5);
  padBot->SetTicks(1, 1);
  padBot->SetTopMargin(0.001);
  padBot->SetBottomMargin(0.14);
  padBot->SetLeftMargin(0.08);
  padBot->SetRightMargin(0.05);
  padBot->Draw();

  // anchor points as dotted lines from yMin to yMax in the current pad
  auto drawAnchorLines = [&](double yMin, double yMax) {
    if (!drawAnchorPoints) {
      return;
    }
    for (size_t i = 0; i < curves.size(); ++i) {
      if (curves[i].anchorPoint > 0) {
        TLine* line = new TLine(curves[i].anchorPoint, yMin, curves[i].anchorPoint, yMax);
        line->SetLineStyle(kDotted);
        line->SetLineWidth(2);
        line->SetLineColor(inputs[i].color);
        line->Draw();
      }
    }
  };

  // ---- top: data and fits ----
  padTop->cd();
  TH1* hData = curves[0].hData;
  hData->SetTitle("");
  hData->SetLineColor(kBlack);
  hData->GetXaxis()->SetRangeUser(0, xMax);
  hData->SetMinimum(0.25);
  hData->SetMaximum(hData->GetMaximum() * 3);
  hData->GetYaxis()->SetTitle("Count");
  hData->GetYaxis()->SetTitleSize(0.07);
  hData->GetYaxis()->SetTitleOffset(0.9);
  hData->GetYaxis()->SetLabelSize(0.05);
  hData->GetYaxis()->SetTickLength(0.015);
  hData->Draw("hist");
  for (size_t i = 0; i < curves.size(); ++i) {
    curves[i].hFit->SetLineColor(inputs[i].color);
    curves[i].hFit->SetLineWidth(2);
    curves[i].hFit->Draw("hist same");
  }
  drawAnchorLines(hData->GetMinimum(), hData->GetMaximum());

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.06);
  latex->DrawLatex(0.12, 0.88, header);

  // left, below the header: the right of the pad is taken by the inset
  TLegend* leg = new TLegend(0.12, 0.84 - 0.06 * (curves.size() + 1 + drawAnchorPoints), 0.55, 0.84);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->SetTextSize(0.05);
  leg->AddEntry(hData, "Data", "l");
  for (size_t i = 0; i < curves.size(); ++i) {
    leg->AddEntry(curves[i].hFit, inputs[i].label, "l");
  }
  if (drawAnchorPoints) {
    TLine* anchorLegend = new TLine(); // only for the legend, in the color of the fit in the plot
    anchorLegend->SetLineStyle(kDotted);
    anchorLegend->SetLineWidth(2);
    anchorLegend->SetLineColor(kGray + 2);
    leg->AddEntry(anchorLegend, "Anchor point", "l");
  }
  leg->Draw();

  // inset: peripheral region
  TPad* insetTop = new TPad("insetTop", "", 0.575724, 0.441648, 0.924, 0.961098);
  insetTop->SetLogy();
  insetTop->SetTicks(1, 1);
  insetTop->SetBottomMargin(0.194);
  insetTop->SetRightMargin(0.03);
  insetTop->SetTopMargin(0.085);
  insetTop->SetLeftMargin(0.11);
  insetTop->SetFillStyle(0); // transparent, so the histograms of the pad below stay visible
  insetTop->SetFrameFillStyle(0);
  insetTop->Draw();
  insetTop->cd();
  TH1* hDataZoomed = dynamic_cast<TH1*>(hData->Clone("hDataZoomed"));
  hDataZoomed->GetXaxis()->SetRangeUser(0, xMaxInsetTop);
  hDataZoomed->SetMinimum(-1111); // reset, then set from the zoomed range
  hDataZoomed->SetMaximum(-1111);
  double fitMaxZoomed = 0;
  for (const auto& curve : curves) {
    curve.hFit->GetXaxis()->SetRangeUser(0, xMaxInsetTop);
    fitMaxZoomed = std::max(fitMaxZoomed, curve.hFit->GetMaximum());
    curve.hFit->GetXaxis()->SetRangeUser(0, xMax);
  }
  hDataZoomed->SetMinimum(fitMaxZoomed * 0.05);
  hDataZoomed->SetMaximum(hDataZoomed->GetMaximum() * 5.5);
  hDataZoomed->GetXaxis()->SetNdivisions(404);
  hDataZoomed->GetXaxis()->SetLabelSize(0.06);
  hDataZoomed->GetYaxis()->SetLabelSize(0.06);
  hDataZoomed->GetYaxis()->SetTitleSize(0.1);
  hDataZoomed->GetYaxis()->SetTitleOffset(1.1);
  hDataZoomed->GetYaxis()->SetTickLength(0.03);
  hDataZoomed->GetYaxis()->SetTitle("");
  hDataZoomed->Draw("hist");
  for (const auto& curve : curves) {
    curve.hFit->Draw("hist same");
  }
  drawAnchorLines(hDataZoomed->GetMinimum(), hDataZoomed->GetMaximum());

  // ---- bottom: data / fit ----
  padBot->cd();
  for (size_t i = 0; i < curves.size(); ++i) {
    for (TH1* hRatio : { curves[i].hRatio, curves[i].hRatioFine }) {
      hRatio->SetTitle("");
      hRatio->SetLineColor(inputs[i].color);
    }
  }
  TH1* hFrame = curves[0].hRatio;
  hFrame->GetXaxis()->SetRangeUser(0, xMax);
  hFrame->SetMinimum(ratioMin);
  hFrame->SetMaximum(ratioMax);
  hFrame->GetYaxis()->SetTitle("Data / Fit");
  hFrame->GetYaxis()->SetTitleSize(0.055);
  hFrame->GetYaxis()->SetTitleOffset(0.7);
  hFrame->GetYaxis()->SetLabelSize(0.045);
  hFrame->GetXaxis()->SetTitle(xTitle);
  hFrame->GetXaxis()->SetTitleSize(0.055);
  hFrame->GetXaxis()->SetLabelSize(0.045);
  hFrame->Draw("hist");
  TLine* unity = new TLine(0, 1, xMax, 1);
  unity->SetLineStyle(7);
  unity->SetLineColor(kGray + 1);
  unity->Draw();
  for (const auto& curve : curves) {
    curve.hRatio->Draw("hist same");
  }

  // inset: peripheral region, in data bins
  TPad* insetBot = new TPad("insetBot", "", 0.104677, 0.47032, 0.558775, 0.965753);
  insetBot->SetTicks(1, 1);
  insetBot->SetBottomMargin(0.194);
  insetBot->SetRightMargin(0.03);
  insetBot->SetTopMargin(0.045);
  insetBot->SetLeftMargin(0.11);
  insetBot->SetFillStyle(0);
  insetBot->SetFrameFillStyle(0);
  insetBot->Draw();
  insetBot->cd();
  TH1* hFrameZoomed = curves[0].hRatioFine;
  hFrameZoomed->GetXaxis()->SetRangeUser(0, xMaxInsetBot);
  hFrameZoomed->SetMinimum(ratioMinInset);
  hFrameZoomed->SetMaximum(ratioMaxInset);
  hFrameZoomed->GetXaxis()->SetLabelSize(0.06);
  hFrameZoomed->GetYaxis()->SetLabelSize(0.06);
  hFrameZoomed->Draw("hist");
  TLine* unityZoomed = new TLine(0, 1, xMaxInsetBot, 1);
  unityZoomed->SetLineStyle(7);
  unityZoomed->SetLineColor(kGray + 1);
  unityZoomed->Draw();
  drawAnchorLines(ratioMinInset, ratioMaxInset);
  for (const auto& curve : curves) {
    curve.hRatioFine->Draw("hist same");
  }

  canv->SaveAs(outFile);
}

void drawMultipleGlaubers()
{
  // ---- fits to compare: {calibration file, legend label, color} ----
  const std::vector<GlauberInput> inputs = {
    { "../AnalysisResults/LHC23_pass5_small/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_Collisions_basehistos_Pb_Anchor80.0.root", "Glauber, anchor 80%", kBlue + 1 },
    { "../AnalysisResults/LHC23_pass5_small/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_Collisions_basehistos_Pb_Anchor90.0.root", "Glauber, anchor 90%", kRed + 1 },
  };
  drawMultipleFits(inputs, "MultipleGlaubers.pdf");
}
