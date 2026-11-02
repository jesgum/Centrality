#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"

// <Npart>/<Ncoll> vs centrality for the PbPb anchor-point variations, with the
// ratio to the reference (90%) anchor point in the bottom pad
void doNpartNCollVsAnchorPoint(TH1F* hBase, std::vector<TH1F*> histos, TLegend* leg, const char* outfile, const char* titleYaxis = "<N_{part}>")
{
  TCanvas* canv = new TCanvas("canv", "", 1200, 800);
  canv->SetTopMargin(0.0);
  canv->SetBottomMargin(0.0);

  TPad* padTop = new TPad("padTop", "", 0, 0.3, 1, 1);
  padTop->SetTicks(1, 1);
  padTop->SetLogy(); // PbPb values span ~3 orders of magnitude from central to peripheral
  padTop->SetTopMargin(0.04);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.1);
  padTop->SetRightMargin(0.03);
  padTop->Draw();
  padTop->cd();
  hBase->GetXaxis()->SetLabelSize(0);
  hBase->GetYaxis()->SetLabelSize(0.05);
  hBase->GetYaxis()->SetTitle(titleYaxis);
  hBase->GetYaxis()->SetTitleSize(0.07);
  hBase->GetYaxis()->SetTitleOffset(0.6);
  hBase->SetMinimum(0.5 * hBase->GetMinimum());
  hBase->SetMaximum(2.0 * hBase->GetMaximum());
  hBase->Draw("pe");
  for (const auto& hist : histos) {
    hist->Draw("same pe");
  }
  leg->Draw();

  canv->cd();
  TPad* padBot = new TPad("padBot", "", 0, 0, 1, 0.3);
  padBot->SetTicks(1, 1);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.1);
  padBot->SetRightMargin(0.03);
  padBot->Draw();
  padBot->cd();

  std::vector<TH1F*> hRatios;
  double globalMin = 1e9, globalMax = -1e9;

  for (size_t i = 0; i < histos.size(); ++i) {
    TH1F* hRatio = (TH1F*)histos[i]->Clone(Form("hRatio_%zu", i));
    hRatio->Divide(hRatio, hBase, 1, 1, "B");
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Ratio to 90%");
    hRatio->GetXaxis()->SetTitle("Centrality (%)");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);

    // scan bin content +/- error to get the true drawn extent
    for (int b = 1; b <= hRatio->GetNbinsX(); ++b) {
      double content = hRatio->GetBinContent(b);
      double error = hRatio->GetBinError(b);
      if (content == 0 && error == 0) continue; // skip empty bins
      globalMin = std::min(globalMin, content - error);
      globalMax = std::max(globalMax, content + error);
    }

    hRatios.push_back(hRatio);
  }

  // add a bit of padding so points/error bars aren't flush with the pad edge
  double margin = (globalMax - globalMin) * 0.1;
  if (margin <= 0) margin = 0.01; // fallback in case all values are identical
  globalMin -= margin;
  globalMax += margin;

  for (auto* hRatio : hRatios) {
    hRatio->GetYaxis()->SetRangeUser(globalMin, globalMax);
  }

  hRatios[0]->Draw("pe");
  for (size_t i = 1; i < hRatios.size(); ++i) hRatios[i]->Draw("pe same");

  TLine* lUnity = new TLine(hRatios[0]->GetXaxis()->GetXmin(), 1, hRatios[0]->GetXaxis()->GetXmax(), 1);
  lUnity->SetLineStyle(kDashed);
  lUnity->SetLineColor(kGray + 2);
  lUnity->Draw();

  canv->SaveAs(outfile);
  for (auto* hRatio : hRatios) delete hRatio;
  delete lUnity;
  delete canv;
}

// Data and the Glauber + NBD fit stored by runCalibration (hGlauber: the fit
// rescaled to the data around the anchor point), plus the anchor point itself
struct GlauberFit {
  TH1* hData = nullptr;
  TH1* hGlauber = nullptr;
  float anchorPoint = -1;
  int anchorPercentage = 0;
  int color = kBlack;

  GlauberFit(const char* path, int percentage, int col)
    : anchorPercentage(percentage), color(col)
  {
    TFile* file = TFile::Open(Form("../AnalysisResults/%s", path), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: cannot open " << path << "!" << std::endl;
      return;
    }
    hData = dynamic_cast<TH1*>(file->Get("hData"));
    hGlauber = dynamic_cast<TH1*>(file->Get("hGlauber"));
    TH1* hAnchorPoint = dynamic_cast<TH1*>(file->Get("hAnchorPoint"));
    if (!hData || !hGlauber || !hAnchorPoint) {
      std::cerr << "Error: hData, hGlauber or hAnchorPoint missing in " << path << "!" << std::endl;
      hData = hGlauber = nullptr;
      file->Close();
      return;
    }
    hData->SetName(Form("hData_%d", percentage));
    hGlauber->SetName(Form("hGlauber_%d", percentage));
    hData->SetDirectory(nullptr);
    hGlauber->SetDirectory(nullptr);
    hGlauber->SetLineColor(color);
    hGlauber->SetLineWidth(2);
    anchorPoint = hAnchorPoint->GetBinContent(1);
    file->Close();
    delete file;
  }

  // data / fit, both rebinned by `rebin`
  TH1* ratio(int rebin, const char* tag) const
  {
    TH1* hNum = (TH1*)hData->Clone(Form("hRatio%s_%d", tag, anchorPercentage));
    TH1* hDen = (TH1*)hGlauber->Clone(Form("hDen%s_%d", tag, anchorPercentage));
    hNum->SetDirectory(nullptr);
    hDen->SetDirectory(nullptr);
    hNum->Rebin(rebin);
    hDen->Rebin(rebin);
    hNum->Divide(hDen);
    delete hDen;
    hNum->SetLineColor(color);
    hNum->SetLineWidth(2);
    hNum->SetTitle("");
    return hNum;
  }
};

void drawAnchorLines(const std::vector<GlauberFit>& fits, double ylo, double yhi)
{
  for (const auto& fit : fits) {
    TLine* line = new TLine(fit.anchorPoint, ylo, fit.anchorPoint, yhi);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->SetLineColorAlpha(fit.color, 0.8);
    line->Draw();
  }
}

// Same as in drawSummaryPlots.cc: the amplitude above which lPercentileRequested % of the events lie
Double_t GetBoundaryForPercentile(TH1* histo, Double_t lPercentileRequested)
{
  Double_t lReturnValue = 0.0;
  Double_t lPercentile = 100.0 - lPercentileRequested;

  const Long_t lNBins = histo->GetNbinsX();
  Double_t lCountDesired = lPercentile * histo->GetEntries() / 100;
  Long_t lCount = 0;
  for (Long_t ibin = 1; ibin < lNBins; ibin++) {
    lCount += histo->GetBinContent(ibin);
    if (lCount >= lCountDesired) {
      Double_t lWidth = histo->GetBinWidth(ibin);
      Double_t lLeftPercentile = 100. * (lCount - histo->GetBinContent(ibin)) / histo->GetEntries();
      Double_t lRightPercentile = 100. * lCount / histo->GetEntries();
      Double_t lProportion = (lPercentile - lLeftPercentile) / (lRightPercentile - lLeftPercentile);
      lReturnValue = histo->GetBinLowEdge(ibin) + lProportion * lWidth;
      break;
    }
  }
  return lReturnValue;
}

// Data with every fit overlaid and data / fit below, with the canvas, pads,
// insets and axes laid out as in drawSummaryPlots.cc. The anchor point of each
// fit is marked in both insets, which is where the variations differ.
void drawGlauberFits(const std::vector<GlauberFit>& fits, const char* outfile, const char* runNumber = "544122")
{
  const TString lXaxisTitle = "FT0C Amplitude";
  const int rebinRatio = 20;
  const Float_t lLoRangeRatio = 0.62;
  const Float_t lHiRangeRatio = 2.495;

  TH1* hV0M = (TH1*)fits.front().hData->Clone("hV0MDraw");
  hV0M->SetDirectory(nullptr);
  const Double_t lFitRangeMax = GetBoundaryForPercentile(hV0M, 0.002);
  const Double_t lZoomMaxTop = 590 * lFitRangeMax / 59000.;
  const Double_t lZoomMaxRatio = 1150 * lFitRangeMax / 59000.;
  double glauberMax = 0;
  for (const auto& fit : fits) {
    glauberMax = std::max(glauberMax, fit.hGlauber->GetMaximum());
  }

  TCanvas* c1 = new TCanvas("canvGlauber", "", 900, 900);
  c1->SetFrameFillStyle(0);
  c1->SetFillStyle(0);
  c1->Divide(1, 2);

  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1, 1);
  c1->cd(1)->SetPad(0, 0.5, 1, 1);
  c1->cd(2)->SetPad(0, 0.0, 1, .5);

  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.05);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.08);

  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.05);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.08);
  c1->cd(2)->SetTicks(1, 1);

  // ---- top: data and fits
  c1->cd(1);
  hV0M->GetXaxis()->SetRangeUser(0, lFitRangeMax);
  hV0M->GetYaxis()->SetRangeUser(0.25, hV0M->GetMaximum() * 3);
  hV0M->SetLineColor(kBlack);
  hV0M->GetYaxis()->SetTitleSize(0.07);
  hV0M->GetYaxis()->SetLabelSize(0.05);
  hV0M->GetYaxis()->SetTitle("Count");
  hV0M->GetYaxis()->SetTitleOffset(0.9);
  hV0M->GetXaxis()->SetLabelSize(0.05);
  hV0M->GetXaxis()->SetTitleSize(0.06);
  hV0M->GetXaxis()->SetTitle(lXaxisTitle.Data());
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->SetStats(0);
  hV0M->SetTitle("");
  hV0M->Draw("hist");
  for (const auto& fit : fits) {
    fit.hGlauber->Draw("same hist");
  }

  TLatex* latPer = new TLatex();
  latPer->SetNDC();
  latPer->SetTextFont(42);
  latPer->SetTextSize(0.06);
  latPer->DrawLatex(0.12, 0.88, "ALICE Pb-Pb 5.36 TeV");
  latPer->SetTextSize(0.062);
  latPer->DrawLatex(0.12, 0.81, Form("Run %s", runNumber));

  TLegend* legO = new TLegend(0.12, 0.5, 0.5, 0.78);
  legO->SetBorderSize(0);
  legO->SetFillColorAlpha(0, 0);
  legO->SetTextSize(0.045);
  legO->AddEntry(hV0M, "Data", "l");
  for (const auto& fit : fits) {
    legO->AddEntry(fit.hGlauber, Form("Glauber + NBD fit, anchor %d%% (%.0f)", fit.anchorPercentage, fit.anchorPoint), "l");
  }
  TLine* lLegAnchor = new TLine();
  lLegAnchor->SetLineStyle(7);
  lLegAnchor->SetLineWidth(2);
  lLegAnchor->SetLineColor(kGray + 2);
  legO->AddEntry(lLegAnchor, "Anchor point (insets)", "l");
  legO->Draw();

  TPad* inset1 = new TPad("inset1", "", 0.575724, 0.441648, 0.924, 0.961098);
  inset1->SetLogy();
  inset1->SetTicks(1, 1);
  inset1->SetBottomMargin(0.194);
  inset1->SetRightMargin(0.03);
  inset1->SetTopMargin(0.085);
  inset1->SetLeftMargin(0.11);
  inset1->Draw();
  inset1->cd();

  TH1* hV0Mzoomed = (TH1*)hV0M->Clone("hV0Mzoomed");
  hV0Mzoomed->SetDirectory(nullptr);
  hV0Mzoomed->GetXaxis()->SetRangeUser(0, lZoomMaxTop);
  hV0Mzoomed->GetXaxis()->SetNdivisions(404);
  hV0Mzoomed->GetXaxis()->SetTitleSize(0.1);
  hV0Mzoomed->GetYaxis()->UnZoom();
  const double insetLo = glauberMax * 0.05, insetHi = hV0Mzoomed->GetMaximum() * 5.5;
  hV0Mzoomed->GetYaxis()->SetRangeUser(insetLo, insetHi);
  hV0Mzoomed->GetYaxis()->SetTitleSize(0.1);
  hV0Mzoomed->GetYaxis()->SetTitleOffset(1.1);
  hV0Mzoomed->GetXaxis()->SetLabelSize(0.06);
  hV0Mzoomed->GetYaxis()->SetLabelSize(0.06);
  hV0Mzoomed->GetXaxis()->SetTitle("");
  hV0Mzoomed->GetYaxis()->SetTitle("");
  hV0Mzoomed->Draw("hist");
  for (const auto& fit : fits) {
    fit.hGlauber->Draw("same hist");
  }
  drawAnchorLines(fits, insetLo, insetHi);

  // ---- bottom: data / fit
  c1->cd(2);
  std::vector<TH1*> hRatios;
  for (const auto& fit : fits) {
    TH1* hRatio = fit.ratio(rebinRatio, "Main");
    hRatio->GetYaxis()->SetTitle("Data/Fit");
    hRatio->GetXaxis()->SetTitle(lXaxisTitle.Data());
    hRatio->GetYaxis()->SetTitleSize(0.055);
    hRatio->GetYaxis()->SetTitleOffset(0.7);
    hRatio->GetXaxis()->SetTitleSize(0.055);
    hRatio->GetYaxis()->SetLabelSize(0.045);
    hRatio->GetXaxis()->SetLabelSize(0.045);
    hRatio->GetYaxis()->SetRangeUser(lLoRangeRatio, lHiRangeRatio);
    hRatio->GetXaxis()->SetRangeUser(0, lFitRangeMax);
    hRatio->SetLineWidth(1);
    hRatio->SetStats(0);
    hRatio->Draw(hRatios.empty() ? "hist" : "same hist");
    hRatios.push_back(hRatio);
  }

  TLine* line = new TLine(0, 1, lFitRangeMax, 1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray + 1);
  line->Draw();

  TPad* inset2 = new TPad("inset2", "", 0.104677, 0.47032, 0.558775, 0.965753);
  inset2->SetTicks(1, 1);
  inset2->SetBottomMargin(0.194);
  inset2->SetRightMargin(0.03);
  inset2->SetTopMargin(0.045);
  inset2->SetLeftMargin(0.11);
  inset2->Draw();
  inset2->cd();

  for (size_t ii = 0; ii < hRatios.size(); ++ii) {
    TH1* hRatioZoomed = (TH1*)hRatios[ii]->Clone(Form("hRatioZoomed_%zu", ii));
    hRatioZoomed->SetDirectory(nullptr);
    hRatioZoomed->GetXaxis()->SetRangeUser(0, lZoomMaxRatio);
    hRatioZoomed->GetYaxis()->SetRangeUser(lLoRangeRatio, 1.55);
    hRatioZoomed->GetYaxis()->SetTitleSize(0.1);
    hRatioZoomed->GetXaxis()->SetTitleSize(0.1);
    hRatioZoomed->GetYaxis()->SetTitleOffset(1.1);
    hRatioZoomed->GetXaxis()->SetLabelSize(0.06);
    hRatioZoomed->GetYaxis()->SetLabelSize(0.06);
    hRatioZoomed->GetXaxis()->SetTitle("");
    hRatioZoomed->GetYaxis()->SetTitle("");
    hRatioZoomed->Draw(ii == 0 ? "hist" : "hist same");
  }

  TLine* line2 = new TLine(0, 1, lZoomMaxRatio, 1);
  line2->SetLineStyle(7);
  line2->SetLineColor(kGray + 1);
  line2->Draw();
  drawAnchorLines(fits, lLoRangeRatio, 1.55);

  c1->SaveAs(outfile);
}

void drawPbPbNpartNcollVsAnchorPoint()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const char* prefix = "LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb";
  syst::GlauberParameters anchor90(Form("%s_Anchor90.root", prefix)); // reference
  syst::GlauberParameters anchor89(Form("%s_Anchor89.root", prefix));
  syst::GlauberParameters anchor91(Form("%s_Anchor91.root", prefix));

  ColorManager cm(2);
  TH1F* hBaseNpart = syst::initGlauParHist(anchor90.npart, "hAnchor90Npart", kBlack);
  TH1F* hBaseNcoll = syst::initGlauParHist(anchor90.ncoll, "hAnchor90Ncoll", kBlack);
  TH1F* hAnchor89Npart = syst::initGlauParHist(anchor89.npart, "hAnchor89Npart", cm.getColor(0));
  TH1F* hAnchor91Npart = syst::initGlauParHist(anchor91.npart, "hAnchor91Npart", cm.getColor(1));
  TH1F* hAnchor89Ncoll = syst::initGlauParHist(anchor89.ncoll, "hAnchor89Ncoll", cm.getColor(0));
  TH1F* hAnchor91Ncoll = syst::initGlauParHist(anchor91.ncoll, "hAnchor91Ncoll", cm.getColor(1));

  std::vector<TH1F*> hNpart = { hAnchor89Npart, hAnchor91Npart };
  std::vector<TH1F*> hNcoll = { hAnchor89Ncoll, hAnchor91Ncoll };

  TLegend* leg = new TLegend(0.75, 0.6, 0.95, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hBaseNpart, "Anchor 90 (reference)", "pl");
  leg->AddEntry(hAnchor89Npart, "Anchor 89", "pl");
  leg->AddEntry(hAnchor91Npart, "Anchor 91", "pl");

  doNpartNCollVsAnchorPoint(hBaseNpart, hNpart, leg, "FiguresSystematics/hNpartVsAnchorPbPb.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hBaseNcoll, hNcoll, leg, "FiguresSystematics/hNcollVsAnchorPbPb.pdf", "<N_{coll}>");

  std::vector<GlauberFit> fits;
  fits.emplace_back(Form("%s_Anchor89.root", prefix), 89, kBlue + 1);
  fits.emplace_back(Form("%s_Anchor90.root", prefix), 90, kRed); // reference, red as in drawSummaryPlots
  fits.emplace_back(Form("%s_Anchor91.root", prefix), 91, kGreen + 2);
  for (const auto& fit : fits) {
    if (!fit.hData) {
      return;
    }
  }
  drawGlauberFits(fits, "FiguresSystematics/hGlauberFitsVsAnchorPbPb.pdf");
}
