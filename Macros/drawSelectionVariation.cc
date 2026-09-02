#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "Includes/histHelper.h"
#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"

void doNpartNCollVsAnchorPoint(TH1F* hBase, std::vector<TH1F*> histos, TLegend* leg, const char* outfile, const char* titleYaxis = "<N_{part}>")
{
  std::vector<std::pair<TH1F*, TH1F*>> runPairs;
  TCanvas* canv = new TCanvas("canv", "", 1200, 800);
  canv->SetTopMargin(0.0);
  canv->SetBottomMargin(0.0);

  TPad* padTop = new TPad("padTop", "", 0, 0.3, 1, 1);
  padTop->SetTicks(1, 1);
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
  hBase->Draw("pe");
  for (const auto& hist : histos) {
    hist->Draw("same pe");
    runPairs.push_back({hist, hBase});
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

  for (int i = 0; i < runPairs.size(); ++i) {
    TH1F* hRatio = (TH1F*)runPairs[i].first->Clone(Form("hRatioNpart_%d", i));
    hRatio->Divide(hRatio, runPairs[i].second, 1, 1, "B");
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Ratio");
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
      double error   = hRatio->GetBinError(b);
      if (content == 0 && error == 0) continue; // skip empty bins
      double lo = content - error;
      double hi = content + error;
      if (lo < globalMin) globalMin = lo;
      if (hi > globalMax) globalMax = hi;
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
  for (int i = 1; i < runPairs.size(); ++i) hRatios[i]->Draw("pe same");

  canv->SaveAs(outfile);
  for (auto* hRatio : hRatios) delete hRatio;
  hRatios.clear();
  delete canv;
}

struct Registry {
  Registry(TFile* infile)
  {
    hData = (TH1F*)infile->Get("hV0MUltraFine");
    hGlau = (TH1F*)infile->Get("hGlauberFine");
  }

  TH1F* hData = nullptr;
  TH1F* hGlau = nullptr;
};

// Truncate histograms to a common number of bins (assumes identical bin width & xmin)
void matchBinning(std::vector<TH1F*>& hists)
{
  int nbinsCommon = hists[0]->GetNbinsX();
  for (auto* h : hists)
    nbinsCommon = std::min(nbinsCommon, h->GetNbinsX());

  for (auto*& h : hists) {
    if (h->GetNbinsX() == nbinsCommon)
      continue;
    TString newName = TString(h->GetName()) + "_trunc";
    double xmin = h->GetXaxis()->GetXmin();
    double binWidth = h->GetXaxis()->GetBinWidth(1);
    double xmax = xmin + nbinsCommon * binWidth;
    TH1F* hNew = new TH1F(newName, h->GetTitle(), nbinsCommon, xmin, xmax);
    for (int i = 1; i <= nbinsCommon; ++i) {
      hNew->SetBinContent(i, h->GetBinContent(i));
      hNew->SetBinError(i, h->GetBinError(i));
    }
    h = hNew;
  }
}

void drawSelectionVariation()
{
  gStyle->SetOptStat(0);
  Registry vtxZ = Registry(new TFile("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374_glauberNBD_hFT0M_Collisions_vtxZ.root", "read"));
  Registry noPileup = Registry(new TFile("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374_glauberNBD_hFT0M_Collisions_nopileup.root", "read"));
  Registry noEvSel = Registry(new TFile("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374_glauberNBD_hFT0M_Collisions_nosel.root", "read"));
  Registry recommended = Registry(new TFile("../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root", "read"));

  // Fix mismatched binning (recommended has more/finer bins than the others)
  std::vector<TH1F*> dataHists = { vtxZ.hData, noPileup.hData, noEvSel.hData, recommended.hData };
  matchBinning(dataHists);
  vtxZ.hData = dataHists[0];
  noPileup.hData = dataHists[1];
  noEvSel.hData = dataHists[2];
  recommended.hData = dataHists[3];

  std::vector<TH1F*> glauHists = { vtxZ.hGlau, noPileup.hGlau, noEvSel.hGlau, recommended.hGlau };
  matchBinning(glauHists);
  vtxZ.hGlau = glauHists[0];
  noPileup.hGlau = glauHists[1];
  noEvSel.hGlau = glauHists[2];
  recommended.hGlau = glauHists[3];

  vtxZ.hData->SetLineColor(kBlue + 1);
  noPileup.hData->SetLineColor(kRed + 1);
  noEvSel.hData->SetLineColor(kGreen + 2);
  recommended.hData->SetLineColor(kBlack);

  vtxZ.hData->SetMarkerColor(kBlue + 1);
  noPileup.hData->SetMarkerColor(kRed + 1);
  noEvSel.hData->SetMarkerColor(kGreen + 2);
  recommended.hData->SetMarkerColor(kBlack);

  vtxZ.hGlau->SetLineColor(kBlue + 1);
  noPileup.hGlau->SetLineColor(kRed + 1);
  noEvSel.hGlau->SetLineColor(kGreen + 2);
  recommended.hGlau->SetLineColor(kBlack);

  vtxZ.hGlau->SetMarkerColor(kBlue + 1);
  noPileup.hGlau->SetMarkerColor(kRed + 1);
  noEvSel.hGlau->SetMarkerColor(kGreen + 2);
  recommended.hGlau->SetMarkerColor(kBlack);

  const float integralLow = 5000.f;
  const float integralHigh = 10000.f;

  recommended.hData->Scale(1. / recommended.hData->Integral());
  vtxZ.hData->Scale(1. / vtxZ.hData->Integral());
  noPileup.hData->Scale(1. / noPileup.hData->Integral());
  noEvSel.hData->Scale(1. / noEvSel.hData->Integral());
  recommended.hGlau->Scale(1. / recommended.hGlau->Integral());
  vtxZ.hGlau->Scale(1. / vtxZ.hGlau->Integral());
  noPileup.hGlau->Scale(1. / noPileup.hGlau->Integral());
  noEvSel.hGlau->Scale(1. / noEvSel.hGlau->Integral());

  vtxZ.hData->Rebin(20);
  noPileup.hData->Rebin(20);
  noEvSel.hData->Rebin(20);
  recommended.hData->Rebin(20);
  vtxZ.hGlau->Rebin(20);
  noPileup.hGlau->Rebin(20);
  noEvSel.hGlau->Rebin(20);
  recommended.hGlau->Rebin(20);

  // Full
  const float xmin_full = 0;
  const float xmax_full = 31000;
  const float ymin_full = 2e-9;
  const float ymax_full = 4e-2;

  // Zoom
  const float xmin_zoom = 0;
  const float xmax_zoom = 1000;
  const float ymin_zoom = 2e-3;
  const float ymax_zoom = 4e-2;

  const bool doZoom = false;
  const float xmin = doZoom ? xmin_zoom : xmin_full;
  const float xmax = doZoom ? xmax_zoom : xmax_full;
  const float ymin = doZoom ? ymin_zoom : ymin_full;
  const float ymax = doZoom ? ymax_zoom : ymax_full;

  TCanvas* canvData = new TCanvas("canvData", "", 1200, 1000);
  canvData->SetLogy();
  recommended.hData->GetXaxis()->SetRangeUser(xmin, xmax);
  recommended.hData->GetYaxis()->SetRangeUser(ymin, ymax);

  recommended.hData->Draw("hist");
  vtxZ.hData->Draw("hist same");
  noPileup.hData->Draw("hist same");
  noEvSel.hData->Draw("hist same");
  canvData->SaveAs("hSelectionComparisonData.pdf");

  TCanvas* canvGlau = new TCanvas("canvGlau", "", 1200, 1000);
  canvGlau->SetLogy();
  recommended.hGlau->GetXaxis()->SetRangeUser(xmin, xmax);
  recommended.hGlau->GetYaxis()->SetRangeUser(ymin, ymax);

  recommended.hGlau->Draw("hist");
  vtxZ.hGlau->Draw("hist same");
  noPileup.hGlau->Draw("hist same");
  noEvSel.hGlau->Draw("hist same");
  canvGlau->SaveAs("hSelectionComparisonGlau.pdf");

  const float ratioYmin = 0.5;
  const float ratioYmax = 1.5;

  TH1F* hRatioVtxZData = (TH1F*)vtxZ.hData->Clone("hRatioVtxZData");
  TH1F* hRatioNoPileupData = (TH1F*)noPileup.hData->Clone("hRatioNoPileupData");
  TH1F* hRatioNoEvSelData = (TH1F*)noEvSel.hData->Clone("hRatioNoEvSelData");

  hRatioVtxZData->Divide(recommended.hData);
  hRatioNoPileupData->Divide(recommended.hData);
  hRatioNoEvSelData->Divide(recommended.hData);

  TH1F* hRatioVtxZGlau = (TH1F*)vtxZ.hGlau->Clone("hRatioVtxZGlau");
  TH1F* hRatioNoPileupGlau = (TH1F*)noPileup.hGlau->Clone("hRatioNoPileupGlau");
  TH1F* hRatioNoEvSelGlau = (TH1F*)noEvSel.hGlau->Clone("hRatioNoEvSelGlau");

  hRatioVtxZGlau->Divide(recommended.hGlau);
  hRatioNoPileupGlau->Divide(recommended.hGlau);
  hRatioNoEvSelGlau->Divide(recommended.hGlau);

  TCanvas* canvDataRatio = new TCanvas("canvDataRatio", "", 1200, 1000);
  hRatioVtxZData->GetXaxis()->SetRangeUser(xmin, xmax);
  hRatioVtxZData->GetYaxis()->SetRangeUser(ratioYmin, ratioYmax);
  hRatioVtxZData->GetYaxis()->SetTitle("variation / recommended");
  hRatioVtxZData->SetTitle("");

  hRatioVtxZData->Draw("hist");
  hRatioNoPileupData->Draw("hist same");
  hRatioNoEvSelData->Draw("hist same");

  TLine* lineData = new TLine(xmin, 1, xmax, 1);
  lineData->SetLineStyle(2);
  lineData->SetLineColor(kGray + 2);
  lineData->Draw("same");

  TLegend* legData = new TLegend(0.6, 0.7, 0.88, 0.88);
  legData->SetBorderSize(0);
  legData->SetFillStyle(0);
  legData->AddEntry(hRatioVtxZData, "vtxZ / recommended", "l");
  legData->AddEntry(hRatioNoPileupData, "noPileup / recommended", "l");
  legData->AddEntry(hRatioNoEvSelData, "noEvSel / recommended", "l");
  legData->Draw("same");

  canvDataRatio->SaveAs("hSelectionComparisonDataRatio.pdf");

  TCanvas* canvGlauRatio = new TCanvas("canvGlauRatio", "", 1200, 1000);
  hRatioVtxZGlau->GetXaxis()->SetRangeUser(xmin, xmax);
  hRatioVtxZGlau->GetYaxis()->SetRangeUser(ratioYmin, ratioYmax);
  hRatioVtxZGlau->GetYaxis()->SetTitle("variation / recommended");
  hRatioVtxZGlau->SetTitle("");

  hRatioVtxZGlau->Draw("hist");
  hRatioNoPileupGlau->Draw("hist same");
  hRatioNoEvSelGlau->Draw("hist same");

  TLine* lineGlau = new TLine(xmin, 1, xmax, 1);
  lineGlau->SetLineStyle(2);
  lineGlau->SetLineColor(kGray + 2);
  lineGlau->Draw("same");

  TLegend* legGlau = new TLegend(0.6, 0.7, 0.88, 0.88);
  legGlau->SetBorderSize(0);
  legGlau->SetFillStyle(0);
  legGlau->AddEntry(hRatioVtxZGlau, "vtxZ / recommended", "l");
  legGlau->AddEntry(hRatioNoPileupGlau, "noPileup / recommended", "l");
  legGlau->AddEntry(hRatioNoEvSelGlau, "noEvSel / recommended", "l");
  legGlau->Draw("same");

  canvGlauRatio->SaveAs("hSelectionComparisonGlauRatio.pdf");

  
  syst::GlauberParameters glauBase("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters glauVtxZ("LHC25ae_pass2_extra2/AR_564374_calibration_hFT0M_Collisions_vtxZ_LightIonDef.root");
  syst::GlauberParameters glauNopileUp("LHC25ae_pass2_extra2/AR_564374_calibration_hFT0M_Collisions_nopileup_LightIonDef.root");
  syst::GlauberParameters glauNoSel("LHC25ae_pass2_extra2/AR_564374_calibration_hFT0M_Collisions_nosel_LightIonDef.root");

  ColorManager cm(3);
  TH1F* hBaseNpart = syst::initGlauParHist(glauBase.npart, "hBaseNpart", kBlack);
  TH1F* hVtxZNpart = syst::initGlauParHist(glauVtxZ.npart, "hVtxZNpart", cm.getColor(0));
  TH1F* hNoPileUpNpart = syst::initGlauParHist(glauNopileUp.npart, "hNoPileUpNpart", cm.getColor(1));
  TH1F* hNoSelNpart = syst::initGlauParHist(glauNoSel.npart, "hNoSelNpart", cm.getColor(2));
  TH1F* hBaseNcoll = syst::initGlauParHist(glauBase.ncoll, "hBaseNcoll", kBlack);
  TH1F* hVtxZNcoll = syst::initGlauParHist(glauVtxZ.ncoll, "hVtxZNcoll", cm.getColor(0));
  TH1F* hNoPileUpNcoll = syst::initGlauParHist(glauNopileUp.ncoll, "hNoPileUpNcoll", cm.getColor(1));
  TH1F* hNoSelNcoll = syst::initGlauParHist(glauNoSel.ncoll, "hNoSelNcoll", cm.getColor(2));

  // std::vector<TH1F*> hNpart = { hVtxZNpart, hNoPileUpNpart, hNoSelNpart };
  // std::vector<TH1F*> hNcoll = { hVtxZNcoll, hNoPileUpNcoll, hNoSelNcoll };
  std::vector<TH1F*> hNpart = { hVtxZNpart, hNoPileUpNpart };
  std::vector<TH1F*> hNcoll = { hVtxZNcoll, hNoPileUpNcoll };

  TLegend* leg = new TLegend(0.75, 0.65, 0.95, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hBaseNpart, "Recommended selections", "pl");
  leg->AddEntry(hVtxZNpart, "IsGoodZvtxFT0vsPV", "pl");
  leg->AddEntry(hNoPileUpNpart, "NoSameBunchPileUp", "pl");
  // leg->AddEntry(hNoSelNpart, "No selection", "pl");

  doNpartNCollVsAnchorPoint(hBaseNpart, hNpart, leg, "hNpartComparison.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hBaseNcoll, hNcoll, leg, "hNcollComparison.pdf", "<N_{coll}>");

}