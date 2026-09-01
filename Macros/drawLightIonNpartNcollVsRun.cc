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

#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"

static const std::vector<float> NpartVal = { 26.34, 24.08, 21.03, 16.88, 13.05, 9.88, 7.39, 5.51, 4.15, 3.17, 2.42 };
static const std::vector<float> NpartErr = { 0.59, 0.69, 0.84, 0.89, 0.75, 0.61, 0.48, 0.36, 0.26, 0.23, 0.16 };
static const std::vector<float> NpartRelErr = { 0.022399, 0.028654, 0.039943, 0.052720, 0.057463, 0.061743, 0.064954, 0.065336, 0.062651, 0.072555, 0.066116};
static const std::vector<float> NcollVal = { 42.31, 35.82, 28.77, 20.58, 14.19, 9.64, 6.51, 4.41, 3.00, 2.07, 1.38 };
static const std::vector<float> NcollErr = { 2.61, 2.70, 2.73, 2.35, 1.75, 1.22, 0.82, 0.54, 0.35, 0.26, 0.16 };
static const std::vector<float> NcollRelErr = { 0.061688, 0.075377, 0.094876, 0.114188, 0.123326, 0.126535, 0.125960, 0.122473, 0.116667, 0.125604, 0.115942};

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


void drawLightIonNpartNcollVsRun()
{
  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  syst::GlauberParameters baseNeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564472("LHC25af_pass2_systematics/AR_564472_calibration_hFT0M_Collisions_LightIonDef.root");

  syst::GlauberParameters base("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters base0("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_0_LightIonDef.root");
  syst::GlauberParameters base1("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_1_LightIonDef.root");
  syst::GlauberParameters run564356("LHC25ae_pass2/AR_564356_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564359("LHC25ae_pass2/AR_564359_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564373("LHC25ae_pass2/AR_564373_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564387("LHC25ae_pass2/AR_564387_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564400("LHC25ae_pass2/AR_564400_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564414("LHC25ae_pass2/AR_564414_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564430("LHC25ae_pass2/AR_564430_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564445("LHC25ae_pass2/AR_564445_calibration_hFT0M_Collisions_LightIonDef.root");

  ColorManager cm3(8);
  TH1F* hRun564374Npart = syst::initGlauParHist(base.npart, "hRun564374Npart", kBlack);
  TH1F* hRun564356Npart = syst::initGlauParHist(run564356.npart, "hRun564356Npart", cm3.getColor(0));
  TH1F* hRun564359Npart = syst::initGlauParHist(run564359.npart, "hRun564359Npart", cm3.getColor(1));
  TH1F* hRun564373Npart = syst::initGlauParHist(run564373.npart, "hRun564373Npart", cm3.getColor(2));
  TH1F* hRun564387Npart = syst::initGlauParHist(run564387.npart, "hRun564387Npart", cm3.getColor(3));
  TH1F* hRun564400Npart = syst::initGlauParHist(run564400.npart, "hRun564400Npart", cm3.getColor(4));
  TH1F* hRun564414Npart = syst::initGlauParHist(run564414.npart, "hRun564414Npart", cm3.getColor(5));
  TH1F* hRun564430Npart = syst::initGlauParHist(run564430.npart, "hRun564430Npart", cm3.getColor(6));
  TH1F* hRun564445Npart = syst::initGlauParHist(run564445.npart, "hRun564445Npart", cm3.getColor(7));
  TH1F* hRun564374Ncoll = syst::initGlauParHist(base.ncoll, "hRun564374Ncoll", kBlack);
  TH1F* hRun564356Ncoll = syst::initGlauParHist(run564356.ncoll, "hRun564356Ncoll", cm3.getColor(0));
  TH1F* hRun564359Ncoll = syst::initGlauParHist(run564359.ncoll, "hRun564359Ncoll", cm3.getColor(1));
  TH1F* hRun564373Ncoll = syst::initGlauParHist(run564373.ncoll, "hRun564373Ncoll", cm3.getColor(2));
  TH1F* hRun564387Ncoll = syst::initGlauParHist(run564387.ncoll, "hRun564387Ncoll", cm3.getColor(3));
  TH1F* hRun564400Ncoll = syst::initGlauParHist(run564400.ncoll, "hRun564400Ncoll", cm3.getColor(4));
  TH1F* hRun564414Ncoll = syst::initGlauParHist(run564414.ncoll, "hRun564414Ncoll", cm3.getColor(5));
  TH1F* hRun564430Ncoll = syst::initGlauParHist(run564430.ncoll, "hRun564430Ncoll", cm3.getColor(6));
  TH1F* hRun564445Ncoll = syst::initGlauParHist(run564445.ncoll, "hRun564445Ncoll", cm3.getColor(7));

  TH1F* hBaseNpart2 = syst::initGlauParHist(base.npart, "hBaseNpart", kBlack);
  TH1F* hBaseNcoll2 = syst::initGlauParHist(base.ncoll, "hBaseNcoll", kBlack);
  TH1F* hBaseNpart0 = syst::initGlauParHist(base0.npart, "hBaseNpart", kRed + 1);
  TH1F* hBaseNpart1 = syst::initGlauParHist(base1.npart, "hBaseNpart", kBlue + 1);
  TH1F* hBaseNcoll0 = syst::initGlauParHist(base0.ncoll, "hBaseNcoll", kRed + 1);
  TH1F* hBaseNcoll1 = syst::initGlauParHist(base1.ncoll, "hBaseNcoll", kBlue + 1);

  std::vector<TH1F*> hColNpart = { hRun564374Npart, hRun564356Npart, hRun564359Npart, hRun564373Npart, hRun564387Npart, hRun564400Npart, hRun564414Npart, hRun564430Npart, hRun564445Npart };
  std::vector<TH1F*> hColNcoll = { hRun564374Ncoll, hRun564356Ncoll, hRun564359Ncoll, hRun564373Ncoll, hRun564387Ncoll, hRun564400Ncoll, hRun564414Ncoll, hRun564430Ncoll, hRun564445Ncoll };

  TH1F* hRun564468Npart = syst::initGlauParHist(baseNeNe.npart, "hRun564468Npart", kBlack);
  TH1F* hRun564468Ncoll = syst::initGlauParHist(baseNeNe.ncoll, "hRun564468Ncoll", kBlack);
  TH1F* hRun564472Npart = syst::initGlauParHist(run564472.npart, "hRun564472Npart", kBlue + 2);
  TH1F* hRun564472Ncoll = syst::initGlauParHist(run564472.ncoll, "hRun564472Ncoll", kBlue + 2);
  std::vector<TH1F*> hColNpartNeNe = { hRun564472Npart };
  std::vector<TH1F*> hColNcollNeNe = { hRun564472Ncoll };

  TLegend* legRun = new TLegend(0.75, 0.35, 0.95, 0.9);
  legRun->SetBorderSize(0);
  legRun->SetFillColorAlpha(0, 0);
  legRun->AddEntry(hRun564374Npart, "564374 (Ref)", "pl");
  legRun->AddEntry(hRun564356Npart, "564356", "pl");
  legRun->AddEntry(hRun564359Npart, "564359", "pl");
  legRun->AddEntry(hRun564373Npart, "564373", "pl");
  legRun->AddEntry(hRun564387Npart, "564387", "pl");
  legRun->AddEntry(hRun564400Npart, "564400", "pl");
  legRun->AddEntry(hRun564414Npart, "564414", "pl");
  legRun->AddEntry(hRun564430Npart, "564430", "pl");
  legRun->AddEntry(hRun564445Npart, "564445", "pl");

  TLegend* legRunNeNe = new TLegend(0.75, 0.35, 0.95, 0.9);
  legRunNeNe->SetBorderSize(0);
  legRunNeNe->SetFillColorAlpha(0, 0);
  legRunNeNe->AddEntry(hRun564468Npart, "564468 (Ref)", "pl");
  legRunNeNe->AddEntry(hRun564472Npart, "564472", "pl");


  std::vector<TH1F*> hBaseNpart = {hBaseNpart0, hBaseNpart1};
  std::vector<TH1F*> hBaseNcoll = {hBaseNcoll0, hBaseNcoll1};

  TLegend* legBase = new TLegend(0.75, 0.55, 0.95, 0.9);
  legBase->SetBorderSize(0);
  legBase->SetFillColorAlpha(0, 0);
  legBase->AddEntry(hBaseNpart0, "Mode 0", "pl");
  legBase->AddEntry(hBaseNpart1, "Mode 1", "pl");
  legBase->AddEntry(hBaseNcoll2, "Mode 2 (Ref)", "pl");

  doNpartNCollVsAnchorPoint(hRun564374Npart, hColNpart, legRun, "FiguresSystematics/hNpartOOVsRun.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hRun564374Ncoll, hColNcoll, legRun, "FiguresSystematics/hNcollOOVsRun.pdf", "<N_{coll}>");
  doNpartNCollVsAnchorPoint(hRun564468Npart, hColNpartNeNe, legRunNeNe, "FiguresSystematics/hNpartNeNeVsRun.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hRun564468Ncoll, hColNcollNeNe, legRunNeNe, "FiguresSystematics/hNcollNeNeVsRun.pdf", "<N_{coll}>");
  doNpartNCollVsAnchorPoint(hBaseNpart2, hBaseNpart, legBase, "FiguresSystematics/hNpartVsBase.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hBaseNcoll2, hBaseNcoll, legBase, "FiguresSystematics/hNcollVsBase.pdf", "<N_{coll}>");
  
}
