#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"

#include <array>
#include <iostream>

#include "centralityHelper.h"

void drawLightIonCollisionSelections()
{
  gStyle->SetOptStat(0);
  const float xmin = 0;
  const float xmax = 50000;
  TFile* fileCalibSel = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564374.root");
  TH1D* hFT0MCalib = (TH1D*)fileCalibSel->Get("centrality-study/hFT0M_Collisions");
  hFT0MCalib->SetLineColor(kBlack);
  hFT0MCalib->GetXaxis()->SetRangeUser(xmin, xmax);

  TFile* fileSystSel = new TFile("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374.root");
  TH1D* hFT0MSyst = (TH1D*)fileSystSel->Get("centrality-study/hFT0M_Collisions");
  hFT0MSyst->SetLineColor(kRed);

  std::array<float, 2> normRange = {2000, 10000};
  TH1D* hPercentiles = new TH1D("hPercentiles", "", 999, 0 , 99);
  for (int bin = 0; bin < hPercentiles->GetNbinsX(); ++bin) {
    const double percentile = hPercentiles->GetXaxis()->GetBinCenter(bin + 1);
    const double valCalib = cent::getBoundaryForPercentile(hFT0MCalib, percentile);
    const double valSyst = cent::getBoundaryForPercentile(hFT0MSyst, percentile);
    const double ratio = valSyst / valCalib;
    if (std::isnan(ratio)) {
      continue;
    }
    hPercentiles->SetBinContent(bin + 1, ratio);
  }

  const double calib = cent::getBoundaryForPercentile(hFT0MCalib, 1);
  const double syst = cent::getBoundaryForPercentile(hFT0MSyst, 1);


  hFT0MCalib->Scale(1. / hFT0MCalib->Integral(normRange[0], normRange[1]));
  hFT0MSyst->Scale(1. / hFT0MSyst->Integral(normRange[0], normRange[1]));

  TLine* lCalib = new TLine(calib, hFT0MCalib->GetMinimum(), calib, hFT0MCalib->GetMaximum());
  TLine* lSyst = new TLine(syst, hFT0MSyst->GetMinimum(), syst, hFT0MSyst->GetMaximum());
  lCalib->SetLineColor(kBlack);
  lCalib->SetLineStyle(7);
  lSyst->SetLineColor(kRed + 1);
  lSyst->SetLineStyle(7);

  TCanvas* canvCompare = new TCanvas("canvCompare", "", 1200, 1000);
  canvCompare->SetLogy();
  hFT0MCalib->Draw("hist");
  hFT0MSyst->Draw("hist same");
  lCalib->Draw();
  lSyst->Draw();
  canvCompare->SaveAs("hSelectionComparison.pdf");

  TCanvas* canvRatio = new TCanvas("canvRatio", "", 1200, 1000);
  hPercentiles->GetYaxis()->SetRangeUser(0.995, 1.18);
  hPercentiles->Draw("");
  canvRatio->SaveAs("hPercentileRatios.pdf");
  
}