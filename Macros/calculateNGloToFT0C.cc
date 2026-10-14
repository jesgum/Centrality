#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

#include <iostream>

double getBoundaryForPercentile(TH1* histo, double lPercentileRequested ) {
  //This function returns the boundary for a specific percentile.
  double lReturnValue = 0.0;
  double lPercentile = 100.0 - lPercentileRequested;
  
  const Long_t lNBins = histo->GetNbinsX();
  double lCountDesired = lPercentile * histo->GetEntries()/100;
  Long_t lCount = 0;
  for(Long_t ibin=1;ibin<lNBins;ibin++){
    lCount += histo->GetBinContent(ibin);
    if(lCount >= lCountDesired){
      //Found bin I am looking for!
      double lWidth = histo->GetBinWidth(ibin);
      double lLeftPercentile = 100.*(lCount - histo->GetBinContent(ibin))/histo->GetEntries();
      double lRightPercentile = 100.*lCount / histo->GetEntries();
      
      double lProportion = (lPercentile - lLeftPercentile)/(lRightPercentile-lLeftPercentile);
      
      lReturnValue = histo->GetBinLowEdge(ibin) + lProportion*lWidth;
      break;
    }
  }
  return lReturnValue;
}

void calculateNGloToFT0C()
{
  gStyle->SetOptStat(0);
  TFile* infile = new TFile("../AnalysisResults/LHC23_pass5/AR_544122.root", "root");
  if (!infile) {
    std::cerr << "Could not find infile!" << std::endl;
    return;
  }

  TH1D* hFT0C = (TH1D*)infile->Get("centrality-study/hFT0C_Collisions");
  TH1D* hNGlo = (TH1D*)infile->Get("centrality-study/hNGlobalTracks");
  if (!hFT0C) {
    std::cerr << "Could not find hFT0C!" << std::endl;
    return;
  }
  if (!hNGlo) {
    std::cerr << "Could not find hNGlo!" << std::endl;
    return;
  }

  TCanvas* canv = new TCanvas("canv", "", 800, 600);
  canv->SetTopMargin(0.03);
  canv->SetRightMargin(0.03);
  canv->SetLeftMargin(0.13);
  canv->SetBottomMargin(0.13);
  canv->SetLogy();

  hFT0C->SetTitle("");
  hFT0C->GetXaxis()->SetRangeUser(0, 54000);


  // Find scaling factor such that the 0.008% hNGlo boundary maps to 50000
  const double target = 50000.;
  double maxNGlo = getBoundaryForPercentile(hNGlo, 0.008);
  double scale = target / maxNGlo;
  std::cout << "hNGlo boundary at 0.008%: " << maxNGlo << std::endl;
  std::cout << "Scaling factor NGlo -> FT0C: " << scale << std::endl;

  // Rebuild hNGlo with scaled x-axis
  TAxis* ax = hNGlo->GetXaxis();
  TH1D* hNGloScaled = new TH1D("hNGloScaled", "", ax->GetNbins(), ax->GetXmin() * scale, ax->GetXmax() * scale);
  for (int i = 0; i <= ax->GetNbins() + 1; i++) {
    hNGloScaled->SetBinContent(i, hNGlo->GetBinContent(i));
    hNGloScaled->SetBinError(i, hNGlo->GetBinError(i));
  }
  hNGloScaled->SetLineColor(kRed);

  // Normalize in x range 10000 to 40000
  const double normMin = 10000.;
  const double normMax = 40000.;
  double intFT0C = hFT0C->Integral(hFT0C->FindBin(normMin), hFT0C->FindBin(normMax));
  double intNGlo = hNGloScaled->Integral(hNGloScaled->FindBin(normMin), hNGloScaled->FindBin(normMax));
  hFT0C->Scale(1. / intFT0C, "width");
  hNGloScaled->Scale(1. / intNGlo, "width");

  hFT0C->Draw("hist");
  hNGloScaled->Draw("hist same");

  // Inset with zoomed x-axis
  TPad* inset = new TPad("inset", "", 0.4, 0.6, 0.95, 0.95);
  inset->SetTopMargin(0.03);
  inset->SetRightMargin(0.03);
  inset->SetLeftMargin(0.15);
  inset->SetBottomMargin(0.15);
  inset->SetLogy();
  inset->Draw();
  inset->cd();

  TH1D* hFT0CZoom = (TH1D*)hFT0C->Clone("hFT0CZoom");
  TH1D* hNGloZoom = (TH1D*)hNGloScaled->Clone("hNGloZoom");
  hFT0CZoom->GetXaxis()->SetRangeUser(0, 1000);
  hFT0CZoom->Draw("hist");
  hNGloZoom->Draw("hist same");

  canv->cd();
  canv->SaveAs("hNGloToFT0C.pdf");
}