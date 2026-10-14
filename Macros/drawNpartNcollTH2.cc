#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"

void drawMarkerForPercentile(TCanvas* canv, TH2D* h2d)
{
  canv->cd();
  TH1D* hAvg = new TH1D(Form("hAvg_%s", h2d->GetName()), "", 100, 0, 100);
  for (int ii = 0; ii < h2d->GetNbinsX(); ++ii) {
    TH1D* hProj = (TH1D*)h2d->ProjectionY(Form("hProj_%d", ii + 1), ii + 1, ii +1);
    hAvg->SetBinContent(ii + 1, hProj->GetMean());
    hAvg->SetBinError(ii + 1, hProj->GetMeanError());
    hAvg->SetMarkerColor(kBlack);
    hAvg->SetLineColor(kBlack);
    hAvg->SetMarkerStyle(kFullCircle);
    hAvg->Draw("same");
  }
}

void drawNpartNcollTH2()
{
  gStyle->SetOptStat(0);
  TFile* infile = new TFile("..//AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root", "read");
  TH2D* h2dNpart = (TH2D*)infile->Get("h2dNpart");
  TH2D* h2dNcoll = (TH2D*)infile->Get("h2dNcoll");

  h2dNcoll->RebinY(3);
  h2dNpart->SetMinimum(1e-10);
  h2dNcoll->SetMinimum(1e-10);
  h2dNpart->GetYaxis()->SetRangeUser(0, 35);
  h2dNcoll->GetYaxis()->SetRangeUser(0, 100);
  h2dNpart->GetYaxis()->SetTitle("#it{N}_{part}");
  h2dNcoll->GetYaxis()->SetTitle("#it{N}_{coll}");
  h2dNpart->GetXaxis()->SetTitle("FT0M percentile");
  h2dNcoll->GetXaxis()->SetTitle("FT0M percentile");
  h2dNpart->GetYaxis()->SetTitleSize(0.06);
  h2dNcoll->GetYaxis()->SetTitleSize(0.06);
  h2dNpart->GetYaxis()->SetTitleOffset(1.05);
  h2dNcoll->GetYaxis()->SetTitleOffset(1.05);
  h2dNpart->GetXaxis()->SetTitleSize(0.05);
  h2dNcoll->GetXaxis()->SetTitleSize(0.05);

  TCanvas* canvNpart = new TCanvas("canvNpart", "", 800, 800);
  canvNpart->SetLeftMargin(0.13);
  canvNpart->SetRightMargin(0.13);
  canvNpart->SetBottomMargin(0.13);
  canvNpart->SetTopMargin(0.13);
  canvNpart->SetLogz();
  h2dNpart->Draw("colz");
  drawMarkerForPercentile(canvNpart, h2dNpart);
  canvNpart->SaveAs("h2dNpartOO.pdf");

  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 800, 800);
  canvNcoll->SetLeftMargin(0.13);
  canvNcoll->SetRightMargin(0.13);
  canvNcoll->SetBottomMargin(0.13);
  canvNcoll->SetTopMargin(0.13);
  canvNcoll->SetLogz();
  h2dNcoll->Draw("colz");
  drawMarkerForPercentile(canvNcoll, h2dNcoll);
  canvNcoll->SaveAs("h2dNcollOO.pdf");
}