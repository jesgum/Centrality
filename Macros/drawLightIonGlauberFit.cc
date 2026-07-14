#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <string>  

#include "colorManager.h"

struct GlauberInfo
{
  TFile* file = nullptr;
  GlauberInfo(const char* path)
  {
    file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
      throw std::runtime_error(std::string("Cannot open file: ") + path);
    }

    hData = (TH1F*)file->Get("hData");
    hGlauber = (TH1F*)file->Get("hGlauberFine");

    hInsetData = (TH1F*)hData->Clone(Form("hInsetData_%s", path));
    hInsetGlauber = (TH1F*)hGlauber->Clone(Form("hInsetGlauber_%s", path));
    hRatio = (TH1F*)hGlauber->Clone(Form("hRatio_%s", path));
    hInsetRatio = (TH1F*)hInsetGlauber->Clone(Form("hInsetRatio_%s", path));
    TH1F* hDataClone = (TH1F*)hData->Clone("hDataClone");
    TH1F* hInsetGlauberClone = (TH1F*)hInsetGlauber->Clone("hInsetGlauberClone");

    const int rebin = 100;
    hRatio->Rebin(rebin);
    hDataClone->Rebin(rebin);
    hInsetGlauberClone->Rebin(rebin);

    hRatio->Divide(hDataClone, hInsetGlauberClone, 1, 1, "B");
    hInsetRatio->Divide(hInsetGlauber, hInsetData, 1, 1, "B");
    delete hDataClone; delete hInsetGlauberClone;
  }

  ~GlauberInfo()
  {
    if (file) {
      file->Close(); delete file;
      hData = nullptr;
      hGlauber = nullptr;
    }
  }

  TH1F* hData = nullptr;
  TH1F* hGlauber = nullptr;
  TH1F* hRatio = nullptr;

  TH1F* hInsetData = nullptr;
  TH1F* hInsetGlauber = nullptr;
  TH1F* hInsetRatio = nullptr;
};



void drawLightIonGlauberFit()
{
  const float xmax = 34000;
  const float ymax = 3e+9;
  const float insetXmax = 2000;

  float legX = 0.12;
  const float legY = 0.7;
  const float legWidth = 0.2;
  const float legHeight = 0.2;

  gStyle->SetOptStat(0);
  GlauberInfo bcs("../AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_BCs_Anchor90.root");
  GlauberInfo coll("../AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor90.root");
  // GlauberInfo bcs("../AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_BCs_LightIonDef.root");
  // GlauberInfo coll("../AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");

  TCanvas* canvGlauber = new TCanvas("canvGlauber", "", 1600, 1200);
  canvGlauber->SetTopMargin(0.0);
  canvGlauber->SetBottomMargin(0.0);

  bcs.hData->SetLineColor(kBlack);
  coll.hData->SetLineColor(kBlack);
  bcs.hInsetData->SetLineColor(kBlack);
  coll.hInsetData->SetLineColor(kBlack);

  bcs.hGlauber->SetLineColor(kBlue+1);
  coll.hGlauber->SetLineColor(kRed+1);
  bcs.hRatio->SetLineColor(kBlue+1);
  coll.hRatio->SetLineColor(kRed+1);

  bcs.hInsetGlauber->SetLineColor(kBlue+1);
  coll.hInsetGlauber->SetLineColor(kRed+1);
  bcs.hInsetRatio->SetLineColor(kBlue+1);
  coll.hInsetRatio->SetLineColor(kRed+1);

  bcs.hData->SetMarkerColor(kBlack);
  coll.hData->SetMarkerColor(kBlack);
  bcs.hInsetData->SetMarkerColor(kBlack);
  coll.hInsetData->SetMarkerColor(kBlack);

  bcs.hGlauber->SetMarkerColor(kBlue+1);
  coll.hGlauber->SetMarkerColor(kRed+1);
  bcs.hRatio->SetMarkerColor(kBlue+1);
  coll.hRatio->SetMarkerColor(kRed+1);

  bcs.hInsetGlauber->SetMarkerColor(kBlue+1);
  coll.hInsetGlauber->SetMarkerColor(kRed+1);
  bcs.hInsetRatio->SetMarkerColor(kBlue+1);
  coll.hInsetRatio->SetMarkerColor(kRed+1);
  bcs.hData->GetXaxis()->SetRangeUser(0, xmax);
  coll.hData->GetXaxis()->SetRangeUser(0, xmax);

  TPad* padTop = new TPad("padTop", "", 0, 0.64, 1, 1);
  padTop->SetTicks(1, 1);
  padTop->SetBottomMargin(0.0);
  padTop->SetRightMargin(0.03);
  padTop->SetLeftMargin(0.08);
  padTop->SetTopMargin(0.05);
  padTop->Draw();
  padTop->cd();
  padTop->SetLogy();
  coll.hData->SetMaximum(ymax);
  coll.hData->GetYaxis()->SetTitleOffset(0.7);
  coll.hData->GetXaxis()->SetTitleOffset(1.2);
  coll.hData->GetYaxis()->SetTitle("Count");
  coll.hData->Draw("hist");
  coll.hGlauber->Draw("hist same");

  TLegend* legTop = new TLegend(legX, legY, legX + legWidth, legY + legHeight);
  legTop->SetBorderSize(0);
  legTop->SetFillColorAlpha(0, 0);
  legTop->SetTextFont(42);
  legTop->AddEntry(coll.hData, "Coll data");
  legTop->AddEntry(coll.hGlauber, "Glauber fit");
  legTop->Draw();

  padTop->cd();
  // Coordinates are relative to padTop: (xlow, ylow, xhigh, yhigh)
  TPad* insetTop = new TPad("insetTop", "", 0.53, 0.44, 0.93, 0.88);
  insetTop->SetLogy();
  insetTop->SetTicks(1, 1);
  insetTop->SetRightMargin(0.02);
  insetTop->SetTopMargin(0.02);
  insetTop->SetLeftMargin(0.12);
  insetTop->SetBottomMargin(0.12);
  insetTop->Draw();
  insetTop->cd();

  // coll.hInsetData->GetYaxis()->SetRangeUser(2e+5, 6e+7);
  coll.hInsetData->GetXaxis()->SetRangeUser(0, insetXmax);
  coll.hInsetData->GetYaxis()->SetLabelSize(0.045);
  coll.hInsetData->GetXaxis()->SetLabelSize(0.05);
  coll.hInsetData->GetYaxis()->SetTitleSize(0.06);
  coll.hInsetData->GetXaxis()->SetTitleSize(0.06);
  coll.hInsetData->GetXaxis()->SetTitle("FT0M Amplitude");
  coll.hInsetData->GetYaxis()->SetTitleOffset(0.6);
  coll.hInsetData->GetXaxis()->SetTitleOffset(1.1);
  coll.hInsetData->GetYaxis()->SetTitle("Count");
  coll.hInsetData->Draw("hist");
  coll.hInsetGlauber->Draw("hist same");

  // Return up the chain
  padTop->cd();
  canvGlauber->cd();
  TPad* padMid = new TPad("padMid", "", 0, 0.3, 1, 0.64);
  padMid->SetTicks(1, 1);
  padMid->SetTopMargin(0.0);
  padMid->SetBottomMargin(0.0);
  padMid->SetRightMargin(0.03);
  padMid->SetLeftMargin(0.08);
  padMid->Draw();
  padMid->cd();
  padMid->SetLogy();
  bcs.hData->SetMaximum(ymax);
  bcs.hData->GetYaxis()->SetTitleOffset(0.7);
  bcs.hData->GetXaxis()->SetTitleOffset(1.2);
  bcs.hData->GetYaxis()->SetTitle("Count");
  bcs.hData->Draw("hist");
  bcs.hGlauber->Draw("hist same");

  TLegend* legMid = new TLegend(legX, legY, legX + legWidth, legY + legHeight);
  legMid->SetBorderSize(0);
  legMid->SetTextFont(42);
  legMid->SetFillColorAlpha(0, 0);
  legMid->AddEntry(bcs.hData, "BC data");
  legMid->AddEntry(bcs.hGlauber, "Glauber fit");
  legMid->Draw();

  padMid->cd();
  // Coordinates are relative to padTop: (xlow, ylow, xhigh, yhigh)
  TPad* insetMid = new TPad("insetMid", "", 0.53, 0.46, 0.93, 0.9);
  insetMid->SetLogy();
  insetMid->SetTicks(1, 1);
  insetMid->SetRightMargin(0.02);
  insetMid->SetTopMargin(0.02);
  insetMid->SetLeftMargin(0.12);
  insetMid->SetBottomMargin(0.12);
  insetMid->Draw();
  insetMid->cd();

  bcs.hInsetData->GetYaxis()->SetRangeUser(2e+3, 6e+5);
  bcs.hInsetData->GetXaxis()->SetRangeUser(0, insetXmax);
  bcs.hInsetData->GetYaxis()->SetLabelSize(0.045);
  bcs.hInsetData->GetXaxis()->SetLabelSize(0.05);
  bcs.hInsetData->GetYaxis()->SetTitleSize(0.06);
  bcs.hInsetData->GetXaxis()->SetTitleSize(0.06);
  bcs.hInsetData->GetXaxis()->SetTitle("FT0M Amplitude");
  bcs.hInsetData->GetYaxis()->SetTitleOffset(0.6);
  bcs.hInsetData->GetXaxis()->SetTitleOffset(1.1);
  bcs.hInsetData->GetYaxis()->SetTitle("Count");
  bcs.hInsetData->Draw("hist");
  bcs.hInsetGlauber->Draw("hist same");

  // Return up the chain
  padMid->cd();
  canvGlauber->cd();
  TPad* padBot = new TPad("padBot", "", 0, 0, 1, 0.3);
  padBot->Draw();
  padBot->SetTopMargin(0.0);
  padBot->SetRightMargin(0.03);
  padBot->SetLeftMargin(0.08);
  padBot->SetBottomMargin(0.15);
  padBot->SetTicks(1, 1);
  padBot->Draw();
  padBot->cd();

  coll.hRatio->GetYaxis()->SetRangeUser(0.6, 3.2);
  coll.hRatio->GetXaxis()->SetRangeUser(0, xmax);
  coll.hRatio->GetYaxis()->SetLabelSize(0.045);
  coll.hRatio->GetXaxis()->SetLabelSize(0.05);
  coll.hRatio->GetYaxis()->SetTitleSize(0.06);
  coll.hRatio->GetXaxis()->SetTitleSize(0.06);
  coll.hRatio->GetXaxis()->SetTitle("FT0M Amplitude");
  coll.hRatio->GetYaxis()->SetTitleOffset(0.6);
  coll.hRatio->GetXaxis()->SetTitleOffset(1.1);
  coll.hRatio->GetYaxis()->SetTitle("Fit / Data");

  coll.hRatio->Draw("hist same");
  bcs.hRatio->Draw("hist same");

  legX = 0.75;
  TLegend* legBot = new TLegend(legX, 0.5, legX + legWidth, 0.5 + legHeight);
  legBot->SetBorderSize(0);
  legBot->SetTextFont(42);
  legBot->SetFillColorAlpha(0, 0);
  legBot->AddEntry(coll.hRatio, "Coll ratio", "l");
  legBot->AddEntry(bcs.hRatio, "BC ratio", "l");
  legBot->Draw();

  TLine* line = new TLine(0, 1, xmax, 1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(7);
  line->Draw();

  padBot->cd();
  // Coordinates are relative to padTop: (xlow, ylow, xhigh, yhigh)
  TPad* insetBot = new TPad("insetBot", "", 0.13, 0.46, 0.53, 0.9);
  insetBot->SetTicks(1, 1);
  insetBot->SetRightMargin(0.02);
  insetBot->SetTopMargin(0.02);
  insetBot->SetLeftMargin(0.12);
  insetBot->SetBottomMargin(0.12);
  insetBot->Draw();
  insetBot->cd();

  bcs.hInsetRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
  bcs.hInsetRatio->GetXaxis()->SetRangeUser(0, insetXmax);
  bcs.hInsetRatio->GetYaxis()->SetLabelSize(0.045);
  bcs.hInsetRatio->GetXaxis()->SetLabelSize(0.05);
  bcs.hInsetRatio->GetYaxis()->SetTitleSize(0.06);
  bcs.hInsetRatio->GetXaxis()->SetTitleSize(0.06);
  bcs.hInsetRatio->GetXaxis()->SetTitle("FT0M Amplitude");
  bcs.hInsetRatio->GetYaxis()->SetTitleOffset(0.6);
  bcs.hInsetRatio->GetXaxis()->SetTitleOffset(1.1);
  bcs.hInsetRatio->GetYaxis()->SetTitle("Count");
  bcs.hInsetRatio->Draw("hist");
  coll.hInsetRatio->Draw("hist same");

  TLine* lineInset = new TLine(0, 1, insetXmax, 1);
  lineInset->SetLineColor(kBlack);
  lineInset->SetLineStyle(7);
  lineInset->Draw();

  // Return up the chain
  padBot->cd();
  canvGlauber->cd();
  canvGlauber->SaveAs("hGlauber.pdf");
}