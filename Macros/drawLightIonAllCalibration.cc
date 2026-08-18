#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <iostream>
#include <vector>
#include <string>

#include "colorManager.h"

void styleTH1(TH1F* hist, int color)
{
  hist->SetLineColor(color);
  hist->GetXaxis()->SetRangeUser(0, 14000);
  hist->GetXaxis()->SetTitle("FT0M Amplitude");
  hist->GetYaxis()->SetTitle("Centrality precentile");
}

void styleTGraph(TGraph* gr, int color)
{
  gr->SetMarkerStyle(kFullCircle);
  gr->SetMarkerColor(color);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Centrality percentile");
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitle("FT0M Amplitude");
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetRangeUser(0, 100);
}

void drawLightIonAllCalibration()
{
  gStyle->SetOptStat(0);
  const int anchorPoint = 70;
  TFile* fCollCalib564356 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564356_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564359 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564359_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564373 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564373_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564374 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564387 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564387_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564400 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564400_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564414 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564414_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564430 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564430_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fCollCalib564445 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564445_calibration_hFT0M_Collisions_Anchor70.root", "read");
  TFile* fBcCalib564445 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564445_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564430 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564430_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564414 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564414_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564400 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564400_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564387 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564387_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564374 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564373 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564373_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564359 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564359_calibration_hFT0M_BCs_Anchor70.root", "read");
  TFile* fBcCalib564356 = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564356_calibration_hFT0M_BCs_Anchor70.root", "read");

  TH1F* hCollCalib564445 = (TH1F*)fCollCalib564356->Get("hCalib");
  TH1F* hCollCalib564430 = (TH1F*)fCollCalib564359->Get("hCalib");
  TH1F* hCollCalib564414 = (TH1F*)fCollCalib564373->Get("hCalib");
  TH1F* hCollCalib564400 = (TH1F*)fCollCalib564374->Get("hCalib");
  TH1F* hCollCalib564387 = (TH1F*)fCollCalib564387->Get("hCalib");
  TH1F* hCollCalib564374 = (TH1F*)fCollCalib564400->Get("hCalib");
  TH1F* hCollCalib564373 = (TH1F*)fCollCalib564414->Get("hCalib");
  TH1F* hCollCalib564359 = (TH1F*)fCollCalib564430->Get("hCalib");
  TH1F* hCollCalib564356 = (TH1F*)fCollCalib564445->Get("hCalib");
  TH1F* hBcCalib564445 = (TH1F*)fBcCalib564445->Get("hCalib");
  TH1F* hBcCalib564430 = (TH1F*)fBcCalib564430->Get("hCalib");
  TH1F* hBcCalib564414 = (TH1F*)fBcCalib564414->Get("hCalib");
  TH1F* hBcCalib564400 = (TH1F*)fBcCalib564400->Get("hCalib");
  TH1F* hBcCalib564387 = (TH1F*)fBcCalib564387->Get("hCalib");
  TH1F* hBcCalib564374 = (TH1F*)fBcCalib564374->Get("hCalib");
  TH1F* hBcCalib564373 = (TH1F*)fBcCalib564373->Get("hCalib");
  TH1F* hBcCalib564359 = (TH1F*)fBcCalib564359->Get("hCalib");
  TH1F* hBcCalib564356 = (TH1F*)fBcCalib564356->Get("hCalib");

  ColorManager cm(9);
  hCollCalib564445->SetLineColor(0);
  hCollCalib564430->SetLineColor(1);
  hCollCalib564414->SetLineColor(2);
  hCollCalib564400->SetLineColor(3);
  hCollCalib564387->SetLineColor(4);
  hCollCalib564374->SetLineColor(5);
  hCollCalib564373->SetLineColor(6);
  hCollCalib564359->SetLineColor(7);
  hCollCalib564356->SetLineColor(8);
  hBcCalib564445->SetLineColor(0);
  hBcCalib564430->SetLineColor(1);
  hBcCalib564414->SetLineColor(2);
  hBcCalib564400->SetLineColor(3);
  hBcCalib564387->SetLineColor(4);
  hBcCalib564374->SetLineColor(5);
  hBcCalib564373->SetLineColor(6);
  hBcCalib564359->SetLineColor(7);
  hBcCalib564356->SetLineColor(8);

  TCanvas* canvColl = new TCanvas("canvColl", "", 1200, 1000);
  hCollCalib564445->GetXaxis()->SetRangeUser(0, 10000);
  hCollCalib564445->Draw("hist");
  hCollCalib564430->Draw("hist same");
  hCollCalib564414->Draw("hist same");
  hCollCalib564400->Draw("hist same");
  hCollCalib564387->Draw("hist same");
  hCollCalib564374->Draw("hist same");
  hCollCalib564373->Draw("hist same");
  hCollCalib564359->Draw("hist same");
  hCollCalib564356->Draw("hist same");
  canvColl->SaveAs("hCalibCrossCheckColl.pdf");

  TCanvas* canvBc = new TCanvas("canvBc", "", 1200, 1000);
  hBcCalib564445->GetXaxis()->SetRangeUser(0, 10000);
  hBcCalib564445->Draw("hist");
  hBcCalib564430->Draw("hist same");
  hBcCalib564414->Draw("hist same");
  hBcCalib564400->Draw("hist same");
  hBcCalib564387->Draw("hist same");
  hBcCalib564374->Draw("hist same");
  hBcCalib564373->Draw("hist same");
  hBcCalib564359->Draw("hist same");
  hBcCalib564356->Draw("hist same");
  canvBc->SaveAs("hCalibCrossCheckBc.pdf");

}