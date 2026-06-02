#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"

#include <array>
#include <iostream>

void styleHistogram(TH1F* hist, int color)
{
  std::array<int, 2> normRange = {350, 55000};
  std::array<int, 2> normBins = {hist->FindBin(normRange[0]), hist->FindBin(normRange[1])};
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  hist->Scale(1. / hist->Integral(normBins[0], normBins[1]));
  hist->GetXaxis()->SetTitle("FT0C Amplitude");
  hist->GetYaxis()->SetTitle("Normalized count");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->SetTitle("");

  // Adjust axis range
  hist->GetXaxis()->SetRangeUser(0, 125);
  hist->GetYaxis()->SetRangeUser(1e-7, 5e+0);
}

void styleCanvas(TCanvas* canv)
{
  // canv->SetGrid();
  canv->SetTicks(1, 1);
  canv->SetLeftMargin(0.13);
  canv->SetBottomMargin(0.13);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.03);
}

void drawCollisionVsBC(const char* dataset = "LHC25_pass1", const char* ar = "AR_567905.root")
{
  gStyle->SetOptStat(0);
  const char* filePath = Form("../AnalysisResults/%s/%s", dataset, ar);
  TFile* file = new TFile(filePath, "read");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open file " << filePath << std::endl;
    return;
  }

  TH1F* hCoT0C = (TH1F*)file->Get("centrality-study/hFT0C_Collisions");
  TH1F* hBcT0C = (TH1F*)file->Get("centrality-study/hFT0C_BCs");
  // TH1F* hCoT0C = (TH1F*)file->Get("centrality-study_rejectpileup/hFT0C_Collisions");
  // TH1F* hBcT0C = (TH1F*)file->Get("centrality-study_rejectpileup/hFT0C_BCs");

  if (!hCoT0C || !hBcT0C) {
    std::cerr << "Error: Could not retrieve histograms from file " << filePath << std::endl;
    return;
  }

  styleHistogram(hCoT0C, kRed);
  styleHistogram(hBcT0C, kBlue);

  TLegend* legend = new TLegend(0.6, 0.7, 0.92, 0.92);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0, 0);
  legend->AddEntry(hCoT0C, "Collisions", "l");
  legend->AddEntry(hBcT0C, "BCs", "l");

  TCanvas* canvCollVsBC = new TCanvas("canvCollVsBC", "Collision vs BC", 800, 600);
  styleCanvas(canvCollVsBC);
  canvCollVsBC->SetLogy();
  hCoT0C->Draw("hist ");
  hBcT0C->Draw("hist same");
  legend->Draw();
  canvCollVsBC->SaveAs("hCollVsBC.pdf");
}

// ../AnalysisResults/LHC25_pass1/AR_567905.root
//   /AnalysisResults/LHC25_pass1/AR_567905.root
