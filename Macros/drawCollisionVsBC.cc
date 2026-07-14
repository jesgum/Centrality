#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"

#include <array>
#include <iostream>

static std::array<int, 2> normRange = { 2000, 12000 };

void styleHistogram(TH1F* hist, int color)
{
  std::array<int, 2> normBins = { hist->FindBin(normRange[0]), hist->FindBin(normRange[1]) };
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  hist->Scale(1. / hist->Integral(normBins[0], normBins[1]));
  hist->GetXaxis()->SetTitle("FT0M Amplitude");
  hist->GetYaxis()->SetTitle("Normalized count");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->SetTitle("");
}

void styleRatioHistogram(TH1F* hist)
{
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(2);
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle("FT0M Amplitude");
  hist->GetYaxis()->SetTitle("Coll. / BC");
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.12);
  hist->GetYaxis()->SetTitleSize(0.10);
  hist->GetXaxis()->SetLabelSize(0.10);
  hist->GetYaxis()->SetLabelSize(0.10);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleOffset(0.55);
}

void styleCanvas(TCanvas* canv)
{
  canv->SetTicks(1, 1);
  canv->SetLeftMargin(0.13);
  canv->SetBottomMargin(0.13);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.03);
}

// OO           AR_564374.root
// NeNe         AR_564468.root
// PbPb 2023    AR_545210.root
// PbPb 2024    AR_560089.root
// PbPb 2025    AR_567905.root
void drawCollisionVsBC(const char* dataset = "LHC25ae_pass2", const char* ar = "AR_564374.root", const char* system = "OO")
{
  gStyle->SetOptStat(0);
  const char* filePath = Form("../AnalysisResults/%s/%s", dataset, ar);
  TFile* file = new TFile(filePath, "read");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open file " << filePath << std::endl;
    return;
  }

  // TH1F* hCoT0 = (TH1F*)file->Get("centrality-study_no_rejectpileup/hFT0M_Collisions");
  // TH1F* hBcT0 = (TH1F*)file->Get("centrality-study_no_rejectpileup/hFT0M_BCs");
  TH1F* hCoT0 = (TH1F*)file->Get("centrality-study/hFT0M_Collisions");
  TH1F* hBcT0 = (TH1F*)file->Get("centrality-study/hFT0M_BCs");
  TH1F* hCoT0_Zoom = (TH1F*)hCoT0->Clone("hCoT0_Zoom");
  TH1F* hBcT0_Zoom = (TH1F*)hBcT0->Clone("hBcT0_Zoom");
  const float xmin = 0;
  const float xmax = 46000;
  const float ymax = 5e-2;
  const float ymin = 5e-10;
  const float xmin_zoom = 0;
  const float xmax_zoom = 4000;
  const float ymax_zoom = 5e-2;
  const float ymin_zoom = 5e-7;

  if (!hCoT0 || !hBcT0) {
    std::cerr << "Error: Could not retrieve histograms from file " << filePath << std::endl;
    return;
  }

  styleHistogram(hCoT0, kRed);
  styleHistogram(hBcT0, kBlue);
  styleHistogram(hCoT0_Zoom, kRed);
  styleHistogram(hBcT0_Zoom, kBlue);

  TH1F* hRatio = (TH1F*)hCoT0->Clone("hRatio");
  hRatio->Divide(hBcT0);
  styleRatioHistogram(hRatio);

  TH1F* hRatio_Zoom = (TH1F*)hCoT0_Zoom->Clone("hRatio_Zoom");
  hRatio_Zoom->Divide(hBcT0_Zoom);
  styleRatioHistogram(hRatio_Zoom);

  TLegend* legend = new TLegend(0.6, 0.7, 0.92, 0.92);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0, 0);
  legend->AddEntry(hCoT0, "Collisions", "l");
  legend->AddEntry(hBcT0, "BCs", "l");

  const float splitY = 0.30;
  TCanvas* canvCollVsBC = new TCanvas("canvCollVsBC", "Collision vs BC", 800, 800);
  TPad* padUp = new TPad("padUp", "", 0, splitY, 1, 1);
  padUp->SetLeftMargin(0.13);
  padUp->SetRightMargin(0.03);
  padUp->SetTopMargin(0.04);
  padUp->SetBottomMargin(0.02);
  padUp->SetTicks(1, 1);
  padUp->SetLogy();
  padUp->Draw();
  padUp->cd();
  hCoT0->GetXaxis()->SetLabelSize(0);
  hCoT0->GetXaxis()->SetTitleSize(0);
  hBcT0->GetXaxis()->SetLabelSize(0);
  hBcT0->GetXaxis()->SetTitleSize(0);
  hCoT0->GetYaxis()->SetTitleSize(0.05);
  hCoT0->GetYaxis()->SetLabelSize(0.05);
  hCoT0->GetYaxis()->SetTitleOffset(1.2);
  hCoT0->GetXaxis()->SetRangeUser(xmin, xmax);
  hBcT0->GetXaxis()->SetRangeUser(xmin, xmax);
  hCoT0->GetYaxis()->SetRangeUser(ymin, ymax);
  hCoT0->Draw("hist");
  hBcT0->Draw("hist same");
  legend->Draw();

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.08);
  latex->DrawLatex(0.175, 0.80, Form("%s", system));
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.175, 0.75, Form("Normalization range: [%d, %d]", normRange[0], normRange[1]));

  canvCollVsBC->cd();

  TPad* padLo = new TPad("padLo", "", 0, 0, 1, splitY);
  padLo->SetLeftMargin(0.13);
  padLo->SetRightMargin(0.03);
  padLo->SetTopMargin(0.02);
  padLo->SetBottomMargin(0.3);
  padLo->SetTicks(1, 1);
  padLo->Draw();
  padLo->cd();
  hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
  hRatio->GetYaxis()->SetRangeUser(0.01, 1.21);
  hRatio->Draw("hist");

  TLine* line1 = new TLine(xmin, 1, xmax, 1);
  line1->SetLineStyle(2);
  line1->SetLineColor(kGray + 1);
  line1->Draw();
  canvCollVsBC->SaveAs(Form("hCollVsBC_%s.pdf", system));

  // <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  TCanvas* canvCollVsBC_Zoom = new TCanvas("canvCollVsBC_Zoom", "Collision vs BC (Zoomed)", 800, 800);
  TPad* padUp2 = new TPad("padUp2", "", 0, splitY, 1, 1);
  padUp2->SetLeftMargin(0.13);
  padUp2->SetRightMargin(0.03);
  padUp2->SetTopMargin(0.04);
  padUp2->SetBottomMargin(0.02);
  padUp2->SetTicks(1, 1);
  padUp2->SetLogy();
  padUp2->Draw();
  padUp2->cd();
  hCoT0_Zoom->GetXaxis()->SetLabelSize(0);
  hCoT0_Zoom->GetXaxis()->SetTitleSize(0);
  hBcT0_Zoom->GetXaxis()->SetLabelSize(0);
  hBcT0_Zoom->GetXaxis()->SetTitleSize(0);
  hCoT0_Zoom->GetYaxis()->SetTitleSize(0.05);
  hCoT0_Zoom->GetYaxis()->SetLabelSize(0.05);
  hCoT0_Zoom->GetYaxis()->SetTitleOffset(1.2);
  hCoT0_Zoom->GetXaxis()->SetRangeUser(xmin_zoom, xmax_zoom);
  hBcT0_Zoom->GetXaxis()->SetRangeUser(xmin_zoom, xmax_zoom);
  hCoT0_Zoom->GetYaxis()->SetRangeUser(ymin_zoom, ymax_zoom);
  hCoT0_Zoom->Draw("hist");
  hBcT0_Zoom->Draw("hist same");
  legend->Draw();
  latex->SetTextSize(0.08);
  latex->DrawLatex(0.175, 0.80, Form("%s", system));
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.175, 0.75, Form("Normalization range: [%d, %d]", normRange[0], normRange[1]));
  canvCollVsBC_Zoom->cd();

  TPad* padLo2 = new TPad("padLo2", "", 0, 0, 1, splitY);
  padLo2->SetLeftMargin(0.13);
  padLo2->SetRightMargin(0.03);
  padLo2->SetTopMargin(0.02);
  padLo2->SetBottomMargin(0.3);
  padLo2->SetTicks(1, 1);
  padLo2->Draw();
  padLo2->cd();
  hRatio_Zoom->GetXaxis()->SetRangeUser(xmin_zoom, xmax_zoom);
  hRatio_Zoom->GetYaxis()->SetRangeUser(0, 1.1);
  hRatio_Zoom->Draw("hist");

  TLine* line2 = new TLine(xmin_zoom, 1, xmax_zoom, 1);
  line2->SetLineStyle(2);
  line2->SetLineColor(kGray + 1);
  line2->Draw();
  canvCollVsBC_Zoom->SaveAs(Form("hCollVsBC_Zoom_%s.pdf", system));
}