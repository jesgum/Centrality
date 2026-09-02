#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <string>

#include "Includes/colorManager.h"

static std::array<int, 2> normRange = { 2000, 12000 };

void styleHistogram(TH1F* hist)
{
  std::array<int, 2> normBins = { hist->FindBin(normRange[0]), hist->FindBin(normRange[1]) };
  hist->SetLineWidth(2);
  hist->Scale(1. / hist->Integral(normBins[0], normBins[1]));
  hist->GetXaxis()->SetTitle("FT0M Amplitude");
  hist->GetYaxis()->SetTitle("Normalized count");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetRangeUser(0, 31000);
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

void styleCanvas(TCanvas* canv, bool log = true)
{
  if (log) {
    canv->SetLogy();
  }

  canv->SetTicks(1, 1);
  canv->SetLeftMargin(0.13);
  canv->SetBottomMargin(0.13);
  canv->SetRightMargin(0.2);
  canv->SetTopMargin(0.03);
}

struct CentralityStudy {
  CentralityStudy(const char* path, const char* subwagon, int col0, int col1)
  {
    TFile file(path, "READ");
    if (file.IsZombie()) {
      std::cerr << "Could not open file!" << std::endl;
      return;
    }

    auto* histColOriginal = (TH1F*)file.Get(Form("centrality-study%s/hFT0M_Collisions", subwagon));
    histColOriginal->SetDirectory(0);
    file.Close();

    TFile fileBC("../AnalysisResults/LHC25ae_pass2/AR_564374.root", "READ");
    auto* histBCsOriginal = (TH1F*)fileBC.Get("centrality-study/hFT0M_BCs");
    histBCsOriginal->SetDirectory(0);
    fileBC.Close();

    // Create a histogram with range 0–70000
    auto* histBCs = new TH1F(Form("histBCs%s", subwagon), histBCsOriginal->GetTitle(), 70000, 0, 70000);
    auto* histCol = new TH1F(Form("histCols%s", subwagon), histBCsOriginal->GetTitle(), 70000, 0, 70000);

    // Copy the contents
    for (int bin = 1; bin <= 70000; ++bin) {
      histBCs->SetBinContent(bin, histBCsOriginal->GetBinContent(bin));
      histBCs->SetBinError(bin, histBCsOriginal->GetBinError(bin));
      histCol->SetBinContent(bin, histColOriginal->GetBinContent(bin));
      histCol->SetBinError(bin, histColOriginal->GetBinError(bin));
    }

    if (!histCol || !histBCs) {
      std::cout << Form("centrality-study%s/hFT0M_Collisions", subwagon) << std::endl;
      std::cout << Form("centrality-study/hFT0M_BCs") << std::endl;
      throw std::runtime_error(Form("Could not find histograms for \"%s\"", subwagon));
    }

    hCol = (TH1F*)histCol->Clone(Form("hCol%s", subwagon));
    hBCs = (TH1F*)histBCs->Clone(Form("hBCs%s", subwagon));

    // Make the clones independent of the file
    hCol->SetDirectory(nullptr);
    hBCs->SetDirectory(nullptr);

    hBCs->SetLineColor(col0);
    hCol->SetLineColor(col1);

    styleHistogram(hBCs);
    styleHistogram(hCol);

    TH1F* hBCs_copy = (TH1F*)hBCs->Clone(Form("bccopy_%s", subwagon));
    TH1F* hCol_copy = (TH1F*)hCol->Clone(Form("colcopy_%s", subwagon));
    hBCs_copy->Rebin(20);
    hCol_copy->Rebin(20);

    styleHistogram(hBCs_copy);
    styleHistogram(hCol_copy);

    // Create ratio histogram with same color as hCol
    hRatio = (TH1F*)hCol_copy->Clone(Form("hRatio%s", subwagon));
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(hCol_copy, hBCs_copy, 1, 1, "B");
    hRatio->SetLineColor(col1);
    // styleHistogram(hRatio);

    // Zoomed-in copy of the ratio, for a separate canvas with a tighter y-axis range
    hRatio_Zoom = (TH1F*)hRatio->Clone(Form("hRatioZoom%s", subwagon));
    hRatio_Zoom->SetDirectory(nullptr);
  }

  TH1F* hCol = nullptr;
  TH1F* hBCs = nullptr;
  TH1F* hRatio = nullptr;
  TH1F* hRatio_Zoom = nullptr;
};

void drawLightIonSameSelections()
{
  gStyle->SetOptStat(0);
  ColorManager cm(3);
  CentralityStudy base("../AnalysisResults/LHC25ae_pass2/AR_564374.root", "", kBlack, cm.getColor(0));
  CentralityStudy pileup("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374.root", "_rejectpileup", kBlack, cm.getColor(1));
  CentralityStudy vtxZ("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374.root", "_vtxZ", kBlack, cm.getColor(2));
  CentralityStudy nobc("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374.root", "_rejectnobc", kBlack, cm.getColor(0));
  CentralityStudy nobcpileup("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374.root", "_rejectnobc_rejectpileup", kBlack, cm.getColor(1));
  CentralityStudy nobcvtxZ("../AnalysisResults/LHC25ae_pass2_extra2/AR_564374.root", "_rejectnobc_vtxZ", kBlack, cm.getColor(2));

  TLegend* leg = new TLegend(0.8, 0.35, 1.1, 0.65);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(base.hBCs, "BCs", "l");
  leg->AddEntry(base.hCol, "base", "l");
  leg->AddEntry(pileup.hCol, "pileup", "l");
  leg->AddEntry(vtxZ.hCol, "vtxZ", "l");
  // leg->AddEntry(nobc.hCol, "nobc", "l");
  // leg->AddEntry(nobcpileup.hCol, "nobcpileup", "l");
  // leg->AddEntry(nobcvtxZ.hCol, "nobcvtxZ", "l");

  TCanvas* canv_base = new TCanvas("canv_base", "", 1600, 1000);
  styleCanvas(canv_base);
  base.hBCs->Draw("hist");
  base.hCol->Draw("hist same");
  leg->Draw();
  canv_base->SaveAs("h_data_base.pdf");

  TCanvas* canv_nobc = new TCanvas("canv_nobc", "", 1600, 1000);
  styleCanvas(canv_nobc);
  nobc.hBCs->Draw("hist");
  nobc.hCol->Draw("hist same");
  leg->Draw();
  // canv_nobc->SaveAs("h_data_nobc.pdf");

  TCanvas* canv_pileup = new TCanvas("canv_pileup", "", 1600, 1000);
  styleCanvas(canv_pileup);
  pileup.hBCs->Draw("hist");
  pileup.hCol->Draw("hist same");
  leg->Draw();
  canv_pileup->SaveAs("h_data_pileup.pdf");

  TCanvas* canv_nobcpileup = new TCanvas("canv_nobcpileup", "", 1600, 1000);
  styleCanvas(canv_nobcpileup);
  nobcpileup.hBCs->Draw("hist");
  nobcpileup.hCol->Draw("hist same");
  leg->Draw();
  // canv_nobcpileup->SaveAs("h_data_nobcpileup.pdf");

  TCanvas* canv_vtxZ = new TCanvas("canv_vtxZ", "", 1600, 1000);
  styleCanvas(canv_vtxZ);
  vtxZ.hBCs->Draw("hist");
  vtxZ.hCol->Draw("hist same");
  leg->Draw();
  canv_vtxZ->SaveAs("h_data_vtxZ.pdf");

  TCanvas* canv_nobcvtxZ = new TCanvas("canv_nobcvtxZ", "", 1600, 1000);
  styleCanvas(canv_nobcvtxZ);
  nobcvtxZ.hBCs->Draw("hist");
  nobcvtxZ.hCol->Draw("hist same");
  leg->Draw();
  // canv_nobcvtxZ->SaveAs("h_data_nobcvtxZ.pdf");

  TCanvas* canv_ratio = new TCanvas("canv_ratio", "", 1600, 1000);
  styleCanvas(canv_ratio, false);
  // base.hRatio->GetYaxis()->SetRangeUser(0, 2);
  base.hRatio->GetYaxis()->SetTitle("Collisions / BCs");
  base.hRatio->Draw("hist");
  pileup.hRatio->Draw("hist same");
  vtxZ.hRatio->Draw("hist same");
  // nobc.hRatio->Draw("hist same");
  // nobcpileup.hRatio->Draw("hist same");
  // nobcvtxZ.hRatio->Draw("hist same");
  leg->Draw();
  canv_ratio->SaveAs("h_data_ratio.pdf");

  // Zoomed-in version of the ratio plot: same curves, tighter y-axis window
  // around 1 so the differences between base/nobc/pileup/nobcpileup are
  // actually resolvable (the full 0-2 range above squashes them together).
  TCanvas* canv_ratio_zoom = new TCanvas("canv_ratio_zoom", "", 1600, 1000);
  styleCanvas(canv_ratio_zoom, false);
  base.hRatio_Zoom->GetYaxis()->SetRangeUser(0, 1.2);
  base.hRatio_Zoom->GetXaxis()->SetRangeUser(0, 6000);
  base.hRatio_Zoom->GetYaxis()->SetTitle("Collisions / BCs");
  base.hRatio_Zoom->Draw("hist");
  nobc.hRatio_Zoom->Draw("hist same");
  vtxZ.hRatio_Zoom->Draw("hist same");
  // pileup.hRatio_Zoom->Draw("hist same");
  // nobcpileup.hRatio_Zoom->Draw("hist same");
  // nobcvtxZ.hRatio_Zoom->Draw("hist same");
  leg->Draw();
  canv_ratio_zoom->SaveAs("h_data_ratio_zoom.pdf");
}