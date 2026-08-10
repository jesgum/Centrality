#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TFormula.h"

#include <iostream>
#include <vector>
#include <string>

#include "colorManager.h"

struct Histos {
  static constexpr int NParams = 6;
  Histos(int runNumber, const char* key = "", bool drawBc = false)
  {
    const char* pathCol = Form("centrality-study%s/hFT0M_Collisions", key);
    const char* pathBcs = Form("centrality-study%s/hFT0M_BCs", key);
    const char* pathData = Form("../AnalysisResults/LHC25ae_pass2_extra2/AR_%d.root", runNumber);
    const char* pathMc = Form("../AnalysisResults/LHC25h3c_extra2/AR_%d.root", runNumber);
    const char* pathCalib = Form("../Glauber/mc/FT0M_mcCalib_%d.root", runNumber);

    TFile* fileData = new TFile(pathData, "read");
    TFile* fileMc = new TFile(pathMc, "read");
    TFile* fileCalib = new TFile(pathCalib, "read");

    hColData = (TH1D*)fileData->Get(pathCol);
    hBcData = (TH1D*)fileData->Get(pathBcs);
    hColMc = (TH1D*)fileMc->Get(pathCol);
    hBcMc = (TH1D*)fileMc->Get(pathBcs);
    f1scale = (TFormula*)fileCalib->Get("f1scale");

    hColScaledMc = (TH1D*)hColMc->Clone(Form("hColMc_Clone_%d", runNumber));
    hBcScaledMc = (TH1D*)hBcMc->Clone(Form("hBcMc_Clone_%d", runNumber));
    hColScaledMc->SetDirectory(nullptr);
    hBcScaledMc->SetDirectory(nullptr);
    hColScaledMc->Reset();
    hBcScaledMc->Reset();

    hColMc->SetLineStyle(7);
    hBcMc->SetLineStyle(7);
    hColScaledMc->SetLineStyle(7);
    hBcScaledMc->SetLineStyle(7);

    for (int ixpar = 0; ixpar < NParams; ++ixpar) {
      parameters[ixpar] = f1scale->GetParameter(ixpar);
    }

    auto scaleMC = [](float x, const float pars[6]) {
      float core = ((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4];
      if (core < 0.0f) {
        return 0.0f; // low multiplicity, not mapped — core^pars[5] would be NaN
      }
      return std::pow(core, 1.0f / pars[5]);
    };

    auto fillScaledHisto = [this, &scaleMC](TH1* hSrc, TH1* hDst) {
      for (int bin = 1; bin <= hSrc->GetNbinsX(); ++bin) {
        const float amp = hSrc->GetBinCenter(bin);
        const float count = hSrc->GetBinContent(bin);
        const float scaledAmp = scaleMC(amp, parameters);
        const double scaledBin = hSrc->FindBin(scaledAmp);
        hDst->SetBinContent(scaledBin, count);
      }
    };

    fillScaledHisto(hColMc, hColScaledMc);
    fillScaledHisto(hBcMc, hBcScaledMc);

    static constexpr int Rebin = 20;
    static constexpr int Smooth = 1000;
    hColData->Rebin(Rebin);
    hBcData->Rebin(Rebin);
    hColMc->Rebin(Rebin);
    hBcMc->Rebin(Rebin);
    hColScaledMc->Rebin(Rebin);
    hBcScaledMc->Rebin(Rebin);

    static std::array<int, 2> normRange = { 2000, 12500 };
    auto styleHistogram = [this](TH1D* hist, int color) {
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
    };

    styleHistogram(hColData, kBlack);
    styleHistogram(hBcData, kBlack);
    styleHistogram(hColMc, kRed + 1);
    styleHistogram(hBcMc, kRed + 1);
    styleHistogram(hColScaledMc, kBlue + 1);
    styleHistogram(hBcScaledMc, kBlue + 1);

    // --- Ratios: MC/Data (before scaling) and ScaledMC/Data (after scaling) ---
    auto makeRatio = [this, runNumber](TH1D* hNum, TH1D* hDen, const char* label, int color) {
      TH1D* hRatio = (TH1D*)hNum->Clone(Form("%s_%d", label, runNumber));
      hRatio->SetDirectory(nullptr);
      hRatio->Divide(hDen); // hRatio = hNum / hDen, bin-by-bin
      hRatio->SetLineColor(color);
      hRatio->SetLineWidth(2);
      hRatio->GetYaxis()->SetTitle("MC / Data");
      hRatio->SetTitle("");
      return hRatio;
    };

    hColRatioMc = makeRatio(hColMc, hColData, "hColRatioMc", kRed + 1);
    hBcRatioMc = makeRatio(hBcMc, hBcData, "hBcRatioMc", kRed + 1);
    hColRatioScaledMc = makeRatio(hColScaledMc, hColData, "hColRatioScaledMc", kBlue + 1);
    hBcRatioScaledMc = makeRatio(hBcScaledMc, hBcData, "hBcRatioScaledMc", kBlue + 1);
    hCollBcRatioData = makeRatio(hColData, hBcData, "hCollBcRatioData", kBlack);
    hCollBcRatioMc = makeRatio(hColScaledMc, hBcScaledMc, "hCollBcRatioMc", kBlue + 1);


    // Zoomed variants: same content, restricted y-axis range around 1
    auto makeZoomed = [](TH1D* hSrc, const char* suffix) {
      TH1D* hZoom = (TH1D*)hSrc->Clone(Form("%s%s", hSrc->GetName(), suffix));
      hZoom->SetDirectory(nullptr);
      return hZoom;
    };

    hColRatioMcZoomed = makeZoomed(hColRatioMc, "_Zoomed");
    hBcRatioMcZoomed = makeZoomed(hBcRatioMc, "_Zoomed");
    hColRatioScaledMcZoomed = makeZoomed(hColRatioScaledMc, "_Zoomed");
    hBcRatioScaledMcZoomed = makeZoomed(hBcRatioScaledMc, "_Zoomed");
    hCollBcRatioDataZoom = makeZoomed(hCollBcRatioData, "_Zoomed");
    hCollBcRatioMcZoom = makeZoomed(hCollBcRatioMc, "_Zoomed");

    auto styleCanvas = [](TCanvas* canv) {
      canv->SetTicks(1, 1);
      canv->SetLeftMargin(0.13);
      canv->SetBottomMargin(0.13);
      canv->SetRightMargin(0.03);
      canv->SetTopMargin(0.03);
    };

    hColData->Smooth(Smooth);
    hColMc->Smooth(Smooth);
    hColScaledMc->Smooth(Smooth);
    hBcData->Smooth(Smooth);
    hBcMc->Smooth(Smooth);
    hBcScaledMc->Smooth(Smooth);

    TLegend* legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(hColData, "Data", "l");
    legend->AddEntry(hColMc, "MC", "l");
    legend->AddEntry(hColScaledMc, "MC (scaled)", "l");

    TCanvas* canv_coll = new TCanvas(Form("canv_coll%s", key), "", 1200, 1000);
    styleCanvas(canv_coll);
    canv_coll->SetLogy();
    hColData->GetXaxis()->SetRangeUser(0, 31000);
    hColData->Draw("hist");
    hColMc->Draw("hist same");
    hColScaledMc->Draw("hist same");
    legend->Draw();
    canv_coll->SaveAs(Form("h_coll%s.pdf", key));

    if (drawBc) {
      TCanvas* canv_bc = new TCanvas(Form("canv_bc%s", key), "", 1200, 1000);
      styleCanvas(canv_bc);
      canv_bc->SetLogy();
      hBcData->GetXaxis()->SetRangeUser(0, 31000);

      hBcData->Draw("hist");
      hBcMc->Draw("hist same");
      hBcScaledMc->Draw("hist same");
      legend->Draw();
      canv_bc->SaveAs(Form("h_bc%s.pdf", key));
    }
    
    // --- Ratio canvases: before (red) vs after (blue) scaling, overlaid on hline at 1 ---
    auto drawRatioToMc = [&styleCanvas](TH1D* hBefore, TH1D* hAfter, double xmax, const char* fname) {
      TCanvas* c = new TCanvas(Form("canv_ratio_%s", fname), "", 1200, 1000);
      styleCanvas(c);
      hBefore->GetXaxis()->SetRangeUser(0, xmax);
      hBefore->Smooth(Smooth);
      hAfter->Smooth(Smooth);
      hBefore->Draw("hist");
      hAfter->Draw("hist same");
      TLine* line = new TLine(0, 1, xmax, 1);
      line->SetLineStyle(2);
      line->SetLineColor(kGray + 2);
      line->Draw("same");
      TLegend* leg = new TLegend(0.64, 0.79, 0.92, 0.92);
      leg->SetBorderSize(0);
      leg->AddEntry(hBefore, "MC / Data (unscaled)", "l");
      leg->AddEntry(hAfter, "MC / Data (scaled)", "l");
      leg->Draw();
      c->SaveAs(Form("%s.pdf", fname));
    };
    
    hColRatioMc->GetYaxis()->SetRangeUser(0, 1.3);
    hColRatioMcZoomed->GetYaxis()->SetRangeUser(0.7, 1.3);
    drawRatioToMc(hColRatioMc, hColRatioScaledMc, 31000, Form("h_coll_ratio%s", key));
    drawRatioToMc(hColRatioMcZoomed, hColRatioScaledMcZoomed, 1000, Form("h_coll_ratio_zoomed%s", key));

    if (drawBc) {
      hBcRatioMc->GetYaxis()->SetRangeUser(0, 1.3);
      hBcRatioMcZoomed->GetYaxis()->SetRangeUser(0.7, 1.3);
      drawRatioToMc(hBcRatioMc, hBcRatioScaledMc, 31000, Form("h_bc_ratio%s", key));
      drawRatioToMc(hBcRatioMcZoomed, hBcRatioScaledMcZoomed, 1000, Form("h_bc_ratio_zoomed%s", key));
    }

    auto drawRatioCollToBc = [&styleCanvas](TH1D* hColl, TH1D* hBc, double xmax, const char* fname) {
      TCanvas* c = new TCanvas(Form("canv_ratio_%s", fname), "", 1200, 1000);
      styleCanvas(c);
      hColl->GetXaxis()->SetRangeUser(0, xmax);
      hColl->Smooth(Smooth);
      hBc->Smooth(Smooth);
      hColl->Draw("hist");
      hBc->Draw("hist same");
      TLine* line = new TLine(0, 1, xmax, 1);
      line->SetLineStyle(2);
      line->SetLineColor(kGray + 2);
      line->Draw("same");
      TLegend* leg = new TLegend(0.16, 0.79, 0.44, 0.92);
      leg->SetBorderSize(0);
      leg->AddEntry(hColl, "Data", "l");
      leg->AddEntry(hBc, "MC (scaled)", "l");
      leg->Draw();
      c->SaveAs(Form("%s.pdf", fname));
    };

    hCollBcRatioData->GetYaxis()->SetRangeUser(0, 2.3);
    hCollBcRatioDataZoom->GetYaxis()->SetRangeUser(0, 1.15);
    drawRatioCollToBc(hCollBcRatioData, hCollBcRatioMc, 31000, Form("h_ratio_coll_bc%s", key));
    drawRatioCollToBc(hCollBcRatioDataZoom, hCollBcRatioMcZoom, 1000, Form("h_ratio_coll_bc_zoomed%s", key));

    hColData->SetLineColor(kRed + 1);
    hColScaledMc->SetLineColor(kRed + 1);
    hBcData->SetLineColor(kBlue + 1);
    hBcScaledMc->SetLineColor(kBlue + 1);

    TLegend* legend_master = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend_master->SetBorderSize(0);
    legend_master->AddEntry(hColData, "Data Coll", "l");
    legend_master->AddEntry(hBcData, "Data Bc", "l");
    legend_master->AddEntry(hColScaledMc, "Mc Coll (scaled)", "l");
    legend_master->AddEntry(hBcScaledMc, "Mc Bc (scaled)", "l");

    TCanvas* canv_master = new TCanvas(Form("canv_master%s", key), "", 1200, 1000);
    styleCanvas(canv_master);
    canv_master->SetLogy();
    hColData->GetXaxis()->SetRangeUser(0, 35000);
    hColData->GetYaxis()->SetRangeUser(1e-7, 1e-1);
    hColData->Draw("hist");
    hColScaledMc->Draw("hist same");
    hBcData->Draw("hist same");
    hBcScaledMc->Draw("hist same");
    legend_master->Draw();
    canv_master->SaveAs(Form("h_all%s.pdf", key));
  } // Histos

  TFormula* f1scale = nullptr;
  float parameters[NParams] = { 0.0 };

  TH1D* hColData = nullptr;
  TH1D* hBcData = nullptr;
  TH1D* hColMc = nullptr;
  TH1D* hBcMc = nullptr;
  TH1D* hColScaledMc = nullptr;
  TH1D* hBcScaledMc = nullptr;

  // Ratios before scaling (MC / Data)
  TH1D* hColRatioMc = nullptr;
  TH1D* hBcRatioMc = nullptr;
  TH1D* hColRatioMcZoomed = nullptr;
  TH1D* hBcRatioMcZoomed = nullptr;

  // Ratios after scaling (ScaledMC / Data)
  TH1D* hColRatioScaledMc = nullptr;
  TH1D* hBcRatioScaledMc = nullptr;
  TH1D* hColRatioScaledMcZoomed = nullptr;
  TH1D* hBcRatioScaledMcZoomed = nullptr;

  // Ratios coll to BCs
  TH1D* hCollBcRatioData = nullptr;
  TH1D* hCollBcRatioMc = nullptr;
  TH1D* hCollBcRatioDataZoom = nullptr;
  TH1D* hCollBcRatioMcZoom = nullptr;
};

void drawCalibrationMC()
{
  gStyle->SetOptStat(0);
  Histos base(564374, "", true);
  Histos vtxZ(564374, "_vtxZ");
  Histos pileup(564374, "_rejectpileup");
}