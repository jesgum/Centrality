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
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"

void drawNchNancVsCentrality()
{
  gStyle->SetOptStat(0);

  // Centrality bin edges (%), integers in [0, 100]
  // const std::vector<float> centBins = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
  // const std::vector<float> centBins = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
  std::vector<float> centBins;
  for (int ii = 0; ii < 101; ++ii) {
    centBins.push_back(ii);
  }

  TFile* fileFitToPbPb = new TFile("../AnalysisResults/LHC23_pass5/AR_544122_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb.root", "read");
  TH1F* hGlauParPbPb = (TH1F*)fileFitToPbPb->Get("hGlauberParameters");
  TFile* fileFitToOO = new TFile("../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0C_Collisions.root", "read");
  TH1F* hGlauParOO = (TH1F*)fileFitToOO->Get("hGlauberParameters");

  const float muPbPb = hGlauParPbPb->GetBinContent(1);
  const float muOO = hGlauParOO->GetBinContent(1);

  TLine* lineMuPbPb = new TLine(0, muPbPb, 100, muPbPb);
  lineMuPbPb->SetLineColor(kBlue + 1);
  lineMuPbPb->SetLineStyle(kDashed);
  lineMuPbPb->SetLineWidth(2);
  TLine* lineMuOO = new TLine(0, muOO, 100, muOO);
  lineMuOO->SetLineColor(kRed + 1);
  lineMuOO->SetLineStyle(kDashed);
  lineMuOO->SetLineWidth(2);



  syst::GlauberParameters dataPbPb("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root", false, centBins);
  // syst::GlauberParameters<TH2D> dataOO("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root", false, centBins);
  syst::GlauberParameters<TH2D> dataOO("LHC25ae_pass2/AR_564374_calibration_hFT0C_Collisions_LightIonDef.root", false, centBins);

  TH1F* hNchNancPbPb = (TH1F*)dataPbPb.hNancestor->Clone("hNchNancPbPb");
  hNchNancPbPb->Reset();
  TH1F* hNchNancOO = (TH1F*)dataOO.hNancestor->Clone("hNchNancOO");
  hNchNancOO->Reset();

  const bool useRMS = false;

  // Fill ratio <Nch> / <Nanc> per centrality bin, taken directly from the 2D (centrality vs X) histograms.
  // Each output bin projects the 2D x-bins inside its centrality edges.
  auto fillRatio = [&](TH1* hOut, TH2* h2dMult, TH2* h2dNanc) {
    if (!h2dMult || !h2dNanc) {
      std::cerr << "Error in fillRatio(" << hOut->GetName() << "): h2dMult or h2dNancestor not available!" << std::endl;
      return;
    }
    for (int bin = 1; bin <= hOut->GetNbinsX(); ++bin) {
      const double lo = hOut->GetXaxis()->GetBinLowEdge(bin);
      const double hi = hOut->GetXaxis()->GetBinUpEdge(bin);
      const int firstBin = h2dNanc->GetXaxis()->FindBin(lo + 1e-3);
      const int lastBin = h2dNanc->GetXaxis()->FindBin(hi - 1e-3);

      TH1D* hProjNanc = h2dNanc->ProjectionY(Form("hProjNanc_%s_%d", hOut->GetName(), bin), firstBin, lastBin);
      TH1D* hProjMult = h2dMult->ProjectionY(Form("hProjMult_%s_%d", hOut->GetName(), bin), firstBin, lastBin);

      const double nanc = hProjNanc->GetMean();
      const double mult = hProjMult->GetMean();
      const double nancErr = useRMS ? hProjNanc->GetRMS() : hProjNanc->GetMeanError();
      const double multErr = useRMS ? hProjMult->GetRMS() : hProjMult->GetMeanError();
      delete hProjNanc;
      delete hProjMult;

      if (nanc <= 0 || mult <= 0) {
        continue;
      }
      const double ratio = mult / nanc;
      // relative errors added in quadrature (assumes Nch and Nanc uncorrelated)
      const double ratioErr = ratio * std::sqrt(std::pow(multErr / mult, 2) + std::pow(nancErr / nanc, 2));
      hOut->SetBinContent(bin, ratio);
      hOut->SetBinError(bin, ratioErr);
    }
  };

  fillRatio(hNchNancPbPb, dataPbPb.h2dMult, dataPbPb.h2dNancestor);
  fillRatio(hNchNancOO, dataOO.h2dMult, dataOO.h2dNancestor);

  hNchNancPbPb->SetMarkerSize(1.5);
  hNchNancPbPb->SetMarkerStyle(kFullCircle);
  hNchNancPbPb->SetMarkerColor(kBlue + 1);
  hNchNancPbPb->SetLineColor(kBlue + 1);
  hNchNancOO->SetMarkerSize(1.5);
  hNchNancOO->SetMarkerStyle(kFullCircle);
  hNchNancOO->SetMarkerColor(kRed + 1);
  hNchNancOO->SetLineColor(kRed + 1);

  TCanvas* canvNchNanc = new TCanvas("canvNchNanc", "", 1200, 1000);
  canvNchNanc->SetTopMargin(0.04);
  canvNchNanc->SetLeftMargin(0.15);
  canvNchNanc->SetRightMargin(0.04);
  canvNchNanc->SetBottomMargin(0.15);
  canvNchNanc->SetTicks();
  
  hNchNancPbPb->Scale(1. / muPbPb);
  hNchNancOO->Scale(1. / muOO);
  TH1D* hist = new TH1D("hist", "", 100, 0, 100);
  hist->GetXaxis()->SetTitle("Centrality (%)");
hist->GetYaxis()->SetTitle("1/#it{#mu}_{AA} <#it{N}_{ch}> / <#it{N}_{anc}>");
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.2);
  // hist->GetYaxis()->SetRangeUser(0, 150 / muPbPb);
  hist->GetYaxis()->SetRangeUser(0, 1.9);
  
  TLine* lUnity = new TLine(0, 1, 100, 1);
  lUnity->SetLineColor(kBlack);
  lUnity->SetLineStyle(kDashed);
  lUnity->SetLineWidth(2);
  
  hist->Draw("pe");
  lUnity->Draw();
  hNchNancPbPb->Draw("pe same");
  hNchNancOO->Draw("pe same");
  lineMuPbPb->Draw();
  lineMuOO->Draw();

  const float posX = 0.42;
  const float posY = 0.68 + 0.25 / 2;
  const float width = 0.38;
  const float height = 0.25 / 2;
  TLegend* leg = new TLegend(posX, posY, posX + width, posY + height);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hNchNancPbPb, "Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV", "ple");
  leg->AddEntry(hNchNancOO, "OO #sqrt{#it{s}_{NN}} = 5.36 TeV", "ple");
  // leg->AddEntry(lineMuPbPb, Form("NBD #mu: %.1f", muPbPb), "l");
  // leg->AddEntry(lineMuOO, Form("NBD #mu: %.1f", muOO), "l");
  leg->Draw();



  TLatex* lat = new TLatex();
  // lat->SetTextFont(42);
  lat->SetTextSize(0.07);
  lat->DrawTextNDC(0.18, 0.87, "ALICE");

  canvNchNanc->SaveAs("hNchNancVsCentrality.pdf");
}