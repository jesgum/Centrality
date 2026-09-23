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

// Draws one canvas (top: base vs. Trento variations, bottom: ratio to base)
// for a single quantity (Npart or Ncoll). Only the variations are plotted
// against the base -- the default (un-varied) Trento tune is not shown.
void drawTrentoVariationsCompare(TH1F* base, std::vector<TH1F*>& variationHists, const std::vector<std::string>& variationNames,
                                  ColorManager& colorMgr, const char* yTitle, const char* outFile)
{
  base->SetMarkerColor(kBlack);
  base->SetLineColor(kBlack);
  for (size_t i = 0; i < variationHists.size(); ++i) {
    const int color = colorMgr.getColor(i);
    variationHists[i]->SetMarkerColor(color);
    variationHists[i]->SetLineColor(color);
  }

  TCanvas* canv = new TCanvas(Form("canv_%s", outFile), "", 1200, 800);
  canv->SetTopMargin(0.0);
  canv->SetBottomMargin(0.0);

  TPad* padTop = new TPad(Form("padTop_%s", outFile), "", 0, 0.3, 1, 1);
  padTop->SetTicks(1, 1);
  padTop->SetTopMargin(0.04);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.12);
  padTop->SetRightMargin(0.03);
  padTop->Draw();
  padTop->cd();
  base->SetMinimum(0);
  base->GetXaxis()->SetLabelSize(0);
  base->GetYaxis()->SetTitle(yTitle);
  base->GetYaxis()->SetTitleSize(0.05);
  base->GetYaxis()->SetLabelSize(0.045);
  base->GetYaxis()->SetTitleOffset(1.1);
  base->Draw("pe");
  for (auto* h : variationHists) {
    h->Draw("pe same");
  }

  TLegend* leg = new TLegend(0.65, 0.5, 0.87, 0.94);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(base, "Base", "p");
  for (size_t i = 0; i < variationHists.size(); ++i) {
    leg->AddEntry(variationHists[i], Form("Trento %s", variationNames[i].c_str()), "p");
  }
  leg->Draw();

  canv->cd();
  TPad* padBot = new TPad(Form("padBot_%s", outFile), "", 0, 0, 1, 0.3);
  padBot->SetTicks(1, 1);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.12);
  padBot->SetRightMargin(0.03);
  padBot->Draw();
  padBot->cd();

  std::vector<TH1F*> ratios;
  for (size_t i = 0; i < variationHists.size(); ++i) {
    TH1F* hRatio = dynamic_cast<TH1F*>(variationHists[i]->Clone(Form("hRatio_%zu_%s", i, outFile)));
    hRatio->Divide(base);
    ratios.push_back(hRatio);
  }

  TH1F* hFrame = ratios[0];
  hFrame->SetTitle("");
  hFrame->GetYaxis()->SetTitle("Trento / Base");
  hFrame->GetYaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetTitleSize(0.12);
  hFrame->GetYaxis()->SetTitleOffset(0.45);
  hFrame->GetYaxis()->SetLabelSize(0.10);
  hFrame->GetXaxis()->SetTitle("Centrality (%)");
  hFrame->GetXaxis()->SetTitleSize(0.13);
  hFrame->GetXaxis()->SetLabelSize(0.11);

  double ratioMin = 1, ratioMax = 1;
  for (TH1F* hRatio : ratios) {
    ratioMin = std::min(ratioMin, hRatio->GetMinimum());
    ratioMax = std::max(ratioMax, hRatio->GetMaximum());
  }
  const double ratioMargin = 0.1 * (ratioMax - ratioMin);
  hFrame->SetMinimum(ratioMin - ratioMargin);
  hFrame->SetMaximum(ratioMax + ratioMargin);

  hFrame->Draw("pe");
  for (size_t i = 1; i < ratios.size(); ++i) {
    ratios[i]->Draw("pe same");
  }

  TLine* line = new TLine(hFrame->GetXaxis()->GetXmin(), 1, hFrame->GetXaxis()->GetXmax(), 1);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 1);
  line->Draw("same");

  canv->SaveAs(outFile);
}

void drawNpartNcollTrento()
{
  gStyle->SetOptStat(0);
  syst::GlauberParameters base("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.root");

  const std::vector<std::string> variationNames = { "A_d0.0", "A_d0.4", "A_d0.8", "B_d0.0", "B_d0.4", "B_d0.8" };
  // comsolpc has been split into two variants; produce the same comparison
  // plots for both.
  const std::vector<std::string> comsolFolders = { "comsolpc_80", "comsolpc_90" };

  for (const auto& folder : comsolFolders) {
    std::vector<syst::GlauberParameters<TH2D>> variations;
    for (const auto& name : variationNames) {
      variations.emplace_back(Form("../Trento/%s/trentofit_%s_calib.root", folder.c_str(), name.c_str()), true);
    }

    // === Colors: base stays black, the 6 variations get distinct colors
    // from the ColorManager ===
    ColorManager colorMgr(variations.size());

    std::vector<TH1F*> npartHists, ncollHists;
    for (auto& variation : variations) {
      npartHists.push_back(variation.hNpart);
      ncollHists.push_back(variation.hNcoll);
    }

    const TString suffix = folder.substr(folder.find('_') + 1).c_str();
    drawTrentoVariationsCompare(base.hNpart, npartHists, variationNames, colorMgr, "<N_{part}>", Form("hNpartCompareTrentoBase_%s.pdf", suffix.Data()));
    drawTrentoVariationsCompare(base.hNcoll, ncollHists, variationNames, colorMgr, "<N_{coll}>", Form("hNcollCompareTrentoBase_%s.pdf", suffix.Data()));
  }
}
