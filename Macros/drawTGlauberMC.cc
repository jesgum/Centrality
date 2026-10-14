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

std::vector<TH1F*> doNpartNcoll(TString dataset, TString ar, TString basehistos, const char* name, const std::vector<double>& centBins)
{
  Bool_t fixedMu = true;
  Bool_t fixedK = true;
  Int_t ancesterMode = 2;
  Bool_t printNpartNcoll = true;

  TString mu = (fixedMu) ? "fixedMu" : "freeMu";
  TString k = (fixedK) ? "fixedK" : "freeK";
  TString filePath = Form("../AnalysisResults/%s/%s_calibration_ancestorMode%d_%s_%s_hFT0C_BCs%s.root", dataset.Data(), ar.Data(), ancesterMode, k.Data(), mu.Data(), basehistos.Data());
  TFile* infile = new TFile(filePath.Data(), "read");

  if (!infile || infile->IsZombie()) {
    std::cerr << "Could not find " << filePath << std::endl;
    return {};
  }

  TH2F* h2dNpart = dynamic_cast<TH2F*>(infile->Get("h2dNpart"));
  TH2F* h2dNcoll = dynamic_cast<TH2F*>(infile->Get("h2dNcoll"));

  if (!h2dNpart) {
    std::cerr << "Error: h2dNpart not found or wrong type!" << std::endl;
    return {};
  }

  if (!h2dNcoll) {
    std::cerr << "Error: h2dNcoll not found or wrong type!" << std::endl;
    return {};
  }

  h2dNpart->GetXaxis()->SetTitle("Centrality (%)");
  h2dNcoll->GetXaxis()->SetTitle("Centrality (%)");
  h2dNpart->GetXaxis()->SetTitleSize(0.045);
  h2dNcoll->GetXaxis()->SetTitleSize(0.045);
  h2dNpart->GetXaxis()->SetTitleOffset(1.05);
  h2dNcoll->GetXaxis()->SetTitleOffset(1.05);
  h2dNpart->GetYaxis()->SetTitle("N_{part}");
  h2dNcoll->GetYaxis()->SetTitle("N_{coll}");
  h2dNpart->GetYaxis()->SetTitleSize(0.05);
  h2dNcoll->GetYaxis()->SetTitleSize(0.05);
  h2dNpart->GetYaxis()->SetTitleOffset(1);
  h2dNcoll->GetYaxis()->SetTitleOffset(1);

  const int nCentBins = centBins.size() - 1;
  TH1F* hAvgNpart = new TH1F("hAvgNpart_local", "", nCentBins, centBins.data());
  TH1F* hAvgNcoll = new TH1F("hAvgNcoll_local", "", nCentBins, centBins.data());

  hAvgNpart->SetDirectory(0);
  hAvgNcoll->SetDirectory(0);
  hAvgNpart->SetMarkerStyle(kFullCircle);
  hAvgNpart->SetMarkerColor(kBlack);
  hAvgNcoll->SetMarkerStyle(kFullCircle);
  hAvgNcoll->SetMarkerColor(kBlack);

  if (printNpartNcoll) {
    std::cout << " - - - cent - - - <Npart> - - - <Ncoll> - - - " << std::endl;
  }

  // Each output bin projects the 2D x-bins inside its centrality edges
  for (Int_t icent{0}; icent < nCentBins; ++icent) {
    const double lo = centBins[icent];
    const double hi = centBins[icent + 1];
    const int firstBin = h2dNpart->GetXaxis()->FindBin(lo + 1e-3);
    const int lastBin = h2dNpart->GetXaxis()->FindBin(hi - 1e-3);
    TH1D* hProjectionNpart = dynamic_cast<TH1D*>(h2dNpart->ProjectionY(Form("hProjectionNpart_%d", icent), firstBin, lastBin));
    TH1D* hProjectionNcoll = dynamic_cast<TH1D*>(h2dNcoll->ProjectionY(Form("hProjectionNcoll_%d", icent), firstBin, lastBin));
    hAvgNpart->SetBinContent(icent + 1, hProjectionNpart->GetMean());
    hAvgNpart->SetBinError(icent + 1, hProjectionNpart->GetMeanError());
    hAvgNcoll->SetBinContent(icent + 1, hProjectionNcoll->GetMean());
    hAvgNcoll->SetBinError(icent + 1, hProjectionNcoll->GetMeanError());

    if (printNpartNcoll) {
      std::cout << Form("%g-%g%% | ", lo, hi)
                << hProjectionNpart->GetMean() << " | "
                << hProjectionNcoll->GetMean() << std::endl;
    }
    delete hProjectionNpart; delete hProjectionNcoll;
  }

  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 800);
  canvNpart->SetLogz();
  canvNpart->SetTopMargin(0.03);
  canvNpart->SetLeftMargin(0.1);
  canvNpart->SetRightMargin(0.1);
  canvNpart->SetBottomMargin(0.13);
  h2dNpart->Draw("colz");
  hAvgNpart->Draw("pe same");
  // canvNpart->SaveAs(Form("h2dNpart_%s.pdf", name));

  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 800);
  canvNcoll->SetLogz();
  canvNcoll->SetTopMargin(0.03);
  canvNcoll->SetLeftMargin(0.1);
  canvNcoll->SetRightMargin(0.1);
  canvNcoll->SetBottomMargin(0.13);
  h2dNcoll->Draw("colz");
  hAvgNcoll->Draw("pe same");
  // canvNcoll->SaveAs(Form("h2dNcoll_%s.pdf", name));

  delete canvNpart; delete canvNcoll;
  delete h2dNpart; delete h2dNcoll;
  infile->Close(); delete infile;
  std::cout << std::endl;
  return { hAvgNpart, hAvgNcoll };
}

struct Histos {
  Histos(std::vector<TH1F*> vec, const char* lbl)
  {
    histos = vec;
    label = lbl;
  }
  void setColor(int color)
  {
    for (auto& hist : histos) {
      hist->SetMarkerColor(color);
      hist->SetLineColor(color);
    }
  }

  const char* label;
  std::vector<TH1F*> histos;
};

// Draws <X> vs centrality in a top pad and the ratio to the reference in a bottom pad
void drawWithRatio(int quantity, Histos* reference, const std::vector<Histos*>& variations, bool drawReference,
                   const char* yTitle, const char* ratioTitle, double ratioMin, double ratioMax, const char* outName)
{
  TCanvas* canv = new TCanvas(Form("canv_%s", outName), "", 1000, 900);

  const double padSplit = 0.35;
  TPad* padTop = new TPad(Form("padTop_%s", outName), "", 0, padSplit, 1, 1);
  padTop->SetTopMargin(0.04);
  padTop->SetBottomMargin(0.);
  padTop->SetLeftMargin(0.13);
  padTop->SetRightMargin(0.03);
  padTop->Draw();

  canv->cd();
  TPad* padBottom = new TPad(Form("padBottom_%s", outName), "", 0, 0, 1, padSplit);
  padBottom->SetTopMargin(0.);
  padBottom->SetBottomMargin(0.3);
  padBottom->SetLeftMargin(0.13);
  padBottom->SetRightMargin(0.03);
  padBottom->SetGridy();
  padBottom->Draw();

  // Text sizes scaled so that top and bottom pad labels appear equally large
  const double textSize = 0.045;
  const double scaleTop = 1. / (1. - padSplit);
  const double scaleBottom = 1. / padSplit;

  // Top pad
  padTop->cd();
  TH1F* hFrameTop = (TH1F*)reference->histos[quantity]->Clone(Form("hFrameTop_%s", outName));
  hFrameTop->Reset();
  double yMax = 0.;
  for (auto& hist : variations) {
    yMax = std::max(yMax, hist->histos[quantity]->GetMaximum());
  }
  if (drawReference) {
    yMax = std::max(yMax, reference->histos[quantity]->GetMaximum());
  }
  hFrameTop->GetYaxis()->SetRangeUser(1e-3, 1.1 * yMax);
  hFrameTop->GetYaxis()->SetTitle(yTitle);
  hFrameTop->GetYaxis()->SetTitleSize(textSize * scaleTop);
  hFrameTop->GetYaxis()->SetLabelSize(textSize * scaleTop * 0.85);
  hFrameTop->GetYaxis()->SetTitleOffset(1.4 / scaleTop);
  hFrameTop->GetXaxis()->SetLabelSize(0.);
  hFrameTop->Draw("axis");

  TLegend* legend = new TLegend(0.78, 0.55, 1, 0.93);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0, 0);
  legend->SetTextSize(textSize * scaleTop * 0.8);
  if (drawReference) {
    reference->histos[quantity]->Draw("pe same");
    legend->AddEntry(reference->histos[quantity], reference->label, "pl");
  }
  for (auto& hist : variations) {
    hist->histos[quantity]->Draw("pe same");
    legend->AddEntry(hist->histos[quantity], hist->label, "pl");
  }
  legend->Draw();

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(textSize * scaleTop);
  latex->DrawLatex(0.2, 0.85, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");

  // Bottom pad
  padBottom->cd();
  TH1F* hFrameBottom = (TH1F*)reference->histos[quantity]->Clone(Form("hFrameBottom_%s", outName));
  hFrameBottom->Reset();
  hFrameBottom->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  hFrameBottom->GetYaxis()->SetTitle(ratioTitle);
  hFrameBottom->GetYaxis()->SetTitleSize(textSize * scaleBottom);
  hFrameBottom->GetYaxis()->SetLabelSize(textSize * scaleBottom * 0.85);
  hFrameBottom->GetYaxis()->SetTitleOffset(1.25 / scaleBottom);
  hFrameBottom->GetYaxis()->SetNdivisions(505);
  hFrameBottom->GetYaxis()->CenterTitle();
  hFrameBottom->GetXaxis()->SetTitle("Centrality (%)");
  hFrameBottom->GetXaxis()->SetTitleSize(textSize * scaleBottom);
  hFrameBottom->GetXaxis()->SetLabelSize(textSize * scaleBottom * 0.85);
  hFrameBottom->GetXaxis()->SetTitleOffset(1.);
  hFrameBottom->GetXaxis()->SetTickLength(0.03 * scaleBottom);
  hFrameBottom->Draw("axis");

  TLine* line = new TLine(hFrameBottom->GetXaxis()->GetXmin(), 1, hFrameBottom->GetXaxis()->GetXmax(), 1);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kGray + 1);
  line->Draw();

  for (auto& hist : variations) {
    TH1F* hRatio = (TH1F*)hist->histos[quantity]->Clone(Form("hRatio_%s_%s", outName, hist->label));
    hRatio->Divide(hist->histos[quantity], reference->histos[quantity], 1, 1, "");
    hRatio->Draw("same");
  }

  canv->SaveAs(Form("%s.pdf", outName));
}

void drawTGlauberMC()
{
  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  // Centrality bin edges (%)
  const std::vector<double> centBins = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
  // const std::vector<double> centBins = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
  // std::vector<double> centBins;
  // for (int ii = 0; ii < 101; ++ii) {
  //   centBins.push_back(ii);
  // }

  // Reference (denominator of the ratios)
  const char* dataset = "LHC25_pass1";
  const char* ar = "AR_568117";
  Histos* reference = new Histos(doNpartNcoll(dataset, ar, "", "PbPb25", centBins), "Pb");
  if (reference->histos.empty()) {
    std::cerr << "Error: reference histograms not available!" << std::endl;
    return;
  }
  reference->setColor(kBlack);
  const bool drawReference = true; // draw the reference in the top pad

  // Variations: { basehistos suffix, legend label, enabled }
  struct Variation {
    const char* basehistos;
    const char* label;
    bool enabled;
  };
  const std::vector<Variation> variationConfig = {
    { "_basehistos_Pbstar", "Pb*", true },
    { "_basehistos_Pbpn", "Pbpn", true },
    { "_basehistos_PbHN", "PbHN", true },
    { "_basehistos_Pbrw", "Pbrw", true },
    { "_basehistos_Pbpnrw", "Pbpnrw", true },
  };

  // Colors are assigned by position in the full list, so they stay fixed when toggling entries
  ColorManager cm(variationConfig.size());
  std::vector<Histos*> variations;
  for (size_t i = 0; i < variationConfig.size(); ++i) {
    const auto& var = variationConfig[i];
    if (!var.enabled) {
      continue;
    }
    std::vector<TH1F*> vec = doNpartNcoll(dataset, ar, var.basehistos, "PbPb25", centBins);
    if (vec.empty()) {
      continue;
    }
    Histos* hist = new Histos(vec, var.label);
    hist->setColor(cm.getColor(i));
    variations.push_back(hist);
  }

  drawWithRatio(kNpart, reference, variations, drawReference, "#LT N_{part} #GT", "Ratio to Pb", 0.97, 1.045, "hNpartTGlauberMC");
  drawWithRatio(kNcoll, reference, variations, drawReference, "#LT N_{coll} #GT", "Ratio to Pb", 0.97, 1.09, "hNcollTGlauberMC");
}
