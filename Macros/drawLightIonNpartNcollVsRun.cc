#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "colorManager.h"
#include "centralityHelper.h"

static const std::vector<float> CentBins = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
static const std::vector<float> NpartVal = { 26.34, 24.08, 21.03, 16.88, 13.05, 9.88, 7.39, 5.51, 4.15, 3.17, 2.42 };
static const std::vector<float> NpartErr = { 0.59, 0.69, 0.84, 0.89, 0.75, 0.61, 0.48, 0.36, 0.26, 0.23, 0.16 };
static const std::vector<float> NpartRelErr = { 0.022399, 0.028654, 0.039943, 0.052720, 0.057463, 0.061743, 0.064954, 0.065336, 0.062651, 0.072555, 0.066116};
static const std::vector<float> NcollVal = { 42.31, 35.82, 28.77, 20.58, 14.19, 9.64, 6.51, 4.41, 3.00, 2.07, 1.38 };
static const std::vector<float> NcollErr = { 2.61, 2.70, 2.73, 2.35, 1.75, 1.22, 0.82, 0.54, 0.35, 0.26, 0.16 };
static const std::vector<float> NcollRelErr = { 0.061688, 0.075377, 0.094876, 0.114188, 0.123326, 0.126535, 0.125960, 0.122473, 0.116667, 0.125604, 0.115942};
static const std::vector<double> centBinsD(CentBins.begin(), CentBins.end());
static TH1F* hNpartNcollBase = new TH1F("hhNpartNcollBase", "", static_cast<int>(centBinsD.size()) - 1, centBinsD.data());

struct GlauberParameters {
  TH2F* h2dNpart = nullptr;
  TH2F* h2dNcoll = nullptr;
  std::vector<float> npart, ncoll;
  GlauberParameters(const char* path);
};

GlauberParameters::GlauberParameters(const char* path)
{
  TFile* infile = new TFile(Form("../AnalysisResults/%s", path), "read");
  if (!infile || infile->IsZombie()) {
    std::cerr << "Error: file not found!" << std::endl;
  }

  h2dNpart = dynamic_cast<TH2F*>(infile->Get("h2dNpart"));
  h2dNcoll = dynamic_cast<TH2F*>(infile->Get("h2dNcoll"));

  if (!h2dNpart) {
    std::cerr << "Error: h2dNpart not found!" << std::endl;
  }

  if (!h2dNcoll) {
    std::cerr << "Error: h2dNcoll not found!" << std::endl;
  }

  h2dNpart->SetDirectory(0);
  h2dNcoll->SetDirectory(0);
  h2dNpart->GetYaxis()->SetRangeUser(0, 40);
  h2dNcoll->GetYaxis()->SetRangeUser(0, 100);

  infile->Close();
  delete infile;

  for (int icent{ 0 }; icent < h2dNcoll->GetNbinsX(); ++icent) {
    TH1D* hProjectionNpart = dynamic_cast<TH1D*>(h2dNpart->ProjectionY(Form("hProjectionNpart_%d", icent + 1), icent + 1, icent + 1));
    TH1D* hProjectionNcoll = dynamic_cast<TH1D*>(h2dNcoll->ProjectionY(Form("hProjectionNcoll_%d", icent + 1), icent + 1, icent + 1));
    npart.push_back(hProjectionNpart->GetMean());
    ncoll.push_back(hProjectionNcoll->GetMean());
    delete hProjectionNpart;
    delete hProjectionNcoll;
  }

  // Rebin npart and ncoll according to CentBins (formerly the static rebin() method)
  auto rebinInPlace = [](std::vector<float>& vals, const std::vector<float>& bins) {
    std::vector<float> weights(vals.size(), 1);
    std::vector<float> result;

    for (int ii = 0; ii < static_cast<int>(bins.size()) - 1; ++ii) {
      const int thisBin = bins[ii];
      const int nextBin = bins[ii + 1];

      std::vector<float> theseValues, theseWeights;
      for (int jj = thisBin; jj < nextBin; ++jj) {
        theseValues.push_back(vals[jj]);
        theseWeights.push_back(weights[jj]);
      }

      result.push_back(cent::weightedAverage(theseValues, theseWeights));
    }

    vals = result;
  };

  rebinInPlace(npart, CentBins);
  rebinInPlace(ncoll, CentBins);
}

TH1F* initGlauParHist(const std::vector<float>& vals, const char* title, const int color)
{
  TH1F* hist = dynamic_cast<TH1F*>(hNpartNcollBase->Clone(title));
  hist->SetMarkerStyle(kFullCircle);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin) {
    hist->SetBinContent(ibin + 1, vals[ibin]);
    hist->SetBinError(ibin + 1, 1e-4);
  }
  return hist;
}

void doGlauParQA(TH2F* h2d, TH1F* h1d, const char* name)
{
  TCanvas* canv = new TCanvas(Form("canv%s", name), "", 1200, 1000);
  canv->SetLogz();
  h2d->Draw("colz");
  h1d->Draw("same pe");
  canv->SaveAs(Form("qa%s.pdf", name));
  delete canv;
}

std::vector<std::vector<float>> computeSystematics(GlauberParameters base, const std::vector<GlauberParameters>& variations)
{
  std::vector<std::vector<float>> variationsNpart;
  std::vector<std::vector<float>> variationsNcoll;
  for (const auto& variation : variations) {
    variationsNpart.push_back(cent::calculateDifference(variation.npart, base.npart));
    variationsNcoll.push_back(cent::calculateDifference(variation.ncoll, base.ncoll));
  }

  return { cent::averageOverVariations(variationsNpart), cent::averageOverVariations(variationsNcoll) };
}

std::vector<std::vector<float>> computeRelativeSystematics(GlauberParameters base, const std::vector<GlauberParameters>& variations)
{
  std::vector<std::vector<float>> variationsNpart;
  std::vector<std::vector<float>> variationsNcoll;
  for (const auto& variation : variations) {
    variationsNpart.push_back(cent::calculateRelativeDifference(variation.npart, base.npart));
    variationsNcoll.push_back(cent::calculateRelativeDifference(variation.ncoll, base.ncoll));
  }

  return { cent::averageOverVariations(variationsNpart), cent::averageOverVariations(variationsNcoll) };
}

template <typename... Rest>
inline std::vector<std::vector<float>> combineSystematicsInQuadrature(const std::vector<std::vector<float>>& first, const Rest&... rest)
{
  std::vector<std::vector<std::vector<float>>> allSystematics{ first, rest... };

  const size_t nVars = first.size();
  std::vector<std::vector<float>> result(nVars);

  for (size_t iVar = 0; iVar < nVars; ++iVar) {
    const size_t nBins = first[iVar].size();
    result[iVar].assign(nBins, float{});

    for (const auto& syst : allSystematics) {
      if (syst[iVar].size() != nBins) {
        std::cerr << "Error in combineSystematicsInQuadrature: inconsistent bin count!" << std::endl;
        return {};
      }
      for (size_t ibin = 0; ibin < nBins; ++ibin) {
        result[iVar][ibin] += syst[iVar][ibin] * syst[iVar][ibin];
      }
    }

    for (auto& val : result[iVar]) {
      val = std::sqrt(val);
    }
  }

  return result;
}

void styleCanvas(TCanvas* canv)
{
  canv->SetTicks(1, 1);
  canv->SetTopMargin(0.03);
  canv->SetLeftMargin(0.14);
  canv->SetRightMargin(0.18);
  canv->SetBottomMargin(0.13);
}

void styleHist(TH1F* hist, const char* titleYaxis)
{
  hist->SetMinimum(0);

  hist->GetXaxis()->SetTitle("Centrality (%)");
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.1);

  hist->GetYaxis()->SetTitle(titleYaxis);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.45);
}

void doSystematics()
{
  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  GlauberParameters base("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564356("LHC25ae_pass2/AR_564356_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564359("LHC25ae_pass2/AR_564359_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564373("LHC25ae_pass2/AR_564373_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564387("LHC25ae_pass2/AR_564387_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564400("LHC25ae_pass2/AR_564400_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564414("LHC25ae_pass2/AR_564414_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564430("LHC25ae_pass2/AR_564430_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564445("LHC25ae_pass2/AR_564445_calibration_hFT0M_Collisions_LightIonDef.root");

  ColorManager cm3(8);
  TH1F* hRun564374Npart = initGlauParHist(base.npart, "hRun564374Npart", kBlack);
  TH1F* hRun564356Npart = initGlauParHist(run564356.npart, "hRun564356Npart", cm3.getColor(0));
  TH1F* hRun564359Npart = initGlauParHist(run564359.npart, "hRun564359Npart", cm3.getColor(1));
  TH1F* hRun564373Npart = initGlauParHist(run564373.npart, "hRun564373Npart", cm3.getColor(2));
  TH1F* hRun564387Npart = initGlauParHist(run564387.npart, "hRun564387Npart", cm3.getColor(3));
  TH1F* hRun564400Npart = initGlauParHist(run564400.npart, "hRun564400Npart", cm3.getColor(4));
  TH1F* hRun564414Npart = initGlauParHist(run564414.npart, "hRun564414Npart", cm3.getColor(5));
  TH1F* hRun564430Npart = initGlauParHist(run564430.npart, "hRun564430Npart", cm3.getColor(6));
  TH1F* hRun564445Npart = initGlauParHist(run564445.npart, "hRun564445Npart", cm3.getColor(7));

  TLegend* legRun = new TLegend(0.75, 0.35, 0.95, 0.9);
  legRun->SetBorderSize(0);
  legRun->SetFillColorAlpha(0, 0);
  legRun->AddEntry(hRun564374Npart, "564374 (Ref)", "pl");
  legRun->AddEntry(hRun564356Npart, "564356", "pl");
  legRun->AddEntry(hRun564359Npart, "564359", "pl");
  legRun->AddEntry(hRun564373Npart, "564373", "pl");
  legRun->AddEntry(hRun564387Npart, "564387", "pl");
  legRun->AddEntry(hRun564400Npart, "564400", "pl");
  legRun->AddEntry(hRun564414Npart, "564414", "pl");
  legRun->AddEntry(hRun564430Npart, "564430", "pl");
  legRun->AddEntry(hRun564445Npart, "564445", "pl");

  TCanvas* canvNpartVsRun = new TCanvas("canvNpartVsRun", "", 1200, 800);
  canvNpartVsRun->SetTopMargin(0.0);
  canvNpartVsRun->SetBottomMargin(0.0);

  TPad* padNpartTopVsRun = new TPad("padNpartTopVsRun", "", 0, 0.3, 1, 1);
  padNpartTopVsRun->SetTicks(1, 1);
  padNpartTopVsRun->SetTopMargin(0.04);
  padNpartTopVsRun->SetBottomMargin(0.02);
  padNpartTopVsRun->SetLeftMargin(0.1);
  padNpartTopVsRun->SetRightMargin(0.03);
  padNpartTopVsRun->Draw();
  padNpartTopVsRun->cd();
  hRun564374Npart->GetXaxis()->SetLabelSize(0);
  hRun564374Npart->GetYaxis()->SetLabelSize(0.05);
  hRun564374Npart->GetYaxis()->SetTitle("<N_{part}>");
  hRun564374Npart->GetYaxis()->SetTitleSize(0.07);
  hRun564374Npart->GetYaxis()->SetTitleOffset(0.6);
  hRun564374Npart->Draw("pe");
  hRun564356Npart->Draw("same pe");
  hRun564359Npart->Draw("same pe");
  hRun564373Npart->Draw("same pe");
  hRun564387Npart->Draw("same pe");
  hRun564400Npart->Draw("same pe");
  hRun564414Npart->Draw("same pe");
  hRun564430Npart->Draw("same pe");
  hRun564445Npart->Draw("same pe");
  legRun->Draw();

  canvNpartVsRun->cd();
  TPad* padNpartBotVsRun = new TPad("padNpartBotVsRun", "", 0, 0, 1, 0.3);
  padNpartBotVsRun->SetTicks(1, 1);
  padNpartBotVsRun->SetTopMargin(0.02);
  padNpartBotVsRun->SetBottomMargin(0.35);
  padNpartBotVsRun->SetLeftMargin(0.1);
  padNpartBotVsRun->SetRightMargin(0.03);
  padNpartBotVsRun->Draw();
  padNpartBotVsRun->cd();

  std::vector<std::pair<TH1F*, TH1F*>> runPairsNpartVsRun = {
    {hRun564356Npart, hRun564374Npart},
    {hRun564359Npart, hRun564374Npart},
    {hRun564373Npart, hRun564374Npart},
    {hRun564387Npart, hRun564374Npart},
    {hRun564400Npart, hRun564374Npart},
    {hRun564414Npart, hRun564374Npart},
    {hRun564430Npart, hRun564374Npart},
    {hRun564445Npart, hRun564374Npart}
  };

  std::vector<TH1F*> hRatiosNpartVsRun;
  for (int i = 0; i < runPairsNpartVsRun.size(); ++i) {
    TH1F* hRatio = (TH1F*)runPairsNpartVsRun[i].first->Clone(Form("hRatioNpart_%d", i));
    hRatio->Divide(hRatio, runPairsNpartVsRun[i].second, 1, 1, "B");
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Run / Ref run");
    hRatio->GetXaxis()->SetTitle("Centrality (%)");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->GetYaxis()->SetRangeUser(0.971, 1.014);
    hRatiosNpartVsRun.push_back(hRatio);
  }
  hRatiosNpartVsRun[0]->Draw("pe");
  for (int i = 1; i < runPairsNpartVsRun.size(); ++i) hRatiosNpartVsRun[i]->Draw("pe same");

  canvNpartVsRun->SaveAs("hNpartOOVsRun.pdf");

  TH1F* hRun564374Ncoll = initGlauParHist(base.ncoll, "hRun564374Ncoll", kBlack);
  TH1F* hRun564356Ncoll = initGlauParHist(run564356.ncoll, "hRun564356Ncoll", cm3.getColor(0));
  TH1F* hRun564359Ncoll = initGlauParHist(run564359.ncoll, "hRun564359Ncoll", cm3.getColor(1));
  TH1F* hRun564373Ncoll = initGlauParHist(run564373.ncoll, "hRun564373Ncoll", cm3.getColor(2));
  TH1F* hRun564387Ncoll = initGlauParHist(run564387.ncoll, "hRun564387Ncoll", cm3.getColor(3));
  TH1F* hRun564400Ncoll = initGlauParHist(run564400.ncoll, "hRun564400Ncoll", cm3.getColor(4));
  TH1F* hRun564414Ncoll = initGlauParHist(run564414.ncoll, "hRun564414Ncoll", cm3.getColor(5));
  TH1F* hRun564430Ncoll = initGlauParHist(run564430.ncoll, "hRun564430Ncoll", cm3.getColor(6));
  TH1F* hRun564445Ncoll = initGlauParHist(run564445.ncoll, "hRun564445Ncoll", cm3.getColor(7));

  TCanvas* canvNcollVsRun = new TCanvas("canvNcollVsRun", "", 1200, 800);
  canvNcollVsRun->SetTopMargin(0.0);
  canvNcollVsRun->SetBottomMargin(0.0);

  TPad* padNcollTopVsRun = new TPad("padNcollTopVsRun", "", 0, 0.3, 1, 1);
  padNcollTopVsRun->SetTicks(1, 1);
  padNcollTopVsRun->SetTopMargin(0.04);
  padNcollTopVsRun->SetBottomMargin(0.02);
  padNcollTopVsRun->SetLeftMargin(0.1);
  padNcollTopVsRun->SetRightMargin(0.03);
  padNcollTopVsRun->Draw();
  padNcollTopVsRun->cd();
  hRun564374Ncoll->GetXaxis()->SetLabelSize(0);
  hRun564374Ncoll->GetYaxis()->SetLabelSize(0.05);
  hRun564374Ncoll->GetYaxis()->SetTitle("<N_{part}>");
  hRun564374Ncoll->GetYaxis()->SetTitleSize(0.07);
  hRun564374Ncoll->GetYaxis()->SetTitleOffset(0.6);
  hRun564374Ncoll->Draw("pe");
  hRun564356Ncoll->Draw("same pe");
  hRun564359Ncoll->Draw("same pe");
  hRun564373Ncoll->Draw("same pe");
  hRun564387Ncoll->Draw("same pe");
  hRun564400Ncoll->Draw("same pe");
  hRun564414Ncoll->Draw("same pe");
  hRun564430Ncoll->Draw("same pe");
  hRun564445Ncoll->Draw("same pe");
  legRun->Draw();

  canvNcollVsRun->cd();
  TPad* padNcollBotVsRun = new TPad("padNcollBotVsRun", "", 0, 0, 1, 0.3);
  padNcollBotVsRun->SetTicks(1, 1);
  padNcollBotVsRun->SetTopMargin(0.02);
  padNcollBotVsRun->SetBottomMargin(0.35);
  padNcollBotVsRun->SetLeftMargin(0.1);
  padNcollBotVsRun->SetRightMargin(0.03);
  padNcollBotVsRun->Draw();
  padNcollBotVsRun->cd();

  std::vector<std::pair<TH1F*, TH1F*>> runPairsNcollVsRun = {
    {hRun564356Ncoll, hRun564374Ncoll},
    {hRun564359Ncoll, hRun564374Ncoll},
    {hRun564373Ncoll, hRun564374Ncoll},
    {hRun564387Ncoll, hRun564374Ncoll},
    {hRun564400Ncoll, hRun564374Ncoll},
    {hRun564414Ncoll, hRun564374Ncoll},
    {hRun564430Ncoll, hRun564374Ncoll},
    {hRun564445Ncoll, hRun564374Ncoll}
  };

  std::vector<TH1F*> hRatiosNcollVsRun;
  for (int i = 0; i < runPairsNcollVsRun.size(); ++i) {
    TH1F* hRatio = (TH1F*)runPairsNcollVsRun[i].first->Clone(Form("hRatioNcoll_%d", i));
    hRatio->Divide(hRatio, runPairsNcollVsRun[i].second, 1, 1, "B");
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Run / Ref run");
    hRatio->GetXaxis()->SetTitle("Centrality (%)");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->GetYaxis()->SetRangeUser(0.961, 1.014);
    hRatiosNcollVsRun.push_back(hRatio);
  }
  hRatiosNcollVsRun[0]->Draw("pe");
  for (int i = 1; i < runPairsNcollVsRun.size(); ++i) hRatiosNcollVsRun[i]->Draw("pe same");

  canvNcollVsRun->SaveAs("hNcollOOVsRun.pdf");
}