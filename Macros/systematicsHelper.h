#ifndef SYSTEMATICSHELPER_H_
#define SYSTEMATICSHELPER_H_


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

#include "centralityHelper.h"

namespace syst {

static const std::vector<float> CentBins = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
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

inline TH1F* initGlauParHist(const std::vector<float>& vals, const char* title, const int color)
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

inline void doGlauParQA(TH2F* h2d, TH1F* h1d, const char* name)
{
  TCanvas* canv = new TCanvas(Form("canv%s", name), "", 1200, 1000);
  canv->SetLogz();
  h2d->Draw("colz");
  h1d->Draw("same pe");
  canv->SaveAs(Form("qa%s.pdf", name));
  delete canv;
}

inline std::vector<std::vector<float>> computeSystematics(GlauberParameters base, const std::vector<GlauberParameters>& variations)
{
  std::vector<std::vector<float>> variationsNpart;
  std::vector<std::vector<float>> variationsNcoll;
  for (const auto& variation : variations) {
    variationsNpart.push_back(cent::calculateDifference(variation.npart, base.npart));
    variationsNcoll.push_back(cent::calculateDifference(variation.ncoll, base.ncoll));
  }

  return { cent::averageOverVariations(variationsNpart), cent::averageOverVariations(variationsNcoll) };
}

inline std::vector<std::vector<float>> computeRelativeSystematics(GlauberParameters base, const std::vector<GlauberParameters>& variations)
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

} // namespace syst
#endif // SYSTEMATICSHELPER_H_
