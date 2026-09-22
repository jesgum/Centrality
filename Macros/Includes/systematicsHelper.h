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
#include <initializer_list>

#include "centralityHelper.h"

namespace syst {

static const std::vector<float> CentBins = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
// static const std::vector<float> CentBins = { 0, 100 };
static const std::vector<double> centBinsD(CentBins.begin(), CentBins.end());

// Hist1D is a template parameter (rather than a fixed TH1F/TH1D) because
// GlauberParameters is filled from ROOT files that store these histograms
// with different precisions depending on their origin (see below).
template <typename Hist1D = TH1F>
inline Hist1D* initGlauParHist(const std::vector<float>& vals, const char* title, const int color)
{
  const bool oldAddDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  Hist1D* hist = new Hist1D(title, "", static_cast<int>(centBinsD.size()) - 1, centBinsD.data());
  TH1::AddDirectory(oldAddDirStatus);
  hist->SetMarkerStyle(kFullCircle);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin) {
    hist->SetBinContent(ibin + 1, vals[ibin]);
    hist->SetBinError(ibin + 1, 1e-4);
  }
  return hist;
}

// Hist2D/Hist1D let callers pick TH2F/TH1F or TH2D/TH1D to match the actual
// type of h2dNpart/h2dNcoll in the source file: the Glauber calibration files
// (Glauber/runCalibration*.cc) write TH2F, while the Trento calibration file
// (Trento/runTrentoCalib.cc) writes TH2D. Defaults match the Glauber files.
template <typename Hist2D = TH2F, typename Hist1D = TH1F>
struct GlauberParameters {
  Hist2D* h2dNpart = nullptr;
  Hist2D* h2dNcoll = nullptr;
  Hist1D* hNpart = nullptr;
  Hist1D* hNcoll = nullptr;
  std::vector<float> npart, ncoll;
  GlauberParameters(const char* path, const bool useFullPath = false);
};

template <typename Hist2D, typename Hist1D>
GlauberParameters<Hist2D, Hist1D>::GlauberParameters(const char* path, const bool useFullPath)
{
  const char* fullpath = useFullPath ? path : Form("../AnalysisResults/%s", path);
  TFile* infile = new TFile(fullpath, "read");
  if (!infile || infile->IsZombie()) {
    std::cerr << "Error: file not found!" << std::endl;
  }

  h2dNpart = dynamic_cast<Hist2D*>(infile->Get("h2dNpart"));
  h2dNcoll = dynamic_cast<Hist2D*>(infile->Get("h2dNcoll"));

  if (!h2dNpart) {
    std::cerr << "Error: h2dNpart not found!" << std::endl;
  }

  if (!h2dNcoll) {
    std::cerr << "Error: h2dNcoll not found!" << std::endl;
  }

  h2dNpart->SetDirectory(0);
  h2dNcoll->SetDirectory(0);

  infile->Close();
  delete infile;

  for (int icent{ 0 }; icent < h2dNcoll->GetNbinsX(); ++icent) {
    TH1D* hProjectionNpart = (TH1D*)h2dNpart->ProjectionY(Form("hProjectionNpart_%d", icent + 1), icent + 1, icent + 1);
    TH1D* hProjectionNcoll = (TH1D*)h2dNcoll->ProjectionY(Form("hProjectionNcoll_%d", icent + 1), icent + 1, icent + 1);
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

  TString histName(path);
  histName.ReplaceAll(".root", "");
  histName.ReplaceAll("/", "_");

  hNpart = initGlauParHist<Hist1D>(npart, Form("hNpart_%s", histName.Data()), kBlack);
  hNcoll = initGlauParHist<Hist1D>(ncoll, Form("hNcoll_%s", histName.Data()), kBlack);
}

template <typename Hist2D, typename Hist1D>
inline void doGlauParQA(Hist2D* h2d, Hist1D* h1d, const char* name)
{
  TCanvas* canv = new TCanvas(Form("canv%s", name), "", 1200, 1000);
  canv->SetLogz();
  h2d->Draw("colz");
  h1d->Draw("same pe");
  canv->SaveAs(Form("qa%s.pdf", name));
  delete canv;
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeSystematics(GPBase base, std::initializer_list<GPVar> variations)
{
  std::vector<std::vector<float>> variationsNpart;
  std::vector<std::vector<float>> variationsNcoll;
  for (const auto& variation : variations) {
    variationsNpart.push_back(cent::calculateDifference(variation.npart, base.npart));
    variationsNcoll.push_back(cent::calculateDifference(variation.ncoll, base.ncoll));
  }

  return { cent::averageOverVariations(variationsNpart), cent::averageOverVariations(variationsNcoll) };
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeRelativeSystematics(GPBase base, std::initializer_list<GPVar> variations)
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
  canv->SetRightMargin(0.03);
  canv->SetBottomMargin(0.13);
}

template <typename Hist1D>
inline void styleHist(Hist1D* hist, const char* titleYaxis)
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
