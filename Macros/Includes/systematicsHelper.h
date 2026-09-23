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
inline std::vector<std::vector<float>> computeSystematics(GPBase base, const std::vector<GPVar>& variations)
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
inline std::vector<std::vector<float>> computeSystematics(GPBase base, std::initializer_list<GPVar> variations)
{
  return computeSystematics(base, std::vector<GPVar>(variations));
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeRelativeSystematics(GPBase base, const std::vector<GPVar>& variations)
{
  std::vector<std::vector<float>> variationsNpart;
  std::vector<std::vector<float>> variationsNcoll;
  for (const auto& variation : variations) {
    variationsNpart.push_back(cent::calculateRelativeDifference(variation.npart, base.npart));
    variationsNcoll.push_back(cent::calculateRelativeDifference(variation.ncoll, base.ncoll));
  }

  return { cent::averageOverVariations(variationsNpart), cent::averageOverVariations(variationsNcoll) };
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeRelativeSystematics(GPBase base, std::initializer_list<GPVar> variations)
{
  return computeRelativeSystematics(base, std::vector<GPVar>(variations));
}

inline std::vector<std::vector<float>> combineSystematicsInQuadrature(const std::vector<std::vector<std::vector<float>>>& allSystematics)
{
  if (allSystematics.empty()) {
    return {};
  }

  const size_t nVars = allSystematics.front().size();
  std::vector<std::vector<float>> result(nVars);

  for (size_t iVar = 0; iVar < nVars; ++iVar) {
    const size_t nBins = allSystematics.front()[iVar].size();
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

template <typename... Rest>
inline std::vector<std::vector<float>> combineSystematicsInQuadrature(const std::vector<std::vector<float>>& first, const Rest&... rest)
{
  return combineSystematicsInQuadrature(std::vector<std::vector<std::vector<float>>>{ first, rest... });
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

// A single named source of systematic uncertainty (e.g. "Nuclei shape",
// "Trento", "Anchor shift"). Keeping sources in a vector, each with its own
// `draw`/`includeInTotal` switches, is what makes it possible to add a new
// source or toggle an existing one on/off without touching the plotting code.
struct SystSource {
  std::string label;
  int color = kBlack;
  bool draw = true;          // show this source's curve (and legend entry) on the canvases
  bool includeInTotal = true; // fold this source into the quadrature-summed "Total"
  std::vector<std::vector<float>> syst;    // indexed by cent::GlauberQuantity (kNpart, kNcoll)
  std::vector<std::vector<float>> relSyst; // indexed by cent::GlauberQuantity (kNpart, kNcoll)
  bool asMarkers = false; // draw as markers instead of a line (used for external references, e.g. Run 2)
  int marker = kFullCircle;
  const std::vector<float>& values(int var, bool relative) const
  {
    return relative ? relSyst[var] : syst[var];
  }

  void print(int var, bool relative = true)
  {
    const std::vector<float>& sys = relative ? relSyst[var] : syst[var];
    for (size_t ii = 0; ii < sys.size(); ++ii) {
      std::cout << CentBins[ii] << "-" << CentBins[ii + 1] << "% | " << sys[ii] << std::endl;
    }
  }
};

template <typename GPBase, typename GPVar>
inline SystSource makeSystSource(const char* label, int color, bool draw, bool includeInTotal, GPBase base, const std::vector<GPVar>& variations)
{
  return { label, color, draw, includeInTotal, computeSystematics(base, variations), computeRelativeSystematics(base, variations) };
}

template <typename GPBase, typename GPVar>
inline SystSource makeSystSource(const char* label, int color, bool draw, bool includeInTotal, GPBase base, std::initializer_list<GPVar> variations)
{
  return makeSystSource(label, color, draw, includeInTotal, base, std::vector<GPVar>(variations));
}

// Builds a reference source from externally quoted values and uncertainties
// (e.g. the published Run 2 <Npart>/<Ncoll>), one entry per CentBins interval.
// Bins without a reference value should be given a value <= 0; they are then
// set to 0 and fall outside the log-scale frame, so they are not drawn.
// References are never folded into the total and are drawn as markers.
inline SystSource makeReferenceSource(const char* label, int color, const std::vector<float>& npartVals, const std::vector<float>& npartErrs, const std::vector<float>& ncollVals, const std::vector<float>& ncollErrs)
{
  const size_t nBins = CentBins.size() - 1;
  if (npartVals.size() != nBins || npartErrs.size() != nBins || ncollVals.size() != nBins || ncollErrs.size() != nBins) {
    std::cerr << "Error in makeReferenceSource(" << label << "): expected " << nBins << " values per vector (one per CentBins interval)!" << std::endl;
    return {};
  }

  auto relative = [](const std::vector<float>& vals, const std::vector<float>& errs) {
    std::vector<float> rel;
    for (size_t ii = 0; ii < vals.size(); ++ii) {
      rel.push_back(vals[ii] > 0 ? errs[ii] / vals[ii] : 0.f);
    }
    return rel;
  };
  auto absolute = [](const std::vector<float>& vals, const std::vector<float>& errs) {
    std::vector<float> abs;
    for (size_t ii = 0; ii < vals.size(); ++ii) {
      abs.push_back(vals[ii] > 0 ? errs[ii] : 0.f);
    }
    return abs;
  };

  SystSource ref;
  ref.label = label;
  ref.color = color;
  ref.draw = true;
  ref.includeInTotal = false;
  ref.asMarkers = true;
  ref.syst = { absolute(npartVals, npartErrs), absolute(ncollVals, ncollErrs) };
  ref.relSyst = { relative(npartVals, npartErrs), relative(ncollVals, ncollErrs) };
  return ref;
}

// Quadrature-sums the sources flagged `includeInTotal` into a new SystSource.
inline SystSource combineSourcesInTotal(const char* label, int color, bool draw, const std::vector<SystSource>& sources)
{
  std::vector<std::vector<std::vector<float>>> allSyst, allRelSyst;
  for (const auto& source : sources) {
    if (!source.includeInTotal) {
      continue;
    }
    allSyst.push_back(source.syst);
    allRelSyst.push_back(source.relSyst);
  }

  SystSource total;
  total.label = label;
  total.color = color;
  total.draw = draw;
  total.includeInTotal = false;
  total.syst = combineSystematicsInQuadrature(allSyst);
  total.relSyst = combineSystematicsInQuadrature(allRelSyst);
  return total;
}

// Draws one canvas for a given quantity (cent::kNpart or cent::kNcoll) and
// value kind (absolute or relative), overlaying every source with `draw ==
// true`, plus the total if `total.draw == true`. The (always invisible)
// total histogram is used only to fix the axis range/style, so the frame
// stays stable regardless of which sources are toggled on.
inline void drawSystematics(const std::vector<SystSource>& sources, const SystSource& total, int var, bool relative, const char* titleYaxis, const char* outFile, const float min = 2e-6, const float max = 2e-1)
{
  TH1F* hFrame = initGlauParHist(total.values(var, relative), Form("hFrame_%s", outFile), total.color);
  hFrame->SetMarkerColorAlpha(0, 0);
  hFrame->SetLineColorAlpha(0, 0);

  TCanvas* canv = new TCanvas(Form("canv_%s", outFile), "", 1600, 1200);
  styleCanvas(canv);
  styleHist(hFrame, titleYaxis);
  hFrame->SetMinimum(min);
  hFrame->SetMaximum(max);
  hFrame->Draw();
  canv->SetLogy();

  TLegend* leg = new TLegend(0.18, 0.66, 0.63, 0.94);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);

  TH1F* hTotal = nullptr;
  if (total.draw) {
    hTotal = initGlauParHist(total.values(var, relative), Form("hTotal_%s", outFile), total.color);
    leg->AddEntry(hTotal, total.label.c_str(), "l");
  }

  int isrc = 0;
  for (const auto& source : sources) {
    if (!source.draw) {
      continue;
    }
    TH1F* h = initGlauParHist(source.values(var, relative), Form("hSource%d_%s", isrc++, outFile), source.color);
    if (source.asMarkers) {
      h->SetMarkerSize(1.5);
      h->SetMarkerStyle(source.marker);
      h->Draw("p same");
      leg->AddEntry(h, source.label.c_str(), "pl");
    } else {
      h->Draw("hist same");
      leg->AddEntry(h, source.label.c_str(), "l");
    }
  }

  if (hTotal) {
    hTotal->Draw("hist same");
  }

  leg->Draw();
  canv->SaveAs(outFile);
}

} // namespace syst
#endif // SYSTEMATICSHELPER_H_
