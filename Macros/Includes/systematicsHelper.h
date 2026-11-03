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
#include <algorithm>

#include "centralityHelper.h"

namespace syst {

// Default centrality binning (integer edges, in %) of everything below. It is
// read whenever a GlauberParameters/histogram is created, so a macro can switch
// it (e.g. to { 0, 100 }) by reassigning it before building anything.
inline std::vector<float> CentBins = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };

// Binning used by SystSource::print(), independent of the (plotting) CentBins.
// Integer edges in %. A macro switches it in one place, e.g.
//   syst::PrintBins = { 0, 10, 30, 50, 70, 100 };
inline std::vector<float> PrintBins = CentBins;

// Printable label of centrality interval ii of bins, e.g. "0-5%"
inline TString binLabel(const std::vector<float>& bins, const size_t ii)
{
  return Form("%g-%g%%", bins[ii], bins[ii + 1]);
}

// Hist1D is a template parameter (rather than a fixed TH1F/TH1D) because
// GlauberParameters is filled from ROOT files that store these histograms
// with different precisions depending on their origin (see below).
template <typename Hist1D = TH1F>
inline Hist1D* initGlauParHist(const std::vector<float>& vals, const char* title, const int color, const std::vector<float>& bins = CentBins)
{
  const std::vector<double> binsD(bins.begin(), bins.end());
  const bool oldAddDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  Hist1D* hist = new Hist1D(title, "", static_cast<int>(binsD.size()) - 1, binsD.data());
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

// Averages per-1% values (index = centrality %) into the given integer bin edges
inline std::vector<float> rebinValues(const std::vector<float>& vals, const std::vector<float>& bins)
{
  std::vector<float> result;
  for (int ii = 0; ii < static_cast<int>(bins.size()) - 1; ++ii) {
    const int thisBin = bins[ii];
    const int nextBin = bins[ii + 1];

    std::vector<float> theseValues, theseWeights;
    for (int jj = thisBin; jj < nextBin; ++jj) {
      theseValues.push_back(vals[jj]);
      theseWeights.push_back(1);
    }

    result.push_back(cent::weightedAverage(theseValues, theseWeights));
  }
  return result;
}

// RMS of the equal-weight mixture of the per-1% distributions in each bin
// (same weighting as rebinValues): var = <rms_i^2 + mean_i^2> - <mean_i>^2
inline std::vector<float> rebinRMS(const std::vector<float>& means, const std::vector<float>& rms, const std::vector<float>& bins)
{
  std::vector<float> secondMoments;
  for (size_t ii = 0; ii < means.size(); ++ii) {
    secondMoments.push_back(rms[ii] * rms[ii] + means[ii] * means[ii]);
  }
  const std::vector<float> rebinnedMeans = rebinValues(means, bins);
  const std::vector<float> rebinnedSecondMoments = rebinValues(secondMoments, bins);
  std::vector<float> result;
  for (size_t ii = 0; ii < rebinnedMeans.size(); ++ii) {
    result.push_back(std::sqrt(std::max(0.f, rebinnedSecondMoments[ii] - rebinnedMeans[ii] * rebinnedMeans[ii])));
  }
  return result;
}

// Plain <Npart>/<Ncoll> values (and the RMS of their distributions), shaped
// like GlauberParameters so that computeSystematics() works on both
struct GlauberValues {
  std::vector<float> npart, ncoll;
  std::vector<float> npartRMS, ncollRMS;
};

// Per-1% GlauberValues averaged into the given bins
inline GlauberValues rebinGlauberValues(const GlauberValues& fine, const std::vector<float>& bins)
{
  return { rebinValues(fine.npart, bins), rebinValues(fine.ncoll, bins),
           rebinRMS(fine.npart, fine.npartRMS, bins), rebinRMS(fine.ncoll, fine.ncollRMS, bins) };
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
  Hist2D* h2dNancestor = nullptr; // optional: not produced by every calibration
  Hist1D* hNancestor = nullptr;   // only filled if h2dNancestor is available
  Hist2D* h2dMult = nullptr;      // optional: not produced by every calibration
  Hist1D* hMult = nullptr;        // only filled if h2dMult is available
  std::vector<float> npart, ncoll, nancestor, mult;
  std::vector<float> npartRMS, ncollRMS; // RMS of the Npart/Ncoll distribution in each bin
  GlauberValues fine; // per-1% npart/ncoll before rebinning, to re-evaluate at other binnings
  std::vector<float> bins; // centrality bin edges used for the rebinned vectors and 1D histograms
  GlauberParameters(const char* path, const bool useFullPath = false, const std::vector<float>& centBins = CentBins);
};

template <typename Hist2D, typename Hist1D>
GlauberParameters<Hist2D, Hist1D>::GlauberParameters(const char* path, const bool useFullPath, const std::vector<float>& centBins)
  : bins(centBins)
{
  const TString fullpath = useFullPath ? TString(path) : TString::Format("../AnalysisResults/%s", path);
  TFile* infile = new TFile(fullpath, "read");
  if (!infile || infile->IsZombie()) {
    std::cerr << "Error: file not found!" << std::endl;
  }

  h2dNpart = dynamic_cast<Hist2D*>(infile->Get("h2dNpart"));
  h2dNcoll = dynamic_cast<Hist2D*>(infile->Get("h2dNcoll"));

  // dynamic_cast also fails if the histogram exists with the other precision
  // (TH2F vs TH2D), so name the stored type to make that case obvious
  auto reportMissing = [&](const char* name) {
    TObject* obj = infile->Get(name);
    if (obj) {
      std::cerr << "Error: " << name << " in " << fullpath << " is a " << obj->ClassName() << ", use GlauberParameters<" << obj->ClassName() << "> for it!" << std::endl;
    } else {
      std::cerr << "Error: " << name << " not found in " << fullpath << "!" << std::endl;
    }
  };
  if (!h2dNpart) {
    reportMissing("h2dNpart");
  }
  if (!h2dNcoll) {
    reportMissing("h2dNcoll");
  }

  // optional, silently skipped if absent
  h2dNancestor = dynamic_cast<Hist2D*>(infile->Get("h2dNancestors"));
  h2dMult = dynamic_cast<Hist2D*>(infile->Get("h2dMult"));

  h2dNpart->SetDirectory(0);
  h2dNcoll->SetDirectory(0);
  if (h2dNancestor) {
    h2dNancestor->SetDirectory(0);
  }
  if (h2dMult) {
    h2dMult->SetDirectory(0);
  }

  infile->Close();
  delete infile;

  for (int icent{ 0 }; icent < h2dNcoll->GetNbinsX(); ++icent) {
    TH1D* hProjectionNpart = (TH1D*)h2dNpart->ProjectionY(Form("hProjectionNpart_%d", icent + 1), icent + 1, icent + 1);
    TH1D* hProjectionNcoll = (TH1D*)h2dNcoll->ProjectionY(Form("hProjectionNcoll_%d", icent + 1), icent + 1, icent + 1);
    npart.push_back(hProjectionNpart->GetMean());
    ncoll.push_back(hProjectionNcoll->GetMean());
    npartRMS.push_back(hProjectionNpart->GetRMS());
    ncollRMS.push_back(hProjectionNcoll->GetRMS());
    delete hProjectionNpart;
    delete hProjectionNcoll;
    if (h2dNancestor) {
      TH1D* hProjectionNancestor = (TH1D*)h2dNancestor->ProjectionY(Form("hProjectionNancestor_%d", icent + 1), icent + 1, icent + 1);
      nancestor.push_back(hProjectionNancestor->GetMean());
      delete hProjectionNancestor;
    }
    if (h2dMult) {
      TH1D* hProjectionMult = (TH1D*)h2dMult->ProjectionY(Form("hProjectionMult_%d", icent + 1), icent + 1, icent + 1);
      mult.push_back(hProjectionMult->GetMean());
      delete hProjectionMult;
    }
  }

  fine = { npart, ncoll, npartRMS, ncollRMS };

  const GlauberValues rebinned = rebinGlauberValues(fine, bins);
  npart = rebinned.npart;
  ncoll = rebinned.ncoll;
  npartRMS = rebinned.npartRMS;
  ncollRMS = rebinned.ncollRMS;
  if (h2dNancestor) {
    nancestor = rebinValues(nancestor, bins);
  }
  if (h2dMult) {
    mult = rebinValues(mult, bins);
  }

  TString histName(path);
  histName.ReplaceAll(".root", "");
  histName.ReplaceAll("/", "_");

  hNpart = initGlauParHist<Hist1D>(npart, Form("hNpart_%s", histName.Data()), kBlack, bins);
  hNcoll = initGlauParHist<Hist1D>(ncoll, Form("hNcoll_%s", histName.Data()), kBlack, bins);
  if (h2dNancestor) {
    hNancestor = initGlauParHist<Hist1D>(nancestor, Form("hNancestor_%s", histName.Data()), kBlack, bins);
  }
  if (h2dMult) {
    hMult = initGlauParHist<Hist1D>(mult, Form("hMult_%s", histName.Data()), kBlack, bins);
  }
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

// How the |variation - base| differences of one source are combined, per bin
enum class Combination {
  kAverage, // mean of the differences (e.g. a +-shift pair)
  kRMS      // sqrt(sum d^2 / N), for a set of alternative models (e.g. nuclear shapes)
};

inline std::vector<float> combineVariations(const std::vector<std::vector<float>>& differences, Combination combination)
{
  return combination == Combination::kRMS ? cent::rmsOverVariations(differences) : cent::averageOverVariations(differences);
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeSystematics(GPBase base, const std::vector<GPVar>& variations, Combination combination = Combination::kAverage)
{
  std::vector<std::vector<float>> variationsNpart;
  std::vector<std::vector<float>> variationsNcoll;
  for (const auto& variation : variations) {
    variationsNpart.push_back(cent::calculateDifference(variation.npart, base.npart));
    variationsNcoll.push_back(cent::calculateDifference(variation.ncoll, base.ncoll));
  }

  return { combineVariations(variationsNpart, combination), combineVariations(variationsNcoll, combination) };
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeSystematics(GPBase base, std::initializer_list<GPVar> variations, Combination combination = Combination::kAverage)
{
  return computeSystematics(base, std::vector<GPVar>(variations), combination);
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeRelativeSystematics(GPBase base, const std::vector<GPVar>& variations, Combination combination = Combination::kAverage)
{
  std::vector<std::vector<float>> variationsNpart;
  std::vector<std::vector<float>> variationsNcoll;
  for (const auto& variation : variations) {
    variationsNpart.push_back(cent::calculateRelativeDifference(variation.npart, base.npart));
    variationsNcoll.push_back(cent::calculateRelativeDifference(variation.ncoll, base.ncoll));
  }

  return { combineVariations(variationsNpart, combination), combineVariations(variationsNcoll, combination) };
}

template <typename GPBase, typename GPVar>
inline std::vector<std::vector<float>> computeRelativeSystematics(GPBase base, std::initializer_list<GPVar> variations, Combination combination = Combination::kAverage)
{
  return computeRelativeSystematics(base, std::vector<GPVar>(variations), combination);
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
  std::vector<std::vector<float>> central; // base <Npart>/<Ncoll> the systematics refer to, same indexing
  std::vector<std::vector<float>> rms;     // RMS of the base Npart/Ncoll distributions, same indexing (empty if unknown)
  std::vector<float> bins = CentBins;      // centrality bin edges of all the vectors above
  // What is needed to re-evaluate the source at another binning (see rebinned()):
  GlauberValues baseFine;                    // per-1% base values
  std::vector<GlauberValues> variationsFine; // per-1% values of each variation
  Combination combination = Combination::kAverage; // how the variations are combined
  std::vector<SystSource> components;        // totals only: the sources folded in
  bool asMarkers = false; // draw as markers instead of a line (used for external references, e.g. Run 2)
  int marker = kFullCircle;
  int lineWidth = 1;
  int lineStyle = kSolid;

  void setLineStyle(int col, int width, int style)
  {
    color = col;
    lineWidth = width;
    lineStyle = style;
  }

  const std::vector<float>& values(int var, bool relative) const
  {
    return relative ? relSyst[var] : syst[var];
  }

  static const char* quantityName(int var)
  {
    return var == cent::kNpart ? "Npart" : var == cent::kNcoll ? "Ncoll" : "?";
  }

  // "<value> (RMS <rms>) +- <abs. syst> (<rel. syst>%)" for one quantity and bin
  TString formatBin(int var, size_t ii) const
  {
    TString out = central.empty() ? "" : Form("%8.2f ", central[var][ii]);
    if (!rms.empty() && !rms[var].empty()) {
      out += Form("(RMS %7.2f) ", rms[var][ii]);
    }
    return out + Form("+- %7.3f (%5.2f%%)", syst[var][ii], 100. * relSyst[var][ii]);
  }

  // This source re-evaluated at newBins: the base and variations (or, for a
  // total, its components) are rebinned from 1% first, then compared, so the
  // result is the same as building everything with newBins from the start.
  SystSource rebinned(const std::vector<float>& newBins) const
  {
    if (newBins == bins) {
      return *this;
    }
    SystSource out = *this;
    out.bins = newBins;
    if (!components.empty()) {
      std::vector<std::vector<std::vector<float>>> allSyst, allRelSyst;
      for (auto& component : out.components) {
        component = component.rebinned(newBins);
        allSyst.push_back(component.syst);
        allRelSyst.push_back(component.relSyst);
      }
      out.syst = combineSystematicsInQuadrature(allSyst);
      out.relSyst = combineSystematicsInQuadrature(allRelSyst);
      out.central = out.components.front().central;
      out.rms = out.components.front().rms;
    } else if (!variationsFine.empty()) {
      const GlauberValues base = rebinGlauberValues(baseFine, newBins);
      std::vector<GlauberValues> variations;
      for (const auto& variation : variationsFine) {
        variations.push_back(rebinGlauberValues(variation, newBins));
      }
      out.syst = computeSystematics(base, variations, combination);
      out.relSyst = computeRelativeSystematics(base, variations, combination);
      out.central = { base.npart, base.ncoll };
      out.rms = { base.npartRMS, base.ncollRMS };
    } else {
      std::cerr << "Warning: " << label << " cannot be rebinned (no per-1% values, e.g. a reference), keeping its own binning" << std::endl;
      return *this;
    }
    return out;
  }

  // One quantity (cent::kNpart or cent::kNcoll), one line per bin of printBins
  void print(int var, const std::vector<float>& printBins = PrintBins) const
  {
    const SystSource source = rebinned(printBins);
    std::cout << label << " | <" << quantityName(var) << ">" << std::endl;
    for (size_t ii = 0; ii < source.syst[var].size(); ++ii) {
      std::cout << Form("%-8s | %s", binLabel(source.bins, ii).Data(), source.formatBin(var, ii).Data()) << std::endl;
    }
  }

  // Npart and Ncoll side by side, one line per bin of printBins
  void print(const std::vector<float>& printBins = PrintBins) const
  {
    const SystSource source = rebinned(printBins);
    std::cout << label << std::endl;
    for (size_t ii = 0; ii < source.syst[cent::kNpart].size(); ++ii) {
      std::cout << Form("%-8s | Npart: %s | Ncoll: %s", binLabel(source.bins, ii).Data(), source.formatBin(cent::kNpart, ii).Data(), source.formatBin(cent::kNcoll, ii).Data()) << std::endl;
    }
  }
};

template <typename GPBase, typename GPVar>
inline SystSource makeSystSource(const char* label, int color, bool draw, bool includeInTotal, GPBase base, const std::vector<GPVar>& variations, Combination combination = Combination::kAverage)
{
  SystSource source;
  source.label = label;
  source.color = color;
  source.draw = draw;
  source.includeInTotal = includeInTotal;
  source.combination = combination;
  source.syst = computeSystematics(base, variations, combination);
  source.relSyst = computeRelativeSystematics(base, variations, combination);
  source.central = { base.npart, base.ncoll };
  source.rms = { base.npartRMS, base.ncollRMS };
  source.bins = base.bins;
  source.baseFine = base.fine;
  for (const auto& variation : variations) {
    source.variationsFine.push_back(variation.fine);
  }
  return source;
}

template <typename GPBase, typename GPVar>
inline SystSource makeSystSource(const char* label, int color, bool draw, bool includeInTotal, GPBase base, std::initializer_list<GPVar> variations, Combination combination = Combination::kAverage)
{
  return makeSystSource(label, color, draw, includeInTotal, base, std::vector<GPVar>(variations), combination);
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
  ref.central = { npartVals, ncollVals };
  return ref;
}

// Quadrature-sums the sources flagged `includeInTotal` into a new SystSource.
inline SystSource combineSourcesInTotal(const char* label, int color, bool draw, const std::vector<SystSource>& sources, int lineWidth = 1, int lineStyle = kSolid)
{
  SystSource total;
  std::vector<std::vector<std::vector<float>>> allSyst, allRelSyst;
  const SystSource* first = nullptr; // central values and binning are taken from the first included source
  for (const auto& source : sources) {
    if (!source.includeInTotal) {
      continue;
    }
    if (!first) {
      first = &source;
    } else if (source.bins != first->bins) {
      std::cerr << "Error in combineSourcesInTotal(" << label << "): " << source.label << " has different binning than " << first->label << "!" << std::endl;
    }
    allSyst.push_back(source.syst);
    allRelSyst.push_back(source.relSyst);
    total.components.push_back(source);
  }

  total.label = label;
  total.color = color;
  total.lineWidth = lineWidth;
  total.lineStyle = lineStyle;
  total.draw = draw;
  total.includeInTotal = false;
  total.syst = combineSystematicsInQuadrature(allSyst);
  total.relSyst = combineSystematicsInQuadrature(allRelSyst);
  if (first) {
    total.central = first->central;
    total.rms = first->rms;
    total.bins = first->bins;
  }
  return total;
}

// Draws one canvas for a given quantity (cent::kNpart or cent::kNcoll) and
// value kind (absolute or relative), overlaying every source with `draw ==
// true`, plus every total with `draw == true`. Passing several totals makes it
// possible to compare different ways of combining the sources on one canvas.
// The (always invisible) first total is used only to fix the axis range/style,
// so the frame stays stable regardless of which sources are toggled on.
inline void drawSystematics(const std::vector<SystSource>& sources, const std::vector<SystSource>& totals, int var, bool relative, const char* titleYaxis, const char* outFile, const float min = 2e-6, const float max = 2e-1)
{
  if (totals.empty()) {
    std::cerr << "Error in drawSystematics: no total given!" << std::endl;
    return;
  }

  TH1F* hFrame = initGlauParHist(totals.front().values(var, relative), Form("hFrame_%s", outFile), totals.front().color, totals.front().bins);
  hFrame->SetMarkerColorAlpha(0, 0);
  hFrame->SetLineColorAlpha(0, 0);

  TCanvas* canv = new TCanvas(Form("canv_%s", outFile), "", 1600, 1200);
  styleCanvas(canv);
  styleHist(hFrame, titleYaxis);
  hFrame->SetMinimum(min);
  hFrame->SetMaximum(max);
  hFrame->Draw();
  // canv->SetLogy();

  TLegend* leg = new TLegend(0.18, 0.66, 0.63, 0.94);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);

  auto applyLineStyle = [](TH1F* h, const SystSource& source) {
    h->SetLineColor(source.color);
    h->SetLineWidth(source.lineWidth);
    h->SetLineStyle(source.lineStyle);
  };

  std::vector<TH1F*> hTotals;
  int itot = 0;
  for (const auto& total : totals) {
    if (!total.draw) {
      continue;
    }
    TH1F* h = initGlauParHist(total.values(var, relative), Form("hTotal%d_%s", itot++, outFile), total.color, total.bins);
    applyLineStyle(h, total);
    leg->AddEntry(h, total.label.c_str(), "l");
    hTotals.push_back(h);
  }

  int isrc = 0;
  for (const auto& source : sources) {
    if (!source.draw) {
      continue;
    }
    TH1F* h = initGlauParHist(source.values(var, relative), Form("hSource%d_%s", isrc++, outFile), source.color, source.bins);
    applyLineStyle(h, source);
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

  for (const auto& h : hTotals) {
    h->Draw("hist same");
  }

  leg->Draw();
  canv->SaveAs(outFile);
}

inline void drawSystematics(const std::vector<SystSource>& sources, const SystSource& total, int var, bool relative, const char* titleYaxis, const char* outFile, const float min = 2e-6, const float max = 2e-1)
{
  drawSystematics(sources, std::vector<SystSource>{ total }, var, relative, titleYaxis, outFile, min, max);
}

} // namespace syst
#endif // SYSTEMATICSHELPER_H_
