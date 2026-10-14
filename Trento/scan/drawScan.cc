#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TPRegexp.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "colorManager.h"

// Draws the output of runScan.sh together:
//   scan<tag>_Fits.pdf    data vs. Trento + NBD fit of every scan point (fits normalized to the data in the fit range), bottom: data/fit
//   scan<tag>_Entropy.pdf entropy distributions from saveTrento.cc (normalized)
//   scan<tag>_Npart.pdf   <Npart> vs. centrality, bottom: ratio to the reference point
//   scan<tag>_Ncoll.pdf   <Ncoll> vs. centrality, bottom: ratio to the reference point
//   scan<tag>_Parameters.pdf deviation from the reference and scan parameter values vs. rank (most accurate left)

namespace
{
const std::vector<double> CentBins = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
const int RebinData = 20;
const int RebinEntropy = 10;
const double CentMaxMatch = 90; // centrality bins above this are ignored when looking for the point closest to the reference
const size_t NClosestToPrint = 5;

struct ScanPoint {
  std::string name;
  TString label;           // legend entry
  TH1* hData = nullptr;    // data, rebinned
  TH1* hFit = nullptr;     // Trento + NBD fit, rebinned and normalized to the data in the fit range
  TH1* hEntropy = nullptr; // normalized entropy distribution
  TH1* hNpart = nullptr;   // <Npart> in CentBins
  TH1* hNcoll = nullptr;   // <Ncoll> in CentBins
  double mu = 0, k = 0, chi2 = 0, ndf = 0, fitMin = 0, fitMax = 0, anchorPoint = 0;
  std::string trentoConfig; // contents of trento.conf, "key = value" per line
};

// RMS of the relative deviation from the reference over the centrality bins up to CentMaxMatch
double relativeRms(const TH1* hist, const TH1* ref)
{
  double sum = 0;
  int n = 0;
  for (int ibin = 1; ibin <= ref->GetNbinsX(); ++ibin) {
    if (ref->GetXaxis()->GetBinUpEdge(ibin) > CentMaxMatch + 1e-6 || ref->GetBinContent(ibin) <= 0) {
      continue;
    }
    const double dev = (hist->GetBinContent(ibin) - ref->GetBinContent(ibin)) / ref->GetBinContent(ibin);
    sum += dev * dev;
    ++n;
  }
  return n > 0 ? std::sqrt(sum / n) : 0;
}

// "norm_15_p_0.5" -> {{"norm", 15}, {"p", 0.5}}; empty if the name is not a param_value list
std::vector<std::pair<std::string, double>> parseParameters(const std::string& name)
{
  std::vector<std::pair<std::string, double>> pars;
  TObjArray* tokens = TString(name).Tokenize("_");
  if (tokens->GetEntries() % 2 == 0) {
    for (int i = 0; i + 1 < tokens->GetEntries(); i += 2) {
      const TString value = tokens->At(i + 1)->GetName();
      if (value.IsFloat()) {
        pars.emplace_back(tokens->At(i)->GetName(), value.Atof());
      }
    }
  }
  delete tokens;
  return pars;
}

// "norm_15_p_0.5" -> "norm = 15, p = 0.5" (an underscore is a subscript in TLatex)
TString makeLabel(const std::string& name)
{
  TObjArray* tokens = TString(name).Tokenize("_");
  TString label;
  for (int i = 0; i + 1 < tokens->GetEntries(); i += 2) {
    label += Form("%s%s = %s", label.IsNull() ? "" : ", ", tokens->At(i)->GetName(), tokens->At(i + 1)->GetName());
  }
  if (tokens->GetEntries() % 2 == 1) {
    label = name.c_str(); // not a param_value list, e.g. "default"
  }
  delete tokens;
  return label;
}

TH1* getHist(TFile* file, const char* histName)
{
  TH1* hist = dynamic_cast<TH1*>(file->Get(histName));
  if (!hist) {
    std::cerr << "Error: " << histName << " not found in " << file->GetName() << std::endl;
  }
  return hist;
}

// Mean of y per 1% centrality bin, averaged (equal weights) into CentBins
TH1D* meanVsCentrality(TH2* h2d, const char* name)
{
  TH1D* hist = new TH1D(name, "", static_cast<int>(CentBins.size()) - 1, CentBins.data());
  for (int ibin = 0; ibin < hist->GetNbinsX(); ++ibin) {
    const int first = h2d->GetXaxis()->FindBin(CentBins[ibin] + 1e-3);
    const int last = h2d->GetXaxis()->FindBin(CentBins[ibin + 1] - 1e-3);
    double sum = 0;
    for (int icent = first; icent <= last; ++icent) {
      TH1D* hProjection = h2d->ProjectionY("hProjection", icent, icent);
      sum += hProjection->GetMean();
      delete hProjection;
    }
    hist->SetBinContent(ibin + 1, sum / (last - first + 1));
    hist->SetBinError(ibin + 1, 1e-4);
  }
  return hist;
}

// Reference from a calibration file (Trento or Glauber), only h2dNpart and h2dNcoll are needed
bool loadReference(const TString& path, ScanPoint& point)
{
  TFile* file = TFile::Open(path, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: cannot open reference file " << path << std::endl;
    return false;
  }
  TH2* h2dNpart = dynamic_cast<TH2*>(getHist(file, "h2dNpart"));
  TH2* h2dNcoll = dynamic_cast<TH2*>(getHist(file, "h2dNcoll"));
  if (!h2dNpart || !h2dNcoll) {
    return false;
  }
  point.hNpart = meanVsCentrality(h2dNpart, "hNpart_reference");
  point.hNcoll = meanVsCentrality(h2dNcoll, "hNcoll_reference");
  file->Close();
  return true;
}

bool loadPoint(const TString& dir, ScanPoint& point)
{
  TFile* fileOut = TFile::Open(dir + "/outfile.root", "READ");
  TFile* fileFit = TFile::Open(dir + "/trentofit.root", "READ");
  TFile* fileCalib = TFile::Open(dir + "/trentofit_calib.root", "READ");
  if (!fileOut || !fileFit || !fileCalib) {
    return false;
  }

  TH1* hDataFine = getHist(fileFit, "hV0MUltraFine");
  TH1* hFitFine = getHist(fileFit, "hGlauber");
  TH1* hPars = getHist(fileFit, "hGlauberParameters");
  TH1* hFitRange = getHist(fileFit, "hGlauberFitRange");
  TH1* hAnchor = getHist(fileCalib, "hAnchorPoint");
  TH1* hEntropy = getHist(fileOut, "hEntropy");
  TH2* h2dNpart = dynamic_cast<TH2*>(getHist(fileCalib, "h2dNpart"));
  TH2* h2dNcoll = dynamic_cast<TH2*>(getHist(fileCalib, "h2dNcoll"));
  if (!hDataFine || !hFitFine || !hPars || !hFitRange || !hAnchor || !hEntropy || !h2dNpart || !h2dNcoll) {
    return false;
  }

  const char* name = point.name.c_str();
  point.hData = hDataFine->Rebin(RebinData, Form("hData_%s", name));
  point.hFit = hFitFine->Rebin(RebinData, Form("hFit_%s", name));
  point.hEntropy = hEntropy->Rebin(RebinEntropy, Form("hEntropy_%s", name));
  point.hEntropy->Scale(1. / point.hEntropy->Integral(), "width");
  point.hNpart = meanVsCentrality(h2dNpart, Form("hNpart_%s", name));
  point.hNcoll = meanVsCentrality(h2dNcoll, Form("hNcoll_%s", name));
  point.mu = hPars->GetBinContent(1);
  point.k = hPars->GetBinContent(2);
  point.chi2 = hPars->GetBinContent(6);
  point.ndf = hPars->GetBinContent(7);
  point.fitMin = hFitRange->GetBinContent(1);
  point.fitMax = hFitRange->GetBinContent(2);

  // hGlauber is evaluated on the fine binning but normalized to the rebinned data used in the fit:
  // normalize it to the data inside the fit range instead
  const int binMin = point.hData->FindBin(point.fitMin);
  const int binMax = point.hData->FindBin(point.fitMax);
  point.hFit->Scale(point.hData->Integral(binMin, binMax) / point.hFit->Integral(binMin, binMax));
  point.anchorPoint = hAnchor->GetBinContent(1);

  std::ifstream config((dir + "/trento.conf").Data());
  std::string line;
  while (std::getline(config, line)) {
    point.trentoConfig += "    " + line + "\n";
  }

  fileOut->Close();
  fileFit->Close();
  fileCalib->Close();
  return true;
}

void setStyle(TH1* hist, int color, int marker = kFullCircle)
{
  hist->SetTitle("");
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(0.8);
  hist->SetLineWidth(2);
}

TLegend* makeLegend(size_t nEntries, double x1, double y1, double x2, double y2)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->SetTextSize(0.035);
  leg->SetNColumns(std::max<int>(1, std::ceil(nEntries / 12.)));
  return leg;
}

// Canvas split in a top pad (values) and a bottom pad (ratios), as in Macros/drawNpartNcollTrento.cc
void makePads(const char* name, TPad*& padTop, TPad*& padBot, bool logyTop)
{
  TCanvas* canv = new TCanvas(Form("canv_%s", name), "", 1200, 800);
  padTop = new TPad(Form("padTop_%s", name), "", 0, 0.3, 1, 1);
  padTop->SetTicks(1, 1);
  padTop->SetTopMargin(0.04);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.12);
  padTop->SetRightMargin(0.03);
  padTop->SetLogy(logyTop);
  padTop->Draw();
  canv->cd();
  padBot = new TPad(Form("padBot_%s", name), "", 0, 0, 1, 0.3);
  padBot->SetTicks(1, 1);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.12);
  padBot->SetRightMargin(0.03);
  padBot->Draw();
}

void styleTopFrame(TH1* frame, const char* yTitle)
{
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle(yTitle);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.1);
}

void styleRatioFrame(TH1* frame, const char* xTitle, const char* yTitle, double yMin, double yMax)
{
  frame->SetTitle("");
  frame->SetMinimum(yMin);
  frame->SetMaximum(yMax);
  frame->GetYaxis()->SetTitle(yTitle);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleSize(0.12);
  frame->GetYaxis()->SetTitleOffset(0.45);
  frame->GetYaxis()->SetLabelSize(0.10);
  frame->GetXaxis()->SetTitle(xTitle);
  frame->GetXaxis()->SetTitleSize(0.13);
  frame->GetXaxis()->SetLabelSize(0.11);
}

void drawUnityLine(TH1* frame, double xMax)
{
  TLine* line = new TLine(frame->GetXaxis()->GetXmin(), 1, xMax, 1);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 1);
  line->Draw("same");
}

void drawFits(const std::vector<ScanPoint>& points, const TString& outFile)
{
  TPad *padTop, *padBot;
  makePads("Fits", padTop, padBot, true);

  // data is the same for every point
  TH1* hData = dynamic_cast<TH1*>(points[0].hData->Clone("hDataFrame"));
  const double xMax = hData->GetBinLowEdge(hData->FindLastBinAbove(0) + 1) * 1.02;

  padTop->cd();
  setStyle(hData, kBlack);
  hData->SetMarkerSize(0.5);
  hData->GetXaxis()->SetRangeUser(0, xMax);
  hData->SetMinimum(0.2);
  hData->SetMaximum(hData->GetMaximum() * 100); // room for the legend
  styleTopFrame(hData, "Count");
  hData->Draw("pe");

  // data is drawn in black, so every fit gets a color
  ColorManager colorMgr(points.size());
  TLegend* leg = makeLegend(points.size() + 1, 0.4, 0.45, 0.95, 0.94);
  leg->AddEntry(hData, "Data", "p");
  for (size_t i = 0; i < points.size(); ++i) {
    setStyle(points[i].hFit, colorMgr.getColor(i));
    points[i].hFit->Draw("hist same");
    const TString chi2 = points[i].ndf > 0 ? Form(", #chi^{2}/ndf = %.1f", points[i].chi2 / points[i].ndf) : "";
    leg->AddEntry(points[i].hFit, Form("%s (#mu = %.2f, k = %.2f%s)", points[i].label.Data(), points[i].mu, points[i].k, chi2.Data()), "l");
  }
  hData->Draw("pe same");
  leg->Draw();

  padBot->cd();
  std::vector<TH1*> ratios;
  for (size_t i = 0; i < points.size(); ++i) {
    TH1* hRatio = dynamic_cast<TH1*>(points[i].hData->Clone(Form("hDataOverFit_%zu", i)));
    hRatio->Divide(points[i].hFit);
    setStyle(hRatio, colorMgr.getColor(i));
    hRatio->SetMarkerSize(0.4);
    hRatio->GetXaxis()->SetRangeUser(0, xMax);
    ratios.push_back(hRatio);
  }
  styleRatioFrame(ratios[0], "FT0C Amplitude", "Data / Fit", 0.45, 1.55);
  ratios[0]->Draw("p");
  drawUnityLine(ratios[0], xMax);
  for (TH1* hRatio : ratios) {
    hRatio->Draw("p same");
  }

  padTop->GetCanvas()->SaveAs(outFile);
}

void drawEntropy(const std::vector<ScanPoint>& points, const std::vector<int>& colors, const TString& outFile)
{
  TCanvas* canv = new TCanvas("canvEntropy", "", 1200, 800);
  canv->SetTicks(1, 1);
  canv->SetLogy();
  canv->SetLeftMargin(0.12);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.04);

  double xMax = 0, yMax = 0;
  for (const auto& point : points) {
    xMax = std::max(xMax, point.hEntropy->GetBinLowEdge(point.hEntropy->FindLastBinAbove(0) + 1));
    yMax = std::max(yMax, point.hEntropy->GetMaximum());
  }

  // common frame, the entropy range differs between the scan points
  TH1F* frame = canv->DrawFrame(0, yMax * 1e-6, xMax * 1.02, yMax * 5);
  frame->GetXaxis()->SetTitle("Entropy");
  frame->GetYaxis()->SetTitle("Probability density");

  TLegend* leg = makeLegend(points.size(), 0.6, 0.55, 0.95, 0.94);
  for (size_t i = 0; i < points.size(); ++i) {
    TH1* hist = points[i].hEntropy;
    setStyle(hist, colors[i]);
    hist->SetLineWidth(1);
    hist->Draw("hist same");
    leg->AddEntry(hist, points[i].label, "l");
  }
  leg->Draw();
  canv->SaveAs(outFile);
}

// <Npart> or <Ncoll> vs centrality, ratio to the reference point
void drawVsCentrality(const std::vector<ScanPoint>& points, const std::vector<int>& colors, size_t iRef,
                      TH1* ScanPoint::*member, const char* yTitle, const TString& outFile)
{
  TPad *padTop, *padBot;
  makePads(yTitle, padTop, padBot, false);

  padTop->cd();
  double yMax = 0;
  for (const auto& point : points) {
    yMax = std::max(yMax, (point.*member)->GetMaximum());
  }
  TLegend* leg = makeLegend(points.size(), 0.6, 0.45, 0.95, 0.94);
  for (size_t i = 0; i < points.size(); ++i) {
    TH1* hist = points[i].*member;
    setStyle(hist, colors[i]);
    if (i == 0) {
      hist->SetMinimum(0);
      hist->SetMaximum(yMax * 1.1);
      styleTopFrame(hist, yTitle);
      hist->Draw("p");
    } else {
      hist->Draw("p same");
    }
    leg->AddEntry(hist, Form("%s%s", points[i].label.Data(), i == iRef ? " (ref.)" : ""), "p");
  }
  leg->Draw();

  padBot->cd();
  std::vector<TH1*> ratios;
  double ratioMin = 1, ratioMax = 1;
  for (size_t i = 0; i < points.size(); ++i) {
    if (i == iRef) {
      continue; // the reference is the unity line
    }
    TH1* hRatio = dynamic_cast<TH1*>((points[i].*member)->Clone(Form("%s_ratio", (points[i].*member)->GetName())));
    hRatio->SetMinimum(); // drop the range inherited from the top pad
    hRatio->SetMaximum();
    hRatio->Divide(points[iRef].*member);
    ratioMin = std::min(ratioMin, hRatio->GetMinimum());
    ratioMax = std::max(ratioMax, hRatio->GetMaximum());
    ratios.push_back(hRatio);
  }
  const double margin = std::max(0.1 * (ratioMax - ratioMin), 0.01);
  TH1* frame = dynamic_cast<TH1*>((points[iRef].*member)->Clone(Form("%s_ratioFrame", (points[iRef].*member)->GetName())));
  frame->Reset();
  styleRatioFrame(frame, "Centrality (%)", "Ratio to ref.", ratioMin - margin, ratioMax + margin);
  frame->GetYaxis()->SetTitleSize(0.08);
  frame->GetYaxis()->SetTitleOffset(0.65);
  frame->Draw("axis");
  drawUnityLine(frame, frame->GetXaxis()->GetXmax());
  for (TH1* hRatio : ratios) {
    hRatio->Draw("p same");
  }

  padTop->GetCanvas()->SaveAs(outFile);
}
// One pad per row with a shared x axis: the deviation from the reference on top, then the value of every scan
// parameter. The x axis is the rank in the deviation, most accurate on the left; points right of the dashed
// line are not drawn in the other plots.
void drawParametersVsRank(const std::vector<std::pair<std::string, double>>& ranking, int nDraw, const TString& outFile)
{
  if (ranking.empty()) {
    return;
  }
  std::vector<std::vector<std::pair<std::string, double>>> pars;
  std::vector<std::string> parNames; // union over the scan points, in order of appearance
  for (const auto& entry : ranking) {
    pars.push_back(parseParameters(entry.first));
    for (const auto& par : pars.back()) {
      if (std::find(parNames.begin(), parNames.end(), par.first) == parNames.end()) {
        parNames.push_back(par.first);
      }
    }
  }

  const int nPoints = ranking.size();
  const int nRows = parNames.size() + 1;
  const double rowPx = 160, topPx = 20, bottomPx = 70;
  const double heightPx = nRows * rowPx + topPx + bottomPx;
  TCanvas* canv = new TCanvas("canvParameters", "", 1200, heightPx);
  ColorManager colorMgr(std::max<size_t>(parNames.size(), 1));

  for (int row = 0; row < nRows; ++row) {
    const bool first = row == 0, last = row == nRows - 1;
    const double upPx = first ? 0 : topPx + row * rowPx;
    const double lowPx = topPx + (row + 1) * rowPx + (last ? bottomPx : 0);
    canv->cd();
    TPad* pad = new TPad(Form("padParameters_%d", row), "", 0, 1 - lowPx / heightPx, 1, 1 - upPx / heightPx);
    pad->SetTicks(1, 1);
    pad->SetGridx();
    pad->SetLeftMargin(0.1);
    pad->SetRightMargin(0.03);
    pad->SetTopMargin(first ? topPx / (lowPx - upPx) : 0);
    pad->SetBottomMargin(last ? bottomPx / (lowPx - upPx) : 0);
    pad->Draw();
    pad->cd();

    TGraph* graph = new TGraph();
    for (int i = 0; i < nPoints; ++i) {
      if (first) {
        graph->SetPoint(graph->GetN(), i + 1, 100 * ranking[i].second);
        continue;
      }
      for (const auto& par : pars[i]) {
        if (par.first == parNames[row - 1]) {
          graph->SetPoint(graph->GetN(), i + 1, par.second);
        }
      }
    }
    double yMin = 0, yMax = 1;
    if (graph->GetN() > 0) {
      yMin = *std::min_element(graph->GetY(), graph->GetY() + graph->GetN());
      yMax = *std::max_element(graph->GetY(), graph->GetY() + graph->GetN());
    }
    const double margin = yMax > yMin ? 0.1 * (yMax - yMin) : std::max(0.1 * std::abs(yMax), 0.1);
    yMin -= margin;
    yMax += margin;

    TH1D* frame = new TH1D(Form("hParametersFrame_%d", row), "", nPoints, 0.5, nPoints + 0.5);
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    for (TAxis* axis : { frame->GetXaxis(), frame->GetYaxis() }) {
      axis->SetLabelFont(43);
      axis->SetLabelSize(18);
      axis->SetTitleFont(43);
      axis->SetTitleSize(20);
    }
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetTitleOffset(2.2);
    frame->GetYaxis()->SetTitle(first ? "Deviation (%)" : parNames[row - 1].c_str());
    frame->GetXaxis()->SetTitle("Rank (most accurate left)");
    frame->GetXaxis()->SetTitleOffset(3);
    if (!last) {
      frame->GetXaxis()->SetLabelSize(0);
    }
    frame->Draw("axis");

    if (nDraw > 0 && nDraw < nPoints) {
      TLine* line = new TLine(nDraw + 0.5, yMin, nDraw + 0.5, yMax);
      line->SetLineStyle(2);
      line->SetLineColor(kGray + 1);
      line->Draw();
    }
    const int color = first ? static_cast<int>(kBlack) : colorMgr.getColor(row - 1);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerSize(0.8);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->Draw("p");
  }
  canv->SaveAs(outFile);
}
} // namespace

/// @param outDir output directory of runScan.sh
/// @param filter regular expression on the scan point names, e.g. "^p_" or "norm_18.5_.*d_0.4"; empty = all
/// @param reference reference for the Npart/Ncoll ratios: either the name of a scan point, or a calibration
///        file (*.root, e.g. a Glauber or Trento calibration) containing h2dNpart and h2dNcoll; empty = first point
/// @param referenceLabel legend entry of a reference file
/// @param nDraw number of scan points closest to the reference (combined Npart/Ncoll metric) to draw; <= 0 = all
void drawScan(TString outDir = "output", TString filter = "", TString reference = "", TString referenceLabel = "Reference", int nDraw = 5)
{
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);
  reference = "data/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root";

  // Collect the scan points
  std::vector<std::string> names;
  TSystemDirectory dir(outDir, outDir);
  TList* entries = dir.GetListOfFiles();
  if (!entries) {
    std::cerr << "Error: cannot read directory " << outDir << std::endl;
    return;
  }
  TPRegexp regexp(filter);
  for (TObject* obj : *entries) {
    auto* entry = dynamic_cast<TSystemFile*>(obj);
    const TString name = entry->GetName();
    if (!entry->IsDirectory() || name == "." || name == "..") {
      continue;
    }
    if (!filter.IsNull() && !regexp.Match(name)) {
      continue;
    }
    if (gSystem->AccessPathName(outDir + "/" + name + "/trentofit_calib.root")) {
      std::cout << "Skipping " << name << ": no trentofit_calib.root" << std::endl;
      continue;
    }
    names.push_back(name.Data());
  }
  std::sort(names.begin(), names.end());

  std::vector<ScanPoint> points;
  for (const auto& name : names) {
    ScanPoint point;
    point.name = name;
    point.label = makeLabel(name);
    if (loadPoint(outDir + "/" + name.c_str(), point)) {
      points.push_back(point);
    } else {
      std::cout << "Skipping " << name << ": could not read its output" << std::endl;
    }
  }
  if (points.empty()) {
    std::cerr << "Error: no scan points found in " << outDir << (filter.IsNull() ? TString("") : " matching " + filter) << std::endl;
    return;
  }

  // Npart/Ncoll curves: the scan points, preceded by the reference if it comes from a file
  std::vector<ScanPoint> centPoints;
  size_t iRef = 0;
  const bool referenceFromFile = reference.EndsWith(".root");
  if (referenceFromFile) {
    ScanPoint ref;
    ref.name = "reference";
    ref.label = referenceLabel;
    if (!loadReference(reference, ref)) {
      return;
    }
    centPoints.push_back(ref);
  } else if (!reference.IsNull()) {
    auto it = std::find_if(points.begin(), points.end(), [&](const ScanPoint& p) { return p.name == reference.Data(); });
    if (it == points.end()) {
      std::cerr << "Warning: reference " << reference << " not found, using " << points[0].name << std::endl;
    } else {
      iRef = it - points.begin();
    }
  }
  centPoints.insert(centPoints.end(), points.begin(), points.end());

  // Summary table
  printf("\n%-40s %10s %10s %12s %10s %10s %10s %10s %10s\n", "Scan point", "mu", "k", "chi2/ndf", "fit min", "fit max", "anchor", "Npart 0-5", "Ncoll 0-5");
  for (const auto& point : points) {
    printf("%-40s %10.4f %10.4f %12.2f %10.1f %10.1f %10.1f %10.2f %10.2f\n", point.name.c_str(), point.mu, point.k,
           point.ndf > 0 ? point.chi2 / point.ndf : -1., point.fitMin, point.fitMax, point.anchorPoint,
           point.hNpart->GetBinContent(1), point.hNcoll->GetBinContent(1));
  }
  if (referenceFromFile) {
    printf("%-40s %10s %10s %12s %10s %10s %10s %10.2f %10.2f\n", referenceLabel.Data(), "", "", "", "", "", "",
           centPoints[0].hNpart->GetBinContent(1), centPoints[0].hNcoll->GetBinContent(1));
  }
  printf("(chi2/ndf = -1: not stored, rerun runTrentoFit.cc)\n");
  printf("\n");

  // Scan points closest to the reference in <Npart> and <Ncoll>
  struct Match {
    size_t index;
    double npart, ncoll, combined;
  };
  std::vector<Match> matches;
  std::vector<std::pair<std::string, double>> ranking; // scan point name and combined deviation, most accurate first
  for (size_t i = 0; i < centPoints.size(); ++i) {
    if (i == iRef) {
      continue;
    }
    const double npart = relativeRms(centPoints[i].hNpart, centPoints[iRef].hNpart);
    const double ncoll = relativeRms(centPoints[i].hNcoll, centPoints[iRef].hNcoll);
    matches.push_back({ i, npart, ncoll, std::sqrt(0.5 * (npart * npart + ncoll * ncoll)) });
  }
  if (!matches.empty()) {
    std::sort(matches.begin(), matches.end(), [](const Match& a, const Match& b) { return a.combined < b.combined; });
    for (const Match& m : matches) {
      ranking.emplace_back(centPoints[m.index].name, m.combined);
    }
    printf("Closest to %s (RMS of the relative deviation in centrality bins up to %.0f%%):\n", centPoints[iRef].label.Data(), CentMaxMatch);
    printf("  %-40s %10s %10s %10s\n", "Scan point", "Npart", "Ncoll", "combined");
    for (size_t i = 0; i < std::min(matches.size(), NClosestToPrint); ++i) {
      const Match& m = matches[i];
      printf("  %-40s %9.2f%% %9.2f%% %9.2f%%\n", centPoints[m.index].name.c_str(), 100 * m.npart, 100 * m.ncoll, 100 * m.combined);
    }
    auto bestNpart = std::min_element(matches.begin(), matches.end(), [](const Match& a, const Match& b) { return a.npart < b.npart; });
    auto bestNcoll = std::min_element(matches.begin(), matches.end(), [](const Match& a, const Match& b) { return a.ncoll < b.ncoll; });
    printf("  Best in Npart only: %s, best in Ncoll only: %s\n", centPoints[bestNpart->index].name.c_str(), centPoints[bestNcoll->index].name.c_str());
    const ScanPoint& best = centPoints[matches[0].index];
    printf("\nClosest scan point: %s (mu = %.4f, k = %.4f)\n  trento.conf:\n%s\n", best.name.c_str(), best.mu, best.k,
           best.trentoConfig.empty() ? "    (not found)\n" : best.trentoConfig.c_str());
  }

  // Only draw the reference and the nDraw scan points closest to it, ordered by the combined metric
  std::vector<size_t> selected = { iRef };
  for (size_t i = 0; i < matches.size() && (nDraw <= 0 || i < static_cast<size_t>(nDraw)); ++i) {
    selected.push_back(matches[i].index);
  }
  std::vector<ScanPoint> drawCentPoints;
  for (size_t i : selected) {
    drawCentPoints.push_back(centPoints[i]);
  }
  centPoints = drawCentPoints;
  iRef = 0;
  points.assign(centPoints.begin() + (referenceFromFile ? 1 : 0), centPoints.end()); // scan points only

  // Reference in black, the others from the ColorManager
  ColorManager colorMgr(std::max<size_t>(centPoints.size() - 1, 1));
  std::vector<int> centColors;
  for (size_t i = 0, iColor = 0; i < centPoints.size(); ++i) {
    centColors.push_back(i == iRef ? static_cast<int>(kBlack) : colorMgr.getColor(iColor++));
  }
  const std::vector<int> colors(centColors.begin() + (referenceFromFile ? 1 : 0), centColors.end()); // scan points only

  // file name tag from the filter, keeping only file name friendly characters
  TString tag = filter;
  TPRegexp("[^A-Za-z0-9.-]+").Substitute(tag, "_", "g");
  tag = tag.Strip(TString::kBoth, '_');
  tag = tag.IsNull() ? "" : "_" + tag;
  drawFits(points, outDir + "/scan" + tag + "_Fits.pdf");
  drawEntropy(points, colors, outDir + "/scan" + tag + "_Entropy.pdf");
  drawVsCentrality(centPoints, centColors, iRef, &ScanPoint::hNpart, "<N_{part}>", outDir + "/scan" + tag + "_Npart.pdf");
  drawVsCentrality(centPoints, centColors, iRef, &ScanPoint::hNcoll, "<N_{coll}>", outDir + "/scan" + tag + "_Ncoll.pdf");
  drawParametersVsRank(ranking, nDraw, "scan" + tag + "_Parameters.pdf");
}
