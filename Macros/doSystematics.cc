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
  GlauberParameters baseFT0C("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");
  std::vector<std::vector<float>> systFT0C = { NpartErr, NcollErr };
  std::vector<std::vector<float>> relSystFT0C = { NpartRelErr, NcollRelErr };

  GlauberParameters bcs("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_LightIonDef.root");
  std::vector<std::vector<float>> systBcs = computeSystematics(base, { bcs });
  std::vector<std::vector<float>> relSystBcs = computeRelativeSystematics(base, { bcs });

  GlauberParameters anchor50("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor50.root");
  GlauberParameters anchor60("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor60.root");
  GlauberParameters anchor70("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor70.root");
  GlauberParameters anchor80("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor80.root");
  GlauberParameters anchor90("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor90.root");
  std::vector<std::vector<float>> systAnchor = computeSystematics(base, { anchor50, anchor60, anchor70, anchor80, anchor90 });
  std::vector<std::vector<float>> relSystAnchor = computeRelativeSystematics(base, { anchor50, anchor60, anchor70, anchor80, anchor90 });

  GlauberParameters run564356("LHC25ae_pass2/AR_564356_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564359("LHC25ae_pass2/AR_564359_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564373("LHC25ae_pass2/AR_564373_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564387("LHC25ae_pass2/AR_564387_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564400("LHC25ae_pass2/AR_564400_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564414("LHC25ae_pass2/AR_564414_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564430("LHC25ae_pass2/AR_564430_calibration_hFT0M_Collisions_LightIonDef.root");
  GlauberParameters run564445("LHC25ae_pass2/AR_564445_calibration_hFT0M_Collisions_LightIonDef.root");
  std::vector<std::vector<float>> systRun = computeSystematics(base, { run564356 });
  std::vector<std::vector<float>> relSystRun = computeRelativeSystematics(base, { run564356 });
  // std::vector<std::vector<float>> systRun = computeSystematics(base, { run564356, run564359, run564373, run564387, run564400, run564414, run564430, run564445 });
  // std::vector<std::vector<float>> relSystRun = computeRelativeSystematics(base, { run564356, run564359, run564373, run564387, run564400, run564414, run564430, run564445 });

  std::vector<std::vector<float>> systTotal = combineSystematicsInQuadrature(systAnchor, systBcs, systRun, systFT0C);
  std::vector<std::vector<float>> relSystTotal = combineSystematicsInQuadrature(relSystAnchor, relSystBcs, relSystRun, relSystFT0C);

  ColorManager cm1(4);
  TH1F* hSysShapeNpart = initGlauParHist(systFT0C[kNpart], "hSysShapeNpart", cm1.getColor(0));
  TH1F* hSysShapeNcoll = initGlauParHist(systFT0C[kNcoll], "hSysShapeNcoll", cm1.getColor(0));
  TH1F* hRelSysShapeNpart = initGlauParHist(relSystFT0C[kNpart], "hRelSysShapeNpart", cm1.getColor(0));
  TH1F* hRelSysShapeNcoll = initGlauParHist(relSystFT0C[kNcoll], "hRelSysShapeNcoll", cm1.getColor(0));

  TH1F* hSysBcsNpart = initGlauParHist(systBcs[kNpart], "hSysBcsNpart", cm1.getColor(1));
  TH1F* hSysBcsNcoll = initGlauParHist(systBcs[kNcoll], "hSysBcsNcoll", cm1.getColor(1));
  TH1F* hRelSysBcsNpart = initGlauParHist(relSystBcs[kNpart], "hRelSysBcsNpart", cm1.getColor(1));
  TH1F* hRelSysBcsNcoll = initGlauParHist(relSystBcs[kNcoll], "hRelSysBcsNcoll", cm1.getColor(1));

  TH1F* hSysAnchorNpart = initGlauParHist(systAnchor[kNpart], "hSysAnchorNpart", cm1.getColor(2));
  TH1F* hSysAnchorNcoll = initGlauParHist(systAnchor[kNcoll], "hSysAnchorNcoll", cm1.getColor(2));
  TH1F* hRelSysAnchorNpart = initGlauParHist(relSystAnchor[kNpart], "hRelSysAnchorNpart", cm1.getColor(2));
  TH1F* hRelSysAnchorNcoll = initGlauParHist(relSystAnchor[kNcoll], "hRelSysAnchorNcoll", cm1.getColor(2));

  TH1F* hSysRunNpart = initGlauParHist(systRun[kNpart], "hSysRunNpart", cm1.getColor(3));
  TH1F* hSysRunNcoll = initGlauParHist(systRun[kNcoll], "hSysRunNcoll", cm1.getColor(3));
  TH1F* hRelSysRunNpart = initGlauParHist(relSystRun[kNpart], "hRelSysRunNpart", cm1.getColor(3));
  TH1F* hRelSysRunNcoll = initGlauParHist(relSystRun[kNcoll], "hRelSysRunNcoll", cm1.getColor(3));

  TH1F* hSysTotalNpart = initGlauParHist(systTotal[kNpart], "hSysTotalNpart", kBlack);
  TH1F* hSysTotalNcoll = initGlauParHist(systTotal[kNcoll], "hSysTotalNcoll", kBlack);
  TH1F* hRelSysTotalNpart = initGlauParHist(relSystTotal[kNpart], "hRelSysTotalNpart", kBlack);
  TH1F* hRelSysTotalNcoll = initGlauParHist(relSystTotal[kNcoll], "hRelSysTotalNcoll", kBlack);

  TH1F* hNpartFT0C = initGlauParHist(NpartVal, "hNpartFT0C", kBlack);
  TH1F* hNcollFT0C = initGlauParHist(NcollVal, "hNpartFT0C", kBlack);


  TLegend* leg = new TLegend(0.82, 0.36, 1.12, 0.64);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hSysTotalNpart, "Total", "l");
  leg->AddEntry(hSysShapeNpart, "Shape", "l");
  leg->AddEntry(hSysAnchorNpart, "Anchor", "l");
  leg->AddEntry(hSysRunNpart, "Run", "l");
  leg->AddEntry(hSysBcsNpart, "Bcs", "l");

  const bool drawShape = true;
  const bool drawAnchor = true;
  const bool drawRun = true;
  const bool drawBcs = true;
  const bool drawTotal = true;

  TH1F* hSysTotalNpartCopy = dynamic_cast<TH1F*>(hSysTotalNpart->Clone("hSysTotalNpartCopy"));
  hSysTotalNpartCopy->SetMarkerColorAlpha(0, 0);
  hSysTotalNpartCopy->SetLineColorAlpha(0, 0);
  TH1F* hSysTotalNcollCopy = dynamic_cast<TH1F*>(hSysTotalNcoll->Clone("hSysTotalNcollCopy"));
  hSysTotalNcollCopy->SetMarkerColorAlpha(0, 0);
  hSysTotalNcollCopy->SetLineColorAlpha(0, 0);
  TH1F* hRelSysTotalNpartCopy = dynamic_cast<TH1F*>(hRelSysTotalNpart->Clone("hRelSysTotalNpartCopy"));
  hRelSysTotalNpartCopy->SetMarkerColorAlpha(0, 0);
  hRelSysTotalNpartCopy->SetLineColorAlpha(0, 0);
  TH1F* hRelSysTotalNcollCopy = dynamic_cast<TH1F*>(hRelSysTotalNcoll->Clone("hRelSysTotalNcollCopy"));
  hRelSysTotalNcollCopy->SetMarkerColorAlpha(0, 0);
  hRelSysTotalNcollCopy->SetLineColorAlpha(0, 0);

  const int canvWidth = 1600;
  const int canvHeight = 1000;
  
  TCanvas* canvNpartSystematics = new TCanvas("canvNpartSystematics", "", canvWidth, canvHeight);
  styleCanvas(canvNpartSystematics);
  styleHist(hSysTotalNpartCopy, "<N_{part}> systematic uncertainty");
  hSysTotalNpartCopy->Draw("");
  if (drawShape) hSysShapeNpart->Draw("hist same");
  if (drawAnchor) hSysAnchorNpart->Draw("hist same");
  if (drawRun) hSysRunNpart->Draw("hist same");
  if (drawBcs) hSysBcsNpart->Draw("hist same");
  if (drawTotal) hSysTotalNpart->Draw("hist same");
  leg->Draw();
  canvNpartSystematics->SaveAs("hNpartSystematics.pdf");

  TCanvas* canvNcollSystematics = new TCanvas("canvNcollSystematics", "", canvWidth, canvHeight);
  styleCanvas(canvNcollSystematics);
  styleHist(hSysTotalNcollCopy, "<N_{coll}> systematic uncertainty");
  hSysTotalNcollCopy->Draw("");
  if (drawShape) hSysShapeNcoll->Draw("hist same");
  if (drawAnchor) hSysAnchorNcoll->Draw("hist same");
  if (drawRun) hSysRunNcoll->Draw("hist same");
  if (drawBcs) hSysBcsNcoll->Draw("hist same");
  if (drawTotal) hSysTotalNcoll->Draw("hist same");
  leg->Draw();
  canvNcollSystematics->SaveAs("hNcollSystematics.pdf");

  TCanvas* canvRelNpartSystematics = new TCanvas("canvRelNpartSystematics", "", canvWidth, canvHeight);
  styleCanvas(canvRelNpartSystematics);
  styleHist(hRelSysTotalNpartCopy, "<N_{part}> relative systematic uncertainty");
  hRelSysTotalNpartCopy->Draw("");
  if (drawShape) hRelSysShapeNpart->Draw("hist same");
  if (drawAnchor) hRelSysAnchorNpart->Draw("hist same");
  if (drawRun) hRelSysRunNpart->Draw("hist same");
  if (drawBcs) hRelSysBcsNpart->Draw("hist same");
  if (drawTotal) hRelSysTotalNpart->Draw("hist same");
  leg->Draw();
  canvRelNpartSystematics->SaveAs("hNpartRelSystematics.pdf");

  TCanvas* canvRelNcollSystematics = new TCanvas("canvRelNcollSystematics", "", canvWidth, canvHeight);
  styleCanvas(canvRelNcollSystematics);
  styleHist(hRelSysTotalNcollCopy, "<N_{coll}> relative systematic uncertainty");
  hRelSysTotalNcollCopy->Draw("");
  if (drawShape) hRelSysShapeNcoll->Draw("hist same");
  if (drawAnchor) hRelSysAnchorNcoll->Draw("hist same");
  if (drawRun) hRelSysRunNcoll->Draw("hist same");
  if (drawBcs) hRelSysBcsNcoll->Draw("hist same");
  if (drawTotal) hRelSysTotalNcoll->Draw("hist same");
  leg->Draw();
  canvRelNcollSystematics->SaveAs("hNcollRelSystematics.pdf");


  /* Extra QA */
  ColorManager cm2(7);
  TH1F* hNpartFT0M = initGlauParHist(base.npart, "hNpartFT0M", kBlue);
  TH1F* hNpartAnchor50 = initGlauParHist(anchor50.npart, "hNpartAnchor50", cm2.getColor(0));
  TH1F* hNpartAnchor60 = initGlauParHist(anchor60.npart, "hNpartAnchor60", cm2.getColor(1));
  TH1F* hNpartAnchor70 = initGlauParHist(anchor70.npart, "hNpartAnchor70", cm2.getColor(2));
  TH1F* hNpartAnchor80 = initGlauParHist(anchor80.npart, "hNpartAnchor80", cm2.getColor(3));
  TH1F* hNpartAnchor90 = initGlauParHist(anchor90.npart, "hNpartAnchor90", cm2.getColor(4));
  TH1F* hNpartBcs = initGlauParHist(bcs.npart, "hNpartBcs", cm2.getColor(5));
  TH1F* hNpartRun = initGlauParHist(bcs.npart, "hNpartRun", cm2.getColor(6));

  TH1F* hNcollFT0M = initGlauParHist(base.ncoll, "hNcollFT0M", kBlue);
  TH1F* hNcollAnchor50 = initGlauParHist(anchor50.ncoll, "hNcollAnchor50", cm2.getColor(0));
  TH1F* hNcollAnchor60 = initGlauParHist(anchor60.ncoll, "hNcollAnchor60", cm2.getColor(1));
  TH1F* hNcollAnchor70 = initGlauParHist(anchor70.ncoll, "hNcollAnchor70", cm2.getColor(2));
  TH1F* hNcollAnchor80 = initGlauParHist(anchor80.ncoll, "hNcollAnchor80", cm2.getColor(3));
  TH1F* hNcollAnchor90 = initGlauParHist(anchor90.ncoll, "hNcollAnchor90", cm2.getColor(4));
  TH1F* hNcollBcs = initGlauParHist(bcs.ncoll, "hNcollBcs", cm2.getColor(5));
  TH1F* hNcollRun = initGlauParHist(bcs.ncoll, "hNcollRun", cm2.getColor(6));

  doGlauParQA(base.h2dNpart, hNpartFT0M, "BaseNpartFT0M");
  doGlauParQA(base.h2dNcoll, hNcollFT0M, "BaseNcollFT0M");

  std::vector<float> errNpart;
  std::vector<float> errNcoll;

  const bool publishValues = true;
  for (int ii = 0; ii < CentBins.size() - 1; ++ii) {
    errNpart.push_back(NpartVal[ii] * relSystTotal[kNpart][ii]);
    errNcoll.push_back(NcollVal[ii] * relSystTotal[kNcoll][ii]);
    if (publishValues) {
      std::cout << Form("%.f-%.f%% | Npart: %.2f +- %.2f | Ncoll: %.2f +- %.2f", CentBins[ii], CentBins[ii + 1], NpartVal[ii], errNpart[ii], NcollVal[ii], errNcoll[ii]) << std::endl;
    }
  }

  TH1F* hOldNpartErr = initGlauParHist(NpartErr, "hOldNpartErr", kBlack);
  TH1F* hNewNpartErr = initGlauParHist(errNpart, "hNewNpartErr", kBlue);
  TH1F* hOldNcollErr = initGlauParHist(NcollErr, "hOldNcollErr", kBlack);
  TH1F* hNewNcollErr = initGlauParHist(errNcoll, "hNewNcollErr", kBlue);

  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 1000);
  styleCanvas(canvNpart);
  canvNpart->SetRightMargin(0.03);
  styleHist(hOldNpartErr, "Uncertainty");
  hOldNpartErr->SetMaximum(1.2);
  hOldNpartErr->SetMinimum(0);
  hOldNpartErr->Draw("hist");
  hNewNpartErr->Draw("hist same");
  canvNpart->SaveAs("qaNpart.pdf");
  
  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 1000);
  styleCanvas(canvNcoll);
  canvNcoll->SetRightMargin(0.03);
  styleHist(hOldNcollErr, "Uncertainty");
  hOldNcollErr->SetMaximum(3.3);
  hOldNcollErr->SetMinimum(0);
  hOldNcollErr->Draw("hist");
  hNewNcollErr->Draw("hist same");
  canvNcoll->SaveAs("qaNcoll.pdf");


}