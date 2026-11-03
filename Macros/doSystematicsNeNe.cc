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

#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"
#include "Includes/saveHelper.h"

//  NeNe Npart            NeNe Ncoll
//  32.10    0.70         53.86    4.82
//  29.39    1.04         45.18    4.42
//  25.35    1.07         35.40    3.49
//  19.84    0.74         24.36    2.21
//  14.97    0.51         16.26    1.48
//  11.06    0.48         10.73    1.08
//  8.07     0.43         7.07     0.73
//  5.88     0.30         4.68     0.43
//  4.34     0.15         3.14     0.22
//  3.26     0.15         2.13     0.18
// //  2.45     0.14         1.41     0.14
// Published values per centrality class (CentClasses) and for 0-100% (MB)
static const std::vector<float> CentClasses = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
static const std::vector<float> NpartValCent = { 32.10, 29.39, 25.35, 19.84, 14.97, 11.06, 8.07, 5.88, 4.34, 3.26, 2.45 };
static const std::vector<float> NpartErrCent = { 0.70, 1.04, 1.07, 0.74, 0.51, 0.48, 0.43, 0.30, 0.15, 0.15, 0.14 };
static const std::vector<float> NcollValCent = { 53.86, 45.18, 35.40, 24.36, 16.26, 10.73, 7.07, 4.68, 3.14, 2.13, 1.41 };
static const std::vector<float> NcollErrCent = { 4.82, 4.42, 3.49, 2.21, 1.48, 1.08, 0.73, 0.43, 0.22, 0.18, 0.14 };
static const std::vector<float> NpartValMB = { 12.60 };
static const std::vector<float> NpartErrMB = { 0.54 };
static const std::vector<float> NcollValMB = { 15.47 };
static const std::vector<float> NcollErrMB = { 1.43 };

// err / val, bin by bin
static std::vector<float> relativeErrors(const std::vector<float>& errs, const std::vector<float>& vals)
{
  std::vector<float> rel;
  for (size_t ii = 0; ii < vals.size(); ++ii) {
    rel.push_back(errs[ii] / vals[ii]);
  }
  return rel;
}


void doSystematicsNeNe()
{
  gStyle->SetOptStat(0);
  // Centrality binning: false = centrality classes, true = a single 0-100% bin
  const bool minBias = true;
  syst::CentBins = minBias ? std::vector<float>{ 0, 100 } : CentClasses;
  const std::vector<float>& NpartVal = minBias ? NpartValMB : NpartValCent;
  const std::vector<float>& NpartErr = minBias ? NpartErrMB : NpartErrCent;
  const std::vector<float>& NcollVal = minBias ? NcollValMB : NcollValCent;
  const std::vector<float>& NcollErr = minBias ? NcollErrMB : NcollErrCent;
  const std::vector<float> NpartRelErr = relativeErrors(NpartErr, NpartVal);
  const std::vector<float> NcollRelErr = relativeErrors(NcollErr, NcollVal);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  syst::GlauberParameters base("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters baseFT0C("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_LightIonDef.root");
  std::vector<std::vector<float>> systFT0C = { NpartErr, NcollErr };
  std::vector<std::vector<float>> relSystFT0C = { NpartRelErr, NcollRelErr };

  syst::GlauberParameters bcAnchor50("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor50.root");
  syst::GlauberParameters bcAnchor60("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor60.root");
  syst::GlauberParameters bcAnchor70("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor70.root");
  syst::GlauberParameters bcAnchor80("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor80.root");
  syst::GlauberParameters bcAnchor90("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor90.root");
  // syst::GlauberParameters bcAnchor70("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_LightIonDef.root");
  std::vector<std::vector<float>> systBcs = syst::computeSystematics(base, { bcAnchor70 });
  std::vector<std::vector<float>> relSystBcs = syst::computeRelativeSystematics(base, { bcAnchor70 });
  // std::vector<std::vector<float>> systBcs = syst::computeSystematics(base, { bcAnchor50, bcAnchor60, bcAnchor70, bcAnchor80, bcAnchor90 });
  // std::vector<std::vector<float>> relSystBcs = syst::computeRelativeSystematics(base, { bcAnchor50, bcAnchor60, bcAnchor70, bcAnchor80, bcAnchor90 });
  syst::GlauberParameters colAnchor50("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor50.root");
  syst::GlauberParameters colAnchor60("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor60.root");
  syst::GlauberParameters colAnchor70("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor70.root");
  syst::GlauberParameters colAnchor80("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor80.root");
  syst::GlauberParameters colAnchor90("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor90.root");
  std::vector<std::vector<float>> systAnchor = syst::computeSystematics(base, { colAnchor70 });
  std::vector<std::vector<float>> relSystAnchor = syst::computeRelativeSystematics(base, { colAnchor70 });
  // std::vector<std::vector<float>> systAnchor = syst::computeSystematics(base, { colAnchor50, colAnchor60, colAnchor70, colAnchor80, colAnchor90 });
  // std::vector<std::vector<float>> relSystAnchor = syst::computeRelativeSystematics(base, { colAnchor50, colAnchor60, colAnchor70, colAnchor80, colAnchor90 });

  syst::GlauberParameters run564472("LHC25af_pass2_systematics/AR_564472_calibration_hFT0M_Collisions_LightIonDef.root");
  std::vector<std::vector<float>> systRun = syst::computeSystematics(base, { run564472 });
  std::vector<std::vector<float>> relSystRun = syst::computeRelativeSystematics(base, { run564472 });

  std::vector<std::vector<float>> systTotal = syst::combineSystematicsInQuadrature(systAnchor, systBcs, systRun, systFT0C);
  std::vector<std::vector<float>> relSystTotal = syst::combineSystematicsInQuadrature(relSystAnchor, relSystBcs, relSystRun, relSystFT0C);

  ColorManager cm1(4);
  TH1F* hSysShapeNpart = syst::initGlauParHist(systFT0C[kNpart], "hSysShapeNpart", cm1.getColor(0));
  TH1F* hSysShapeNcoll = syst::initGlauParHist(systFT0C[kNcoll], "hSysShapeNcoll", cm1.getColor(0));
  TH1F* hRelSysShapeNpart = syst::initGlauParHist(relSystFT0C[kNpart], "hRelSysShapeNpart", cm1.getColor(0));
  TH1F* hRelSysShapeNcoll = syst::initGlauParHist(relSystFT0C[kNcoll], "hRelSysShapeNcoll", cm1.getColor(0));

  TH1F* hSysBcsNpart = syst::initGlauParHist(systBcs[kNpart], "hSysBcsNpart", cm1.getColor(1));
  TH1F* hSysBcsNcoll = syst::initGlauParHist(systBcs[kNcoll], "hSysBcsNcoll", cm1.getColor(1));
  TH1F* hRelSysBcsNpart = syst::initGlauParHist(relSystBcs[kNpart], "hRelSysBcsNpart", cm1.getColor(1));
  TH1F* hRelSysBcsNcoll = syst::initGlauParHist(relSystBcs[kNcoll], "hRelSysBcsNcoll", cm1.getColor(1));

  TH1F* hSysAnchorNpart = syst::initGlauParHist(systAnchor[kNpart], "hSysAnchorNpart", cm1.getColor(2));
  TH1F* hSysAnchorNcoll = syst::initGlauParHist(systAnchor[kNcoll], "hSysAnchorNcoll", cm1.getColor(2));
  TH1F* hRelSysAnchorNpart = syst::initGlauParHist(relSystAnchor[kNpart], "hRelSysAnchorNpart", cm1.getColor(2));
  TH1F* hRelSysAnchorNcoll = syst::initGlauParHist(relSystAnchor[kNcoll], "hRelSysAnchorNcoll", cm1.getColor(2));

  TH1F* hSysRunNpart = syst::initGlauParHist(systRun[kNpart], "hSysRunNpart", cm1.getColor(3));
  TH1F* hSysRunNcoll = syst::initGlauParHist(systRun[kNcoll], "hSysRunNcoll", cm1.getColor(3));
  TH1F* hRelSysRunNpart = syst::initGlauParHist(relSystRun[kNpart], "hRelSysRunNpart", cm1.getColor(3));
  TH1F* hRelSysRunNcoll = syst::initGlauParHist(relSystRun[kNcoll], "hRelSysRunNcoll", cm1.getColor(3));

  TH1F* hSysTotalNpart = syst::initGlauParHist(systTotal[kNpart], "hSysTotalNpart", kBlack);
  TH1F* hSysTotalNcoll = syst::initGlauParHist(systTotal[kNcoll], "hSysTotalNcoll", kBlack);
  TH1F* hRelSysTotalNpart = syst::initGlauParHist(relSystTotal[kNpart], "hRelSysTotalNpart", kBlack);
  TH1F* hRelSysTotalNcoll = syst::initGlauParHist(relSystTotal[kNcoll], "hRelSysTotalNcoll", kBlack);

  TH1F* hNpartFT0C = syst::initGlauParHist(NpartVal, "hNpartFT0C", kBlack);
  TH1F* hNcollFT0C = syst::initGlauParHist(NcollVal, "hNpartFT0C", kBlack);


  TLegend* leg = new TLegend(0.18, 0.66, 0.43, 0.94);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hSysTotalNpart, "Total", "l");
  leg->AddEntry(hSysShapeNpart, "Nuclei shape", "l");
  leg->AddEntry(hSysAnchorNpart, "Coll anchor 70%", "l");
  leg->AddEntry(hSysBcsNpart, "Bc anchor 70%", "l");
  leg->AddEntry(hSysRunNpart, "Run", "l");

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
  
  hRelSysTotalNpartCopy->SetMaximum(0.12);
  hRelSysTotalNcollCopy->SetMaximum(0.22);
  hRelSysTotalNpartCopy->SetMinimum(0);
  hRelSysTotalNcollCopy->SetMinimum(0);

  TCanvas* canvNpartSystematics = new TCanvas("canvNpartSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvNpartSystematics);
  syst::styleHist(hSysTotalNpartCopy, "<N_{part}> systematic uncertainty");
  hSysTotalNpartCopy->Draw("");
  if (drawShape) {
    hSysShapeNpart->Draw("hist same");
  }
  if (drawAnchor) {
    hSysAnchorNpart->Draw("hist same");
  }
  if (drawRun) {
    hSysRunNpart->Draw("hist same");
  }
  if (drawBcs) {
    hSysBcsNpart->Draw("hist same");
  }
  if (drawTotal) {
    hSysTotalNpart->Draw("hist same");
  }
  leg->Draw();
  canvNpartSystematics->SaveAs("hNpartSystematics_NeNe.pdf");
  save::canvasToRoot(canvNpartSystematics, "hNpartSystematics_NeNe.pdf");

  TCanvas* canvNcollSystematics = new TCanvas("canvNcollSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvNcollSystematics);
  syst::styleHist(hSysTotalNcollCopy, "<N_{coll}> systematic uncertainty");
  hSysTotalNcollCopy->Draw("");
  if (drawShape) {
    hSysShapeNcoll->Draw("hist same");
  }
  if (drawAnchor) {
    hSysAnchorNcoll->Draw("hist same");
  }
  if (drawRun) {
    hSysRunNcoll->Draw("hist same");
  }
  if (drawBcs) {
    hSysBcsNcoll->Draw("hist same");
  }
  if (drawTotal) {
    hSysTotalNcoll->Draw("hist same");
  }
  leg->Draw();
  canvNcollSystematics->SaveAs("hNcollSystematics_NeNe.pdf");
  save::canvasToRoot(canvNcollSystematics, "hNcollSystematics_NeNe.pdf");

  TCanvas* canvRelNpartSystematics = new TCanvas("canvRelNpartSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvRelNpartSystematics);
  syst::styleHist(hRelSysTotalNpartCopy, "<N_{part}> relative systematic uncertainty");
  hRelSysTotalNpartCopy->Draw("");
  if (drawShape) {
    hRelSysShapeNpart->Draw("hist same");
  }
  if (drawAnchor) {
    hRelSysAnchorNpart->Draw("hist same");
  }
  if (drawRun) {
    hRelSysRunNpart->Draw("hist same");
  }
  if (drawBcs) {
    hRelSysBcsNpart->Draw("hist same");
  }
  if (drawTotal) {
    hRelSysTotalNpart->Draw("hist same");
  }
  leg->Draw();
  canvRelNpartSystematics->SaveAs("FiguresSystematics/hNpartRelSystematics_NeNe.pdf");
  save::canvasToRoot(canvRelNpartSystematics, "FiguresSystematics/hNpartRelSystematics_NeNe.pdf");

  TCanvas* canvRelNcollSystematics = new TCanvas("canvRelNcollSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvRelNcollSystematics);
  syst::styleHist(hRelSysTotalNcollCopy, "<N_{coll}> relative systematic uncertainty");
  hRelSysTotalNcollCopy->Draw("");
  if (drawShape) {
    hRelSysShapeNcoll->Draw("hist same");
  }
  if (drawAnchor) {
    hRelSysAnchorNcoll->Draw("hist same");
  }
  if (drawRun) {
    hRelSysRunNcoll->Draw("hist same");
  }
  if (drawBcs) {
    hRelSysBcsNcoll->Draw("hist same");
  }
  if (drawTotal) {
    hRelSysTotalNcoll->Draw("hist same");
  }
  leg->Draw();
  canvRelNcollSystematics->SaveAs("FiguresSystematics/hNcollRelSystematics_NeNe.pdf");
  save::canvasToRoot(canvRelNcollSystematics, "FiguresSystematics/hNcollRelSystematics_NeNe.pdf");


  /* Extra QA */
  ColorManager cm2(7);
  TH1F* hNpartFT0M = syst::initGlauParHist(base.npart, "hNpartFT0M", kBlue);
  TH1F* hNpartAnchor50 = syst::initGlauParHist(colAnchor50.npart, "hNpartAnchor50", cm2.getColor(0));
  TH1F* hNpartAnchor60 = syst::initGlauParHist(colAnchor60.npart, "hNpartAnchor60", cm2.getColor(1));
  TH1F* hNpartAnchor70 = syst::initGlauParHist(colAnchor70.npart, "hNpartAnchor70", cm2.getColor(2));
  TH1F* hNpartAnchor80 = syst::initGlauParHist(colAnchor80.npart, "hNpartAnchor80", cm2.getColor(3));
  TH1F* hNpartAnchor90 = syst::initGlauParHist(colAnchor90.npart, "hNpartAnchor90", cm2.getColor(4));
  TH1F* hNpartBcs = syst::initGlauParHist(bcAnchor70.npart, "hNpartBcs", cm2.getColor(5));
  TH1F* hNpartRun = syst::initGlauParHist(bcAnchor70.npart, "hNpartRun", cm2.getColor(6));

  TH1F* hNcollFT0M = syst::initGlauParHist(base.ncoll, "hNcollFT0M", kBlue);
  TH1F* hNcollAnchor50 = syst::initGlauParHist(colAnchor50.ncoll, "hNcollAnchor50", cm2.getColor(0));
  TH1F* hNcollAnchor60 = syst::initGlauParHist(colAnchor60.ncoll, "hNcollAnchor60", cm2.getColor(1));
  TH1F* hNcollAnchor70 = syst::initGlauParHist(colAnchor70.ncoll, "hNcollAnchor70", cm2.getColor(2));
  TH1F* hNcollAnchor80 = syst::initGlauParHist(colAnchor80.ncoll, "hNcollAnchor80", cm2.getColor(3));
  TH1F* hNcollAnchor90 = syst::initGlauParHist(colAnchor90.ncoll, "hNcollAnchor90", cm2.getColor(4));
  TH1F* hNcollBcs = syst::initGlauParHist(bcAnchor70.ncoll, "hNcollBcs", cm2.getColor(5));
  TH1F* hNcollRun = syst::initGlauParHist(bcAnchor70.ncoll, "hNcollRun", cm2.getColor(6));

  syst::doGlauParQA(base.h2dNpart, hNpartFT0M, "BaseNpartFT0M");
  syst::doGlauParQA(base.h2dNcoll, hNcollFT0M, "BaseNcollFT0M");

  std::vector<float> errNpart;
  std::vector<float> errNcoll;

  const bool publishValues = true;
  const bool publishCrossCheck = false;
  for (int ii = 0; ii < syst::CentBins.size() - 1; ++ii) {
    const float errNpartFT0C = NpartVal[ii] * relSystTotal[kNpart][ii];
    const float errNcollFT0C = NcollVal[ii] * relSystTotal[kNcoll][ii];
    // RMS: assume the relative RMS of our calibration (base) also holds for the published values
    const float rmsNpartFT0C = NpartVal[ii] * base.npartRMS[ii] / base.npart[ii];
    const float rmsNcollFT0C = NcollVal[ii] * base.ncollRMS[ii] / base.ncoll[ii];
    const float errNpartFT0M = systTotal[kNpart][ii];
    const float errNcollFT0M = systTotal[kNcoll][ii];

    errNpart.push_back(errNpartFT0C);
    errNcoll.push_back(errNcollFT0C);
    if (publishValues) {
      std::cout << Form("%s | Npart: %.2f (RMS %.2f) +- %.2f | Ncoll: %.2f (RMS %.2f) +- %.2f", syst::binLabel(syst::CentBins, ii).Data(), NpartVal[ii], rmsNpartFT0C, errNpartFT0C, NcollVal[ii], rmsNcollFT0C, errNcollFT0C) << std::endl;
    }
    if (publishCrossCheck) {
      std::cout << Form("%s | Npart: %.4f +- %.4f | Ncoll: %.4f +- %.4f", syst::binLabel(syst::CentBins, ii).Data(), errNpartFT0M, systTotal[kNpart][ii], errNcollFT0M, systTotal[kNcoll][ii]) << std::endl;
    }
  }

  TH1F* hOldNpartErr = syst::initGlauParHist(NpartErr, "hOldNpartErr", kBlack);
  TH1F* hNewNpartErr = syst::initGlauParHist(errNpart, "hNewNpartErr", kBlue);
  TH1F* hOldNcollErr = syst::initGlauParHist(NcollErr, "hOldNcollErr", kBlack);
  TH1F* hNewNcollErr = syst::initGlauParHist(errNcoll, "hNewNcollErr", kBlue);

  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 1000);
  syst::styleCanvas(canvNpart);
  canvNpart->SetRightMargin(0.03);
  syst::styleHist(hOldNpartErr, "Uncertainty");
  hOldNpartErr->SetMaximum(1.2);
  hOldNpartErr->SetMinimum(0);
  hOldNpartErr->Draw("hist");
  hNewNpartErr->Draw("hist same");
  canvNpart->SaveAs("qaNpart_NeNe.pdf");
  save::canvasToRoot(canvNpart, "qaNpart_NeNe.pdf");
  
  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 1000);
  syst::styleCanvas(canvNcoll);
  canvNcoll->SetRightMargin(0.03);
  syst::styleHist(hOldNcollErr, "Uncertainty");
  hOldNcollErr->SetMaximum(3.3);
  hOldNcollErr->SetMinimum(0);
  hOldNcollErr->Draw("hist");
  hNewNcollErr->Draw("hist same");
  canvNcoll->SaveAs("qaNcoll_NeNe.pdf");
  save::canvasToRoot(canvNcoll, "qaNcoll_NeNe.pdf");
}