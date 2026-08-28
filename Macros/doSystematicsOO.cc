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

// OO Npart              OO Ncoll      
// 26.34     0.59        42.31     2.61
// 24.08     0.69        35.82     2.70
// 21.03     0.84        28.77     2.73
// 16.88     0.89        20.58     2.35
// 13.05     0.75        14.19     1.75
// 9.88      0.61        9.64      1.22
// 7.39      0.48        6.51      0.82
// 5.51      0.36        4.41      0.54
// 4.15      0.26        3.00      0.35
// 3.17      0.23        2.07      0.26
// 2.42      0.16        1.38      0.16
static const std::vector<float> NpartVal = { 26.34, 24.08, 21.03, 16.88, 13.05, 9.88, 7.39, 5.51, 4.15, 3.17, 2.42 };
static const std::vector<float> NpartErr = { 0.59, 0.69, 0.84, 0.89, 0.75, 0.61, 0.48, 0.36, 0.26, 0.23, 0.16 };
static const std::vector<float> NpartRelErr = { 0.022399, 0.028654, 0.039943, 0.052720, 0.057463, 0.061743, 0.064954, 0.065336, 0.062651, 0.072555, 0.066116};
static const std::vector<float> NcollVal = { 42.31, 35.82, 28.77, 20.58, 14.19, 9.64, 6.51, 4.41, 3.00, 2.07, 1.38 };
static const std::vector<float> NcollErr = { 2.61, 2.70, 2.73, 2.35, 1.75, 1.22, 0.82, 0.54, 0.35, 0.26, 0.16 };
static const std::vector<float> NcollRelErr = { 0.061688, 0.075377, 0.094876, 0.114188, 0.123326, 0.126535, 0.125960, 0.122473, 0.116667, 0.125604, 0.115942};
// static const std::vector<float> NpartVal = { 10.87 };
// static const std::vector<float> NpartErr = { 0.62 };
// static const std::vector<float> NpartRelErr = { NpartErr[0] / NpartVal[0] };
// static const std::vector<float> NcollVal = { 12.96 };
// static const std::vector<float> NcollErr = { 1.47 };
// static const std::vector<float> NcollRelErr = { NcollErr[0] / NcollVal[0] };

void doSystematicsOO()
{
  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  syst::GlauberParameters base("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");
  std::vector<std::vector<float>> systFT0C = { NpartErr, NcollErr };
  std::vector<std::vector<float>> relSystFT0C = { NpartRelErr, NcollRelErr };

  syst::GlauberParameters bcAnchor50("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor50.root");
  syst::GlauberParameters bcAnchor60("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor60.root");
  syst::GlauberParameters bcAnchor70("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor70.root");
  syst::GlauberParameters bcAnchor80("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor80.root");
  syst::GlauberParameters bcAnchor90("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor90.root");
  // syst::GlauberParameters bcAnchor70("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_LightIonDef.root");
  std::vector<std::vector<float>> systBcs = syst::computeSystematics(base, { bcAnchor70 });
  std::vector<std::vector<float>> relSystBcs = syst::computeRelativeSystematics(base, { bcAnchor70 });
  // std::vector<std::vector<float>> systBcs = syst::computeSystematics(base, { bcAnchor50, bcAnchor60, bcAnchor70, bcAnchor80, bcAnchor90 });
  // std::vector<std::vector<float>> relSystBcs = syst::computeRelativeSystematics(base, { bcAnchor50, bcAnchor60, bcAnchor70, bcAnchor80, bcAnchor90 });

  syst::GlauberParameters colAnchor50("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor50.root");
  syst::GlauberParameters colAnchor60("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor60.root");
  syst::GlauberParameters colAnchor70("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor70.root");
  syst::GlauberParameters colAnchor80("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor80.root");
  syst::GlauberParameters colAnchor90("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor90.root");
  std::vector<std::vector<float>> systAnchor = syst::computeSystematics(base, { colAnchor70 });
  std::vector<std::vector<float>> relSystAnchor = syst::computeRelativeSystematics(base, { colAnchor70 });
  // std::vector<std::vector<float>> systAnchor = syst::computeSystematics(base, { colAnchor50, colAnchor60, colAnchor70, colAnchor80, colAnchor90 });
  // std::vector<std::vector<float>> relSystAnchor = syst::computeRelativeSystematics(base, { colAnchor50, colAnchor60, colAnchor70, colAnchor80, colAnchor90 });

  syst::GlauberParameters run564356("LHC25ae_pass2/AR_564356_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564359("LHC25ae_pass2/AR_564359_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564373("LHC25ae_pass2/AR_564373_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564387("LHC25ae_pass2/AR_564387_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564400("LHC25ae_pass2/AR_564400_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564414("LHC25ae_pass2/AR_564414_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564430("LHC25ae_pass2/AR_564430_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters run564445("LHC25ae_pass2/AR_564445_calibration_hFT0M_Collisions_LightIonDef.root");
  std::vector<std::vector<float>> systRun = syst::computeSystematics(base, { run564356 });
  std::vector<std::vector<float>> relSystRun = syst::computeRelativeSystematics(base, { run564356 });
  // std::vector<std::vector<float>> systRun = syst::computeSystematics(base, { run564356, run564359, run564373, run564387, run564400, run564414, run564430, run564445 });
  // std::vector<std::vector<float>> relSystRun = syst::computeRelativeSystematics(base, { run564356, run564359, run564373, run564387, run564400, run564414, run564430, run564445 });

  std::vector<std::vector<float>> systTotal = syst::combineSystematicsInQuadrature(systAnchor, systBcs, systRun, systFT0C);
  std::vector<std::vector<float>> relSystTotal = syst::combineSystematicsInQuadrature(relSystAnchor, relSystBcs, relSystRun, relSystFT0C);

  ColorManager cm1(4);
  TH1F* hSysShapeNpart = syst::initGlauParHist(systFT0C[kNpart], "hSysShapeNpart", cm1.getColor(0));
  TH1F* hSysShapeNcoll = syst::initGlauParHist(systFT0C[kNcoll], "hSysShapeNcoll", cm1.getColor(0));
  TH1F* hRelSysShapeNpart = syst::initGlauParHist(relSystFT0C[kNpart], "hRelSysShapeNpart", cm1.getColor(0));
  TH1F* hRelSysShapeNcoll = syst::initGlauParHist(relSystFT0C[kNcoll], "hRelSysShapeNcoll", cm1.getColor(0));

  TH1F* hSysAnchorNpart = syst::initGlauParHist(systAnchor[kNpart], "hSysAnchorNpart", cm1.getColor(1));
  TH1F* hSysAnchorNcoll = syst::initGlauParHist(systAnchor[kNcoll], "hSysAnchorNcoll", cm1.getColor(1));
  TH1F* hRelSysAnchorNpart = syst::initGlauParHist(relSystAnchor[kNpart], "hRelSysAnchorNpart", cm1.getColor(1));
  TH1F* hRelSysAnchorNcoll = syst::initGlauParHist(relSystAnchor[kNcoll], "hRelSysAnchorNcoll", cm1.getColor(1));

  TH1F* hSysBcsNpart = syst::initGlauParHist(systBcs[kNpart], "hSysBcsNpart", cm1.getColor(2));
  TH1F* hSysBcsNcoll = syst::initGlauParHist(systBcs[kNcoll], "hSysBcsNcoll", cm1.getColor(2));
  TH1F* hRelSysBcsNpart = syst::initGlauParHist(relSystBcs[kNpart], "hRelSysBcsNpart", cm1.getColor(2));
  TH1F* hRelSysBcsNcoll = syst::initGlauParHist(relSystBcs[kNcoll], "hRelSysBcsNcoll", cm1.getColor(2));

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
  const bool drawBcs = true;
  const bool drawRun = true;
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
  hRelSysTotalNcollCopy->SetMaximum(0.2);
  hRelSysTotalNpartCopy->SetMinimum(0);
  hRelSysTotalNcollCopy->SetMinimum(0);

  TCanvas* canvNpartSystematics = new TCanvas("canvNpartSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvNpartSystematics);
  syst::styleHist(hSysTotalNpartCopy, "<N_{part}> systematic uncertainty");
  hSysTotalNpartCopy->Draw("");
  if (drawShape) hSysShapeNpart->Draw("hist same");
  if (drawAnchor) hSysAnchorNpart->Draw("hist same");
  if (drawRun) hSysRunNpart->Draw("hist same");
  if (drawBcs) hSysBcsNpart->Draw("hist same");
  if (drawTotal) hSysTotalNpart->Draw("hist same");
  leg->Draw();
  canvNpartSystematics->SaveAs("hNpartSystematics_OO.pdf");

  TCanvas* canvNcollSystematics = new TCanvas("canvNcollSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvNcollSystematics);
  syst::styleHist(hSysTotalNcollCopy, "<N_{coll}> systematic uncertainty");
  hSysTotalNcollCopy->Draw("");
  if (drawShape) hSysShapeNcoll->Draw("hist same");
  if (drawAnchor) hSysAnchorNcoll->Draw("hist same");
  if (drawRun) hSysRunNcoll->Draw("hist same");
  if (drawBcs) hSysBcsNcoll->Draw("hist same");
  if (drawTotal) hSysTotalNcoll->Draw("hist same");
  leg->Draw();
  canvNcollSystematics->SaveAs("hNcollSystematics_OO.pdf");

  TCanvas* canvRelNpartSystematics = new TCanvas("canvRelNpartSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvRelNpartSystematics);
  syst::styleHist(hRelSysTotalNpartCopy, "<N_{part}> relative systematic uncertainty");
  hRelSysTotalNpartCopy->Draw("");
  if (drawShape) hRelSysShapeNpart->Draw("hist same");
  if (drawAnchor) hRelSysAnchorNpart->Draw("hist same");
  if (drawRun) hRelSysRunNpart->Draw("hist same");
  if (drawBcs) hRelSysBcsNpart->Draw("hist same");
  if (drawTotal) hRelSysTotalNpart->Draw("hist same");
  leg->Draw();
  canvRelNpartSystematics->SaveAs("FiguresSystematics/hNpartRelSystematics_OO.pdf");

  TCanvas* canvRelNcollSystematics = new TCanvas("canvRelNcollSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvRelNcollSystematics);
  syst::styleHist(hRelSysTotalNcollCopy, "<N_{coll}> relative systematic uncertainty");
  hRelSysTotalNcollCopy->Draw("");
  if (drawShape) hRelSysShapeNcoll->Draw("hist same");
  if (drawAnchor) hRelSysAnchorNcoll->Draw("hist same");
  if (drawRun) hRelSysRunNcoll->Draw("hist same");
  if (drawBcs) hRelSysBcsNcoll->Draw("hist same");
  if (drawTotal) hRelSysTotalNcoll->Draw("hist same");
  leg->Draw();
  canvRelNcollSystematics->SaveAs("FiguresSystematics/hNcollRelSystematics_OO.pdf");


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
    const float errNpartFT0M = systTotal[kNpart][ii];
    const float errNcollFT0M = systTotal[kNcoll][ii];

    errNpart.push_back(errNpartFT0C);
    errNcoll.push_back(errNcollFT0C);
    if (publishValues) {
      std::cout << Form("%.f-%.f%% | Npart: %.2f +- %.2f | Ncoll: %.2f +- %.2f", syst::CentBins[ii], syst::CentBins[ii + 1], NpartVal[ii], errNpartFT0C, NcollVal[ii], errNcollFT0C) << std::endl;
    }
    if (publishCrossCheck) {
      std::cout << Form("%.f-%.f%% | Npart: %.4f +- %.4f | Ncoll: %.4f +- %.4f", syst::CentBins[ii], syst::CentBins[ii + 1], errNpartFT0M, systTotal[kNpart][ii], errNcollFT0M, systTotal[kNcoll][ii]) << std::endl;
    }
  }

  TH1F* hOldNpartErr = syst::initGlauParHist(NpartErr, "hOldNpartErr", kBlack);
  TH1F* hNewNpartErr = syst::initGlauParHist(errNpart, "hNewNpartErr", kBlue);
  TH1F* hOldNcollErr = syst::initGlauParHist(NcollErr, "hOldNcollErr", kBlack);
  TH1F* hNewNcollErr = syst::initGlauParHist(errNcoll, "hNewNcollErr", kBlue);

  TLegend* legSystematics = new TLegend(0.66, 0.7, 0.96, 0.9);
  legSystematics->SetBorderSize(0);
  legSystematics->SetFillColorAlpha(0, 0);
  legSystematics->AddEntry(hOldNpartErr, "Previous", "l");
  legSystematics->AddEntry(hNewNpartErr, "New", "l");

  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 1000);
  syst::styleCanvas(canvNpart);
  canvNpart->SetRightMargin(0.03);
  syst::styleHist(hOldNpartErr, "Uncertainty");
  hOldNpartErr->SetMaximum(1.2);
  hOldNpartErr->SetMinimum(0);
  hOldNpartErr->Draw("hist");
  hNewNpartErr->Draw("hist same");
  legSystematics->Draw();
  canvNpart->SaveAs("qaNpart_OO.pdf");
  
  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 1000);
  syst::styleCanvas(canvNcoll);
  canvNcoll->SetRightMargin(0.03);
  syst::styleHist(hOldNcollErr, "Uncertainty");
  hOldNcollErr->SetMaximum(3.3);
  hOldNcollErr->SetMinimum(0);
  hOldNcollErr->Draw("hist");
  hNewNcollErr->Draw("hist same");
  legSystematics->Draw();
  canvNcoll->SaveAs("qaNcoll_OO.pdf");
}