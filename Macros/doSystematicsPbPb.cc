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

static const std::vector<float> NpartVal = { 382.6, 330.6, 265.0, 194.5, 138.8, 94.7, 60.1, 34.5, 17.4, 7.8, 3.4 };
static const std::vector<float> NcollVal = { 1780.9, 1387.0, 988.2, 618.0, 369.9, 206.1, 103.7, 46.1, 18.1, 6.5, 2.3 };

void doSystematicsPbPb()
{
  gStyle->SetOptStat(0);
  syst::GlauberParameters base("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.root");

  syst::GlauberParameters anchor89("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor89.root");
  syst::GlauberParameters anchor89p5("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor89.5.root");
  syst::GlauberParameters anchor90p5("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.5.root");
  syst::GlauberParameters anchor91("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor91.root");

  syst::GlauberParameters pbHN("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_PbHN_Anchor90.root");
  syst::GlauberParameters pbpn("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbpn_Anchor90.root");
  syst::GlauberParameters pbpnrw("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbpnrw_Anchor90.root");
  syst::GlauberParameters pbrw("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbrw_Anchor90.root");
  syst::GlauberParameters pbstar("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbstar_Anchor90.root");

  syst::GlauberParameters<TH2D> trento("../Trento/trentofit_calib.root", true);

  std::vector<std::vector<float>> systShape = syst::computeSystematics(base, {pbHN, pbpn, pbpnrw, pbrw, pbstar});
  std::vector<std::vector<float>> relSystShape = syst::computeRelativeSystematics(base ,{pbHN, pbpn, pbpnrw, pbrw, pbstar});

  std::vector<std::vector<float>> systAnchorShift0p5 = syst::computeSystematics(base, {anchor89p5, anchor90p5});
  std::vector<std::vector<float>> relSystAnchorShift0p5 = syst::computeRelativeSystematics(base ,{anchor89p5, anchor90p5});
  std::vector<std::vector<float>> systAnchorShift1p0 = syst::computeSystematics(base, {anchor89, anchor91});
  std::vector<std::vector<float>> relSystAnchorShift1p0 = syst::computeRelativeSystematics(base ,{anchor89, anchor91});
  
  std::vector<std::vector<float>> systTrento = syst::computeSystematics(base, {trento});
  std::vector<std::vector<float>> relSystTrento = syst::computeRelativeSystematics(base, {trento});

  std::vector<std::vector<float>> systTotal = syst::combineSystematicsInQuadrature(systShape, systTrento, systAnchorShift0p5);
  std::vector<std::vector<float>> relSystTotal = syst::combineSystematicsInQuadrature(relSystShape, relSystTrento, relSystAnchorShift0p5);

  ColorManager cm(6);

  TH1F* hSysTotalNpart = syst::initGlauParHist(systTotal[cent::kNpart], "hSysTotalNpart", kBlack);
  TH1F* hSysTotalNcoll = syst::initGlauParHist(systTotal[cent::kNcoll], "hSysTotalNcoll", kBlack);
  TH1F* hRelSysTotalNpart = syst::initGlauParHist(relSystTotal[cent::kNpart], "hRelSysTotalNpart", kBlack);
  TH1F* hRelSysTotalNcoll = syst::initGlauParHist(relSystTotal[cent::kNcoll], "hRelSysTotalNcoll", kBlack);

  TH1F* hSysShapeNpart = syst::initGlauParHist(systShape[cent::kNpart], "hSysShapeNpart", cm.getColor(0));
  TH1F* hSysShapeNcoll = syst::initGlauParHist(systShape[cent::kNcoll], "hSysShapeNcoll", cm.getColor(0));
  TH1F* hRelSysShapeNpart = syst::initGlauParHist(relSystShape[cent::kNpart], "hRelSysShapeNpart", cm.getColor(0));
  TH1F* hRelSysShapeNcoll = syst::initGlauParHist(relSystShape[cent::kNcoll], "hRelSysShapeNcoll", cm.getColor(0));

  TH1F* hSysTrentoNpart = syst::initGlauParHist(systTrento[cent::kNpart], "hSysTrentoNpart", cm.getColor(1));
  TH1F* hSysTrentoNcoll = syst::initGlauParHist(systTrento[cent::kNcoll], "hSysTrentoNcoll", cm.getColor(1));
  TH1F* hRelSysTrentoNpart = syst::initGlauParHist(relSystTrento[cent::kNpart], "hRelSysTrentoNpart", cm.getColor(1));
  TH1F* hRelSysTrentoNcoll = syst::initGlauParHist(relSystTrento[cent::kNcoll], "hRelSysTrentoNcoll", cm.getColor(1));

  TH1F* hSysAnchorShift0p5Npart = syst::initGlauParHist(systAnchorShift0p5[cent::kNpart], "hSysAnchorShift0p5Npart", cm.getColor(2));
  TH1F* hSysAnchorShift0p5Ncoll = syst::initGlauParHist(systAnchorShift0p5[cent::kNcoll], "hSysAnchorShift0p5Ncoll", cm.getColor(2));
  TH1F* hRelSysAnchorShift0p5Npart = syst::initGlauParHist(relSystAnchorShift0p5[cent::kNpart], "hRelSysAnchorShift0p5Npart", cm.getColor(2));
  TH1F* hRelSysAnchorShift0p5Ncoll = syst::initGlauParHist(relSystAnchorShift0p5[cent::kNcoll], "hRelSysAnchorShift0p5Ncoll", cm.getColor(2));

  TH1F* hSysAnchorShift1p0Npart = syst::initGlauParHist(systAnchorShift1p0[cent::kNpart], "hSysAnchorShift1p0Npart", cm.getColor(3));
  TH1F* hSysAnchorShift1p0Ncoll = syst::initGlauParHist(systAnchorShift1p0[cent::kNcoll], "hSysAnchorShift1p0Ncoll", cm.getColor(3));
  TH1F* hRelSysAnchorShift1p0Npart = syst::initGlauParHist(relSystAnchorShift1p0[cent::kNpart], "hRelSysAnchorShift1p0Npart", cm.getColor(3));
  TH1F* hRelSysAnchorShift1p0Ncoll = syst::initGlauParHist(relSystAnchorShift1p0[cent::kNcoll], "hRelSysAnchorShift1p0Ncoll", cm.getColor(3));


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

  hSysTotalNpartCopy->SetMinimum(1e-2);
  hSysTotalNcollCopy->SetMinimum(1e-2);
  hRelSysTotalNcollCopy->SetMinimum(1e-2);
  hRelSysTotalNpartCopy->SetMinimum(1e-2);

  const int canvWidth = 1600;
  const int canvHeight = 1000;
  const bool drawShape = true;
  const bool drawTrento = true;
  const bool drawSysAnchorShift0p5 = true;
  const bool drawSysAnchorShift1p0 = true;
  const bool drawTotal = false;

  TLegend* leg = new TLegend(0.18, 0.61, 0.43, 0.94);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hSysTotalNpart, "Total", "l");
  leg->AddEntry(hSysShapeNpart, "Nuclei shape", "l");
  leg->AddEntry(hSysTrentoNpart, "Trento", "l");
  leg->AddEntry(hSysAnchorShift0p5Npart, "0.5% shift", "l");
  leg->AddEntry(hSysAnchorShift1p0Npart, "1.0% shift", "l");

  // hRelSysTotalNpartCopy->SetMaximum(1);
  // hRelSysTotalNcollCopy->SetMaximum(1);

  TCanvas* canvNpartSystematics = new TCanvas("canvNpartSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvNpartSystematics);
  syst::styleHist(hSysTotalNpartCopy, "<N_{part}> systematic uncertainty");
  hSysTotalNpartCopy->Draw("");
  if (drawShape) hSysShapeNpart->Draw("hist same");
  if (drawTrento) hSysTrentoNpart->Draw("hist same");
  if (drawSysAnchorShift0p5) hSysAnchorShift0p5Npart->Draw("hist same");
  if (drawSysAnchorShift1p0) hSysAnchorShift1p0Npart->Draw("hist same");
  if (drawTotal) hSysTotalNpart->Draw("hist same");
  leg->Draw();
  canvNpartSystematics->SaveAs("hNpartSystematics_PbPb.pdf");

  TCanvas* canvNcollSystematics = new TCanvas("canvNcollSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvNcollSystematics);
  syst::styleHist(hSysTotalNcollCopy, "<N_{coll}> systematic uncertainty");
  hSysTotalNcollCopy->Draw("");
  if (drawShape) hSysShapeNcoll->Draw("hist same");
  if (drawTrento) hSysTrentoNcoll->Draw("hist same");
  if (drawSysAnchorShift0p5) hSysAnchorShift0p5Ncoll->Draw("hist same");
  if (drawSysAnchorShift1p0) hSysAnchorShift1p0Ncoll->Draw("hist same");
  if (drawTotal) hSysTotalNcoll->Draw("hist same");
  leg->Draw();
  canvNcollSystematics->SaveAs("hNcollSystematics_PbPb.pdf");

  TCanvas* canvRelNpartSystematics = new TCanvas("canvRelNpartSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvRelNpartSystematics);
  syst::styleHist(hRelSysTotalNpartCopy, "<N_{part}> relative systematic uncertainty");
  hRelSysTotalNpartCopy->Draw("");
  if (drawShape) hRelSysShapeNpart->Draw("hist same");
  if (drawTrento) hRelSysTrentoNpart->Draw("hist same");
  if (drawSysAnchorShift0p5) hRelSysAnchorShift0p5Npart->Draw("hist same");
  if (drawSysAnchorShift1p0) hRelSysAnchorShift1p0Npart->Draw("hist same");
  if (drawTotal) hRelSysTotalNpart->Draw("hist same");
  leg->Draw();
  canvRelNpartSystematics->SaveAs("hNpartRelSystematics_PbPb.pdf");

  TCanvas* canvRelNcollSystematics = new TCanvas("canvRelNcollSystematics", "", canvWidth, canvHeight);
  syst::styleCanvas(canvRelNcollSystematics);
  syst::styleHist(hRelSysTotalNcollCopy, "<N_{coll}> relative systematic uncertainty");
  hRelSysTotalNcollCopy->Draw("");
  if (drawShape) hRelSysShapeNcoll->Draw("hist same");
  if (drawTrento) hRelSysTrentoNcoll->Draw("hist same");
  if (drawSysAnchorShift0p5) hRelSysAnchorShift0p5Ncoll->Draw("hist same");
  if (drawSysAnchorShift1p0) hRelSysAnchorShift1p0Ncoll->Draw("hist same");
  if (drawTotal) hRelSysTotalNcoll->Draw("hist same");
  leg->Draw();
  canvRelNcollSystematics->SaveAs("hNcollRelSystematics_PbPb.pdf");
}