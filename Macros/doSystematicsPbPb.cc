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

// Run 2 reference, one entry per syst::CentBins interval (0-5, 5-10, 10-20, ..., 90-100).
// No value is available for 10-20%, so it is set to -1 and not drawn.
static const std::vector<float> NpartValRun1 = { 381.5, 327.8, 259.3, 186.5, 130.1, 86.69, 54.28, 31.04, 15.81, 7.011, 2.788 };
static const std::vector<float> NpartValRun2 = { 383.4, 331.2, 262.0, 187.9, 130.8, 87.14, 54.34, 30.97, 15.72, 6.973, 2.785 };

static const std::vector<float> NcollValRun1 = { 1619.0, 1269.0, 897.7, 553.7, 325.0, 176.6, 88.21, 39.74, 15.96, 5.657, 1.709 };
static const std::vector<float> NcollValRun2 = { 1763.0, 1382.0, 973.4, 592.7, 343.8, 185.7, 91.41, 40.50, 16.12, 5.667, 1.708 };

static const std::vector<float> NpartErrRun1 = { 1.700, 2.80, 2.70, 2.40, 2.20, 2.000, 1.900, 1.30, 0.780, 0.3300, 0.1200 };
static const std::vector<float> NpartErrRun2 = { 0.568, 1.03, 1.15, 1.34, 1.33, 0.928, 0.802, 0.57, 0.241, 0.0729, 0.0497 };

static const std::vector<float> NcollErrRun1 = { 31.0, 27.0, 21.0, 14.0, 9.70, 6.00, 4.00, 2.40, 0.930, 0.31, 0.0990 };
static const std::vector<float> NcollErrRun2 = { 19.4, 15.7, 11.3, 8.21, 5.76, 3.33, 2.11, 1.03, 0.341, 0.10, 0.0474 };

void doSystematicsPbPb()
{
  gStyle->SetOptStat(0);
  syst::GlauberParameters base("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.root");

  syst::GlauberParameters pbHN("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_PbHN_Anchor90.root");
  syst::GlauberParameters pbpn("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbpn_Anchor90.root");
  syst::GlauberParameters pbpnrw("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbpnrw_Anchor90.root");
  syst::GlauberParameters pbrw("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbrw_Anchor90.root");
  syst::GlauberParameters pbstar("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbstar_Anchor90.root");

  syst::GlauberParameters<TH2D> trento("../Trento/trentofit_calib.root", true);

  syst::GlauberParameters anchor89("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor89.root");
  syst::GlauberParameters anchor89p5("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor89.5.root");
  syst::GlauberParameters anchor90p5("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.5.root");
  syst::GlauberParameters anchor91("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor91.root");

  syst::GlauberParameters nancestor0("LHC23_pass5/AR_544122_calibration_ancestorMode0_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root");
  syst::GlauberParameters nancestor1("LHC23_pass5/AR_544122_calibration_ancestorMode1_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root");
  syst::GlauberParameters globalNorm("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_var1_Anchor90.0.root");

  syst::GlauberParameters collisions("LHC23_pass5_small/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_Collisions_basehistos_Pb_Anchor90.0.root");
  syst::GlauberParameters runDeviation("LHC25_pass1/AR_567905_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root");

  ColorManager cm(9);
  std::vector<syst::SystSource> sources = {
    syst::makeSystSource("Nuclei shape", cm.getColor(0),    /*draw*/ false,    /*includeInTotal*/ true,       base, { pbHN, pbpn, pbpnrw, pbrw, pbstar }),
    syst::makeSystSource("Trento", cm.getColor(1),          /*draw*/ false,    /*includeInTotal*/ true,       base, { trento }),
    syst::makeSystSource("0.5% shift", cm.getColor(2),      /*draw*/ false,    /*includeInTotal*/ false,      base, { anchor89p5, anchor90p5 }),
    syst::makeSystSource("1.0% AP shift", cm.getColor(3),   /*draw*/ false,    /*includeInTotal*/ true,       base, { anchor89, anchor91 }),
    syst::makeSystSource("Variant1", cm.getColor(4),        /*draw*/ false,    /*includeInTotal*/ false,      base, { globalNorm }),
    syst::makeSystSource("Variant2", cm.getColor(5),        /*draw*/ false,    /*includeInTotal*/ false,      base, { nancestor0 }),
    syst::makeSystSource("Nanc Rounded", cm.getColor(5),    /*draw*/ false,    /*includeInTotal*/ false,      base, { nancestor1 }),
    syst::makeSystSource("Collisions", cm.getColor(7),      /*draw*/ false,    /*includeInTotal*/ true,       base, { collisions }),
    syst::makeSystSource("Run deviation", cm.getColor(8),   /*draw*/ false,    /*includeInTotal*/ true,       base, { runDeviation }),
  };

  syst::SystSource run1 = syst::makeReferenceSource("Run 1", kBlue, NpartValRun1, NpartErrRun1, NcollValRun1, NcollErrRun1);
  run1.marker = kFullCircle;
  syst::SystSource run2 = syst::makeReferenceSource("Run 2", kRed, NpartValRun2, NpartErrRun2, NcollValRun2, NcollErrRun2);
  run2.marker = kFullSquare;

  sources.push_back(run1);
  sources.push_back(run2);

  syst::SystSource total = syst::combineSourcesInTotal("Total with TRENTo", kBlack, /*draw*/ true, sources);
  total.print(cent::kNpart);

  // syst::drawSystematics(sources, total, cent::kNpart, /*relative*/ false, "<N_{part}> systematic uncertainty", "hNpartSystematics_PbPb.pdf", 2e-6, 8e-1);
  // syst::drawSystematics(sources, total, cent::kNcoll, /*relative*/ false, "<N_{coll}> systematic uncertainty", "hNcollSystematics_PbPb.pdf", 2e-6, 8e-1);
  syst::drawSystematics(sources, total, cent::kNpart, /*relative*/ true, "<N_{part}> relative systematic uncertainty", "hNpartRelSystematics_PbPb.pdf", 5e-4, 3e+0);
  syst::drawSystematics(sources, total, cent::kNcoll, /*relative*/ true, "<N_{coll}> relative systematic uncertainty", "hNcollRelSystematics_PbPb.pdf", 5e-3, 3e+1);
}
