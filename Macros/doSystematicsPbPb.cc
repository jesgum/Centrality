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
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  syst::GlauberParameters base("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");
  // std::vector<std::vector<float>> systTotal = syst::combineSystematicsInQuadrature(systAnchor, systBcs, systRun, systFT0C);
  // std::vector<std::vector<float>> relSystTotal = syst::combineSystematicsInQuadrature(relSystAnchor, relSystBcs, relSystRun, relSystFT0C);
}