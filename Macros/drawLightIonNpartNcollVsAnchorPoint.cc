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

void doNpartNCollVsAnchorPoint(TH1F* hBase, std::vector<TH1F*> histos, TLegend* leg, const char* outfile, const char* titleYaxis = "<N_{part}>")
{
  std::vector<std::pair<TH1F*, TH1F*>> runPairs;
  TCanvas* canv = new TCanvas("canv", "", 1200, 800);
  canv->SetTopMargin(0.0);
  canv->SetBottomMargin(0.0);

  TPad* padTop = new TPad("padTop", "", 0, 0.3, 1, 1);
  padTop->SetTicks(1, 1);
  padTop->SetTopMargin(0.04);
  padTop->SetBottomMargin(0.02);
  padTop->SetLeftMargin(0.1);
  padTop->SetRightMargin(0.03);
  padTop->Draw();
  padTop->cd();
  hBase->GetXaxis()->SetLabelSize(0);
  hBase->GetYaxis()->SetLabelSize(0.05);
  hBase->GetYaxis()->SetTitle(titleYaxis);
  hBase->GetYaxis()->SetTitleSize(0.07);
  hBase->GetYaxis()->SetTitleOffset(0.6);
  hBase->Draw("pe");
  for (const auto& hist : histos) {
    hist->Draw("same pe");
    runPairs.push_back({hist, hBase});
  }
  leg->Draw();

  canv->cd();
  TPad* padBot = new TPad("padBot", "", 0, 0, 1, 0.3);
  padBot->SetTicks(1, 1);
  padBot->SetTopMargin(0.02);
  padBot->SetBottomMargin(0.35);
  padBot->SetLeftMargin(0.1);
  padBot->SetRightMargin(0.03);
  padBot->Draw();
  padBot->cd();

  std::vector<TH1F*> hRatios;
  double globalMin = 1e9, globalMax = -1e9;

  for (int i = 0; i < runPairs.size(); ++i) {
    TH1F* hRatio = (TH1F*)runPairs[i].first->Clone(Form("hRatioNpart_%d", i));
    hRatio->Divide(hRatio, runPairs[i].second, 1, 1, "B");
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Ratio");
    hRatio->GetXaxis()->SetTitle("Centrality (%)");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);

    // scan bin content +/- error to get the true drawn extent
    for (int b = 1; b <= hRatio->GetNbinsX(); ++b) {
      double content = hRatio->GetBinContent(b);
      double error   = hRatio->GetBinError(b);
      if (content == 0 && error == 0) continue; // skip empty bins
      double lo = content - error;
      double hi = content + error;
      if (lo < globalMin) globalMin = lo;
      if (hi > globalMax) globalMax = hi;
    }

    hRatios.push_back(hRatio);
  }

  // add a bit of padding so points/error bars aren't flush with the pad edge
  double margin = (globalMax - globalMin) * 0.1;
  if (margin <= 0) margin = 0.01; // fallback in case all values are identical
  globalMin -= margin;
  globalMax += margin;

  for (auto* hRatio : hRatios) {
    hRatio->GetYaxis()->SetRangeUser(globalMin, globalMax);
  }

  hRatios[0]->Draw("pe");
  for (int i = 1; i < runPairs.size(); ++i) hRatios[i]->Draw("pe same");

  canv->SaveAs(outfile);
  for (auto* hRatio : hRatios) delete hRatio;
  hRatios.clear();
  delete canv;
}

void drawLightIonNpartNcollVsAnchorPoint()
{
  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  // syst::GlauberParameters base("LHC25ae_pass2/calibOO.root");
  syst::GlauberParameters base("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters col50("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor50.root");
  syst::GlauberParameters col60("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor60.root");
  syst::GlauberParameters col70("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor70.root");
  syst::GlauberParameters col80("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor80.root");
  syst::GlauberParameters col90("LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_Anchor90.root");
  syst::GlauberParameters bc50("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor50.root");
  syst::GlauberParameters bc60("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor60.root");
  syst::GlauberParameters bc70("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor70.root");
  syst::GlauberParameters bc80("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor80.root");
  syst::GlauberParameters bc90("LHC25ae_pass2_systematics/AR_564374_calibration_hFT0M_BCs_Anchor90.root");

  syst::GlauberParameters baseNeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_LightIonDef.root");
  syst::GlauberParameters col50NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor50.root");
  syst::GlauberParameters col60NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor60.root");
  syst::GlauberParameters col70NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor70.root");
  syst::GlauberParameters col80NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor80.root");
  syst::GlauberParameters col90NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor90.root");
  syst::GlauberParameters col91NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor91.root");
  syst::GlauberParameters col92NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor92.root");
  syst::GlauberParameters col93NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor93.root");
  syst::GlauberParameters col94NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor94.root");
  syst::GlauberParameters col95NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor95.root");
  syst::GlauberParameters col96NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor96.root");
  syst::GlauberParameters col97NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor97.root");
  syst::GlauberParameters col98NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor98.root");
  syst::GlauberParameters col99NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_Collisions_Anchor99.root");

  syst::GlauberParameters bc50NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor50.root");
  syst::GlauberParameters bc60NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor60.root");
  syst::GlauberParameters bc70NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor70.root");
  syst::GlauberParameters bc80NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor80.root");
  syst::GlauberParameters bc90NeNe("LHC25af_pass2_systematics/AR_564468_calibration_hFT0M_BCs_Anchor90.root");

  ColorManager cm3(5);
  TH1F* hBaseNpart = syst::initGlauParHist(base.npart, "hBaseNpart", kBlack);
  TH1F* hBaseNcoll = syst::initGlauParHist(base.ncoll, "hBaseNcoll", kBlack);
  TH1F* hCol50Npart = syst::initGlauParHist(col50.npart, "hCol50Npart", cm3.getColor(0));
  TH1F* hCol60Npart = syst::initGlauParHist(col60.npart, "hCol60Npart", cm3.getColor(1));
  TH1F* hCol70Npart = syst::initGlauParHist(col70.npart, "hCol70Npart", cm3.getColor(2));
  TH1F* hCol80Npart = syst::initGlauParHist(col80.npart, "hCol80Npart", cm3.getColor(3));
  TH1F* hCol90Npart = syst::initGlauParHist(col90.npart, "hCol90Npart", cm3.getColor(4));
  TH1F* hBc50Npart = syst::initGlauParHist(bc50.npart, "hBc50Npart", cm3.getColor(0));
  TH1F* hBc60Npart = syst::initGlauParHist(bc60.npart, "hBc60Npart", cm3.getColor(1));
  TH1F* hBc70Npart = syst::initGlauParHist(bc70.npart, "hBc70Npart", cm3.getColor(2));
  TH1F* hBc80Npart = syst::initGlauParHist(bc80.npart, "hBc80Npart", cm3.getColor(3));
  TH1F* hBc90Npart = syst::initGlauParHist(bc90.npart, "hBc90Npart", cm3.getColor(4));
  TH1F* hCol50Ncoll = syst::initGlauParHist(col50.ncoll, "hCol50Ncoll", cm3.getColor(0));
  TH1F* hCol60Ncoll = syst::initGlauParHist(col60.ncoll, "hCol60Ncoll", cm3.getColor(1));
  TH1F* hCol70Ncoll = syst::initGlauParHist(col70.ncoll, "hCol70Ncoll", cm3.getColor(2));
  TH1F* hCol80Ncoll = syst::initGlauParHist(col80.ncoll, "hCol80Ncoll", cm3.getColor(3));
  TH1F* hCol90Ncoll = syst::initGlauParHist(col90.ncoll, "hCol90Ncoll", cm3.getColor(4));
  TH1F* hBc50Ncoll = syst::initGlauParHist(bc50.ncoll, "hBc50Ncoll", cm3.getColor(0));
  TH1F* hBc60Ncoll = syst::initGlauParHist(bc60.ncoll, "hBc60Ncoll", cm3.getColor(1));
  TH1F* hBc70Ncoll = syst::initGlauParHist(bc70.ncoll, "hBc70Ncoll", cm3.getColor(2));
  TH1F* hBc80Ncoll = syst::initGlauParHist(bc80.ncoll, "hBc80Ncoll", cm3.getColor(3));
  TH1F* hBc90Ncoll = syst::initGlauParHist(bc90.ncoll, "hBc90Ncoll", cm3.getColor(4));
  TH1F* hBaseNpartNeNe = syst::initGlauParHist(baseNeNe.npart, "hBaseNpart", kBlack);
  TH1F* hBaseNcollNeNe = syst::initGlauParHist(baseNeNe.ncoll, "hBaseNcoll", kBlack);
  TH1F* hCol50NpartNeNe = syst::initGlauParHist(col50NeNe.npart, "hCol50Npart", cm3.getColor(0));
  TH1F* hCol60NpartNeNe = syst::initGlauParHist(col60NeNe.npart, "hCol60Npart", cm3.getColor(1));
  TH1F* hCol70NpartNeNe = syst::initGlauParHist(col70NeNe.npart, "hCol70Npart", cm3.getColor(2));
  TH1F* hCol80NpartNeNe = syst::initGlauParHist(col80NeNe.npart, "hCol80Npart", cm3.getColor(3));
  TH1F* hCol90NpartNeNe = syst::initGlauParHist(col90NeNe.npart, "hCol90Npart", cm3.getColor(4));
  TH1F* hBc50NpartNeNe = syst::initGlauParHist(bc50NeNe.npart, "hBc50Npart", cm3.getColor(0));
  TH1F* hBc60NpartNeNe = syst::initGlauParHist(bc60NeNe.npart, "hBc60Npart", cm3.getColor(1));
  TH1F* hBc70NpartNeNe = syst::initGlauParHist(bc70NeNe.npart, "hBc70Npart", cm3.getColor(2));
  TH1F* hBc80NpartNeNe = syst::initGlauParHist(bc80NeNe.npart, "hBc80Npart", cm3.getColor(3));
  TH1F* hBc90NpartNeNe = syst::initGlauParHist(bc90NeNe.npart, "hBc90Npart", cm3.getColor(4));
  TH1F* hCol50NcollNeNe = syst::initGlauParHist(col50NeNe.ncoll, "hCol50Ncoll", cm3.getColor(0));
  TH1F* hCol60NcollNeNe = syst::initGlauParHist(col60NeNe.ncoll, "hCol60Ncoll", cm3.getColor(1));
  TH1F* hCol70NcollNeNe = syst::initGlauParHist(col70NeNe.ncoll, "hCol70Ncoll", cm3.getColor(2));
  TH1F* hCol80NcollNeNe = syst::initGlauParHist(col80NeNe.ncoll, "hCol80Ncoll", cm3.getColor(3));
  TH1F* hCol90NcollNeNe = syst::initGlauParHist(col90NeNe.ncoll, "hCol90Ncoll", cm3.getColor(4));
  TH1F* hBc50NcollNeNe = syst::initGlauParHist(bc50NeNe.ncoll, "hBc50Ncoll", cm3.getColor(0));
  TH1F* hBc60NcollNeNe = syst::initGlauParHist(bc60NeNe.ncoll, "hBc60Ncoll", cm3.getColor(1));
  TH1F* hBc70NcollNeNe = syst::initGlauParHist(bc70NeNe.ncoll, "hBc70Ncoll", cm3.getColor(2));
  TH1F* hBc80NcollNeNe = syst::initGlauParHist(bc80NeNe.ncoll, "hBc80Ncoll", cm3.getColor(3));
  TH1F* hBc90NcollNeNe = syst::initGlauParHist(bc90NeNe.ncoll, "hBc90Ncoll", cm3.getColor(4));

  std::vector<TH1F*> hColNpart = { hCol50Npart, hCol60Npart, hCol70Npart, hCol80Npart, hCol90Npart };
  std::vector<TH1F*> hBcNpart = { hBc50Npart, hBc60Npart, hBc70Npart, hBc80Npart, hBc90Npart };
  std::vector<TH1F*> hColNcoll = { hCol50Ncoll, hCol60Ncoll, hCol70Ncoll, hCol80Ncoll, hCol90Ncoll };
  std::vector<TH1F*> hBcNcoll = { hBc50Ncoll, hBc60Ncoll, hBc70Ncoll, hBc80Ncoll, hBc90Ncoll };
  std::vector<TH1F*> hColNpartNeNe = { hCol50NpartNeNe, hCol60NpartNeNe, hCol70NpartNeNe, hCol80NpartNeNe, hCol90NpartNeNe };
  std::vector<TH1F*> hBcNpartNeNe = { hBc50NpartNeNe, hBc60NpartNeNe, hBc70NpartNeNe, hBc80NpartNeNe, hBc90NpartNeNe };
  std::vector<TH1F*> hColNcollNeNe = { hCol50NcollNeNe, hCol60NcollNeNe, hCol70NcollNeNe, hCol80NcollNeNe, hCol90NcollNeNe };
  std::vector<TH1F*> hBcNcollNeNe = { hBc50NcollNeNe, hBc60NcollNeNe, hBc70NcollNeNe, hBc80NcollNeNe, hBc90NcollNeNe };

  TLegend* leg = new TLegend(0.75, 0.35, 0.95, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);
  leg->AddEntry(hBaseNpart, "Reference", "pl");
  leg->AddEntry(hCol50Npart, "Anchor 50", "pl");
  leg->AddEntry(hCol60Npart, "Anchor 60", "pl");
  leg->AddEntry(hCol70Npart, "Anchor 70", "pl");
  leg->AddEntry(hCol80Npart, "Anchor 80", "pl");
  leg->AddEntry(hCol90Npart, "Anchor 90", "pl");

  doNpartNCollVsAnchorPoint(hBaseNpart, hColNpart, leg, "FiguresSystematics/hColNpart.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hBaseNpart, hBcNpart, leg, "FiguresSystematics/hBcNpart.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hBaseNcoll, hColNcoll, leg, "FiguresSystematics/hColNcoll.pdf", "<N_{coll}>");
  doNpartNCollVsAnchorPoint(hBaseNcoll, hBcNcoll, leg, "FiguresSystematics/hBcNcoll.pdf", "<N_{coll}>");
  doNpartNCollVsAnchorPoint(hBaseNpartNeNe, hColNpartNeNe, leg, "FiguresSystematics/hColNpartNeNe.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hBaseNpartNeNe, hBcNpartNeNe, leg, "FiguresSystematics/hBcNpartNeNe.pdf", "<N_{part}>");
  doNpartNCollVsAnchorPoint(hBaseNcollNeNe, hColNcollNeNe, leg, "FiguresSystematics/hColNcollNeNe.pdf", "<N_{coll}>");
  doNpartNCollVsAnchorPoint(hBaseNcollNeNe, hBcNcollNeNe, leg, "FiguresSystematics/hBcNcollNeNe.pdf", "<N_{coll}>");

  ColorManager cm4(14);
  TH1F* hCol50NpartNeNeExtra = syst::initGlauParHist(col50NeNe.npart, "hCol50NpartNeNeExtra", cm4.getColor(0));
  TH1F* hCol60NpartNeNeExtra = syst::initGlauParHist(col60NeNe.npart, "hCol60NpartNeNeExtra", cm4.getColor(1));
  TH1F* hCol70NpartNeNeExtra = syst::initGlauParHist(col70NeNe.npart, "hCol70NpartNeNeExtra", cm4.getColor(2));
  TH1F* hCol80NpartNeNeExtra = syst::initGlauParHist(col80NeNe.npart, "hCol80NpartNeNeExtra", cm4.getColor(3));
  TH1F* hCol90NpartNeNeExtra = syst::initGlauParHist(col90NeNe.npart, "hCol90NpartNeNeExtra", cm4.getColor(4));
  TH1F* hCol91NpartNeNeExtra = syst::initGlauParHist(col91NeNe.npart, "hCol91NpartNeNeExtra", cm4.getColor(5));
  TH1F* hCol92NpartNeNeExtra = syst::initGlauParHist(col92NeNe.npart, "hCol92NpartNeNeExtra", cm4.getColor(6));
  TH1F* hCol93NpartNeNeExtra = syst::initGlauParHist(col93NeNe.npart, "hCol93NpartNeNeExtra", cm4.getColor(7));
  TH1F* hCol94NpartNeNeExtra = syst::initGlauParHist(col94NeNe.npart, "hCol94NpartNeNeExtra", cm4.getColor(8));
  TH1F* hCol95NpartNeNeExtra = syst::initGlauParHist(col95NeNe.npart, "hCol95NpartNeNeExtra", cm4.getColor(9));
  TH1F* hCol96NpartNeNeExtra = syst::initGlauParHist(col96NeNe.npart, "hCol96NpartNeNeExtra", cm4.getColor(10));
  TH1F* hCol97NpartNeNeExtra = syst::initGlauParHist(col97NeNe.npart, "hCol97NpartNeNeExtra", cm4.getColor(11));
  TH1F* hCol98NpartNeNeExtra = syst::initGlauParHist(col98NeNe.npart, "hCol98NpartNeNeExtra", cm4.getColor(12));
  TH1F* hCol99NpartNeNeExtra = syst::initGlauParHist(col99NeNe.npart, "hCol99NpartNeNeExtra", cm4.getColor(13));

  TH1F* hCol50NcollNeNeExtra = syst::initGlauParHist(col50NeNe.ncoll, "hCol50NcollNeNeExtra", cm4.getColor(0));
  TH1F* hCol60NcollNeNeExtra = syst::initGlauParHist(col60NeNe.ncoll, "hCol60NcollNeNeExtra", cm4.getColor(1));
  TH1F* hCol70NcollNeNeExtra = syst::initGlauParHist(col70NeNe.ncoll, "hCol70NcollNeNeExtra", cm4.getColor(2));
  TH1F* hCol80NcollNeNeExtra = syst::initGlauParHist(col80NeNe.ncoll, "hCol80NcollNeNeExtra", cm4.getColor(3));
  TH1F* hCol90NcollNeNeExtra = syst::initGlauParHist(col90NeNe.ncoll, "hCol90NcollNeNeExtra", cm4.getColor(4));
  TH1F* hCol91NcollNeNeExtra = syst::initGlauParHist(col91NeNe.ncoll, "hCol91NcollNeNeExtra", cm4.getColor(5));
  TH1F* hCol92NcollNeNeExtra = syst::initGlauParHist(col92NeNe.ncoll, "hCol92NcollNeNeExtra", cm4.getColor(6));
  TH1F* hCol93NcollNeNeExtra = syst::initGlauParHist(col93NeNe.ncoll, "hCol93NcollNeNeExtra", cm4.getColor(7));
  TH1F* hCol94NcollNeNeExtra = syst::initGlauParHist(col94NeNe.ncoll, "hCol94NcollNeNeExtra", cm4.getColor(8));
  TH1F* hCol95NcollNeNeExtra = syst::initGlauParHist(col95NeNe.ncoll, "hCol95NcollNeNeExtra", cm4.getColor(9));
  TH1F* hCol96NcollNeNeExtra = syst::initGlauParHist(col96NeNe.ncoll, "hCol96NcollNeNeExtra", cm4.getColor(10));
  TH1F* hCol97NcollNeNeExtra = syst::initGlauParHist(col97NeNe.ncoll, "hCol97NcollNeNeExtra", cm4.getColor(11));
  TH1F* hCol98NcollNeNeExtra = syst::initGlauParHist(col98NeNe.ncoll, "hCol98NcollNeNeExtra", cm4.getColor(12));
  TH1F* hCol99NcollNeNeExtra = syst::initGlauParHist(col99NeNe.ncoll, "hCol99NcollNeNeExtra", cm4.getColor(13));


  TLegend* leg2 = new TLegend(0.75, 0.35, 0.95, 0.9);
  leg2->SetBorderSize(0);
  leg2->SetFillColorAlpha(0, 0);
  leg2->AddEntry(hBaseNpart, "Reference", "pl");
  leg2->AddEntry(hCol50NpartNeNeExtra, "Anchor 50", "pl");
  leg2->AddEntry(hCol60NpartNeNeExtra, "Anchor 60", "pl");
  leg2->AddEntry(hCol70NpartNeNeExtra, "Anchor 70", "pl");
  leg2->AddEntry(hCol80NpartNeNeExtra, "Anchor 80", "pl");
  leg2->AddEntry(hCol90NpartNeNeExtra, "Anchor 90", "pl");
  leg2->AddEntry(hCol91NpartNeNeExtra, "Anchor 91", "pl");
  leg2->AddEntry(hCol92NpartNeNeExtra, "Anchor 92", "pl");
  leg2->AddEntry(hCol93NpartNeNeExtra, "Anchor 93", "pl");
  leg2->AddEntry(hCol94NpartNeNeExtra, "Anchor 94", "pl");
  leg2->AddEntry(hCol95NpartNeNeExtra, "Anchor 95", "pl");
  leg2->AddEntry(hCol96NpartNeNeExtra, "Anchor 96", "pl");
  leg2->AddEntry(hCol97NpartNeNeExtra, "Anchor 97", "pl");
  leg2->AddEntry(hCol98NpartNeNeExtra, "Anchor 98", "pl");
  leg2->AddEntry(hCol99NpartNeNeExtra, "Anchor 99", "pl");

  std::vector<TH1F*> hColNpartNeNeExtra = { hCol50NpartNeNeExtra, hCol60NpartNeNeExtra, hCol70NpartNeNeExtra, hCol80NpartNeNeExtra, hCol90NpartNeNeExtra, hCol91NpartNeNeExtra, hCol92NpartNeNeExtra, hCol93NpartNeNeExtra, hCol94NpartNeNeExtra, hCol95NpartNeNeExtra, hCol96NpartNeNeExtra, hCol97NpartNeNeExtra, hCol98NpartNeNeExtra, hCol99NpartNeNeExtra };
  std::vector<TH1F*> hColNcollNeNeExtra = { hCol50NcollNeNeExtra, hCol60NcollNeNeExtra, hCol70NcollNeNeExtra, hCol80NcollNeNeExtra, hCol90NcollNeNeExtra, hCol91NcollNeNeExtra, hCol92NcollNeNeExtra, hCol93NcollNeNeExtra, hCol94NcollNeNeExtra, hCol95NcollNeNeExtra, hCol96NcollNeNeExtra, hCol97NcollNeNeExtra, hCol98NcollNeNeExtra, hCol99NcollNeNeExtra };

  doNpartNCollVsAnchorPoint(hBaseNpartNeNe, hColNpartNeNeExtra, leg2, "hColNcollNeNeExtra.pdf", "<N_{coll}>");
  doNpartNCollVsAnchorPoint(hBaseNcollNeNe, hColNcollNeNeExtra, leg2, "hcolNpartNeNeExtra.pdf", "<N_{part}>");
}