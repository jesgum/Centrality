#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"
#include "Includes/saveHelper.h"

void drawMarkerForPercentile(TCanvas* canv, TH2* h2d, const char* tag)
{
  canv->cd();
  const TAxis* ax = h2d->GetXaxis();
  TH1D* hAvg = new TH1D(Form("hAvg_%s_%s", h2d->GetName(), tag), "", ax->GetNbins(), ax->GetXmin(), ax->GetXmax());
  for (int ii = 0; ii < h2d->GetNbinsX(); ++ii) {
    TH1D* hProj = h2d->ProjectionY(Form("hProj_%s_%d", tag, ii + 1), ii + 1, ii + 1);
    hAvg->SetBinContent(ii + 1, hProj->GetMean());
    hAvg->SetBinError(ii + 1, hProj->GetMeanError());
    delete hProj;
  }
  hAvg->SetMarkerColor(kBlack);
  hAvg->SetLineColor(kBlack);
  hAvg->SetMarkerStyle(kFullCircle);
  hAvg->Draw("same");
}

// Draw the Npart and Ncoll vs centrality 2D histograms of one calibration file.
// tag goes into the canvas and output names (e.g. "OO" -> h2dNpartOO.pdf, canvas canvNpartOO)
void drawSystem(const char* fileName, const char* tag, const char* estimator,
                double npartMax, double ncollMax, int ncollRebin, const char* canvSuffix)
{
  TFile* infile = new TFile(fileName, "read");
  if (!infile || infile->IsZombie()) {
    std::cerr << "ERROR: cannot open " << fileName << std::endl;
    return;
  }
  TH2* h2dNpart = (TH2*)infile->Get("h2dNpart");
  TH2* h2dNcoll = (TH2*)infile->Get("h2dNcoll");
  if (!h2dNpart || !h2dNcoll) {
    std::cerr << "ERROR: h2dNpart / h2dNcoll missing in " << fileName << std::endl;
    return;
  }

  if (ncollRebin > 1) {
    h2dNcoll->RebinY(ncollRebin);
  }
  h2dNpart->SetMinimum(1e-10);
  h2dNcoll->SetMinimum(1e-10);
  h2dNpart->GetYaxis()->SetRangeUser(0, npartMax);
  h2dNcoll->GetYaxis()->SetRangeUser(0, ncollMax);
  h2dNpart->GetYaxis()->SetTitle("#it{N}_{part}");
  h2dNcoll->GetYaxis()->SetTitle("#it{N}_{coll}");
  h2dNpart->GetXaxis()->SetTitle(Form("%s percentile", estimator));
  h2dNcoll->GetXaxis()->SetTitle(Form("%s percentile", estimator));
  h2dNpart->GetYaxis()->SetTitleSize(0.06);
  h2dNcoll->GetYaxis()->SetTitleSize(0.06);
  h2dNpart->GetYaxis()->SetTitleOffset(1.05);
  h2dNcoll->GetYaxis()->SetTitleOffset(1.05);
  h2dNpart->GetXaxis()->SetTitleSize(0.05);
  h2dNcoll->GetXaxis()->SetTitleSize(0.05);

  TCanvas* canvNpart = new TCanvas(Form("canvNpart%s", canvSuffix), "", 800, 800);
  canvNpart->SetLeftMargin(0.13);
  canvNpart->SetRightMargin(0.13);
  canvNpart->SetBottomMargin(0.13);
  canvNpart->SetTopMargin(0.13);
  canvNpart->SetLogz();
  h2dNpart->Draw("colz");
  drawMarkerForPercentile(canvNpart, h2dNpart, tag);
  canvNpart->SaveAs(Form("h2dNpart%s.pdf", tag));
  save::canvasToRoot(canvNpart, Form("h2dNpart%s.pdf", tag));

  TCanvas* canvNcoll = new TCanvas(Form("canvNcoll%s", canvSuffix), "", 800, 800);
  canvNcoll->SetLeftMargin(0.13);
  canvNcoll->SetRightMargin(0.13);
  canvNcoll->SetBottomMargin(0.13);
  canvNcoll->SetTopMargin(0.13);
  canvNcoll->SetLogz();
  h2dNcoll->Draw("colz");
  drawMarkerForPercentile(canvNcoll, h2dNcoll, tag);
  canvNcoll->SaveAs(Form("h2dNcoll%s.pdf", tag));
  save::canvasToRoot(canvNcoll, Form("h2dNcoll%s.pdf", tag));
}

void drawNpartNcollTH2()
{
  gStyle->SetOptStat(0);

  // O-O, run 564374, collisions (canvas names kept as before: canvNpart / canvNcoll)
  drawSystem("..//AnalysisResults/LHC25ae_pass2/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root",
             "OO", "FT0M", 35, 100, 3, "");

  // Pb-Pb 2023, run 544122, BCs, anchor point 90%
  drawSystem("../AnalysisResults/LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.root",
             "PbPb_544122", "FT0C", 420, 2000, 1, "PbPb");
}
