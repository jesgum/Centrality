#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <string>

#include "colorManager.h"

struct GlauberValue {
  double value;
  double error;
};

std::vector<GlauberValue> ooNpart = {
  { 26.34, 0.59 },
  { 24.08, 0.69 },
  { 21.03, 0.84 },
  { 16.88, 0.89 },
  { 13.05, 0.75 },
  { 9.88, 0.61 },
  { 7.39, 0.48 },
  { 5.51, 0.36 },
  { 4.15, 0.26 },
  { 3.17, 0.23 },
  { 2.42, 0.16 }
};

std::vector<GlauberValue> ooNcoll = {
  { 42.31, 2.61 },
  { 35.82, 2.70 },
  { 28.77, 2.73 },
  { 20.58, 2.35 },
  { 14.19, 1.75 },
  { 9.64, 1.22 },
  { 6.51, 0.82 },
  { 4.41, 0.54 },
  { 3.00, 0.35 },
  { 2.07, 0.26 },
  { 1.38, 0.16 }
};

std::vector<TH1F*> doNpartNcoll(const char* name, const int ap)
{
  Bool_t fixedMu = true;
  Bool_t fixedK = true;
  Int_t ancesterMode = 2;
  Bool_t printNpartNcoll = false;

  TString dataset = "LHC25ae_pass2";
  TString ar = "AR_564374";
  TString filePath = Form("../AnalysisResults/%s/%s_calibration_%s_Anchor%d.root", dataset.Data(), ar.Data(), name, ap);
  TFile* infile = new TFile(filePath.Data(), "read");

  if (!infile || infile->IsZombie()) {
    std::cerr << "Could not find " << filePath << std::endl;
    return {};
  }

  TH2F* h2dNpart = dynamic_cast<TH2F*>(infile->Get("h2dNpart"));
  TH2F* h2dNcoll = dynamic_cast<TH2F*>(infile->Get("h2dNcoll"));

  if (!h2dNpart) {
    std::cerr << "Error: h2dNpart not found or wrong type!" << std::endl;
    return {};
  }

  if (!h2dNcoll) {
    std::cerr << "Error: h2dNcoll not found or wrong type!" << std::endl;
    return {};
  }

  h2dNpart->GetXaxis()->SetTitle("Centrality percentile (FT0M)");
  h2dNcoll->GetXaxis()->SetTitle("Centrality percentile (FT0M)");
  h2dNpart->GetXaxis()->SetTitleSize(0.045);
  h2dNcoll->GetXaxis()->SetTitleSize(0.045);
  h2dNpart->GetXaxis()->SetTitleOffset(1.05);
  h2dNcoll->GetXaxis()->SetTitleOffset(1.05);
  h2dNpart->GetYaxis()->SetTitle("N_{part}");
  h2dNcoll->GetYaxis()->SetTitle("N_{coll}");
  h2dNpart->GetYaxis()->SetTitleSize(0.05);
  h2dNcoll->GetYaxis()->SetTitleSize(0.05);
  h2dNpart->GetYaxis()->SetTitleOffset(1);
  h2dNcoll->GetYaxis()->SetTitleOffset(1);

  TH1F* hAvgNpart = new TH1F("hAvgNpart_local", "", 11, 0, 100);
  TH1F* hAvgNcoll = new TH1F("hAvgNcoll_local", "", 11, 0, 100);

  hAvgNpart->SetDirectory(0);
  hAvgNcoll->SetDirectory(0);

  hAvgNpart->SetMarkerStyle(kFullCircle);
  hAvgNpart->SetMarkerColor(kBlack);
  hAvgNcoll->SetMarkerStyle(kFullCircle);
  hAvgNcoll->SetMarkerColor(kBlack);

  if (printNpartNcoll) {
    std::cout << " - - - cent - - - <Npart> - - - <Ncoll> - - - " << std::endl;
  }

  const std::vector<std::pair<int, int>> centRanges = {
    { 0, 5 },
    { 5, 10 },
    { 10, 20 },
    { 20, 30 },
    { 30, 40 },
    { 40, 50 },
    { 50, 60 },
    { 60, 70 },
    { 70, 80 },
    { 80, 90 },
    { 90, 100 }
  };

  for (size_t icent = 0; icent < centRanges.size(); ++icent) {

    int centMin = centRanges[icent].first;
    int centMax = centRanges[icent].second;

    // Convert percentile to histogram bins
    int firstBin = centMin + 1;
    int lastBin = centMax;

    TH1D* hProjectionNpart = dynamic_cast<TH1D*>(
      h2dNpart->ProjectionY(
        Form("hProjectionNpart_%zu", icent),
        firstBin,
        lastBin));

    TH1D* hProjectionNcoll = dynamic_cast<TH1D*>(
      h2dNcoll->ProjectionY(
        Form("hProjectionNcoll_%zu", icent),
        firstBin,
        lastBin));

    hAvgNpart->SetBinContent(icent + 1, hProjectionNpart->GetMean());
    hAvgNcoll->SetBinContent(icent + 1, hProjectionNcoll->GetMean());
    hAvgNpart->SetBinError(icent + 1, hProjectionNpart->GetMeanError());
    hAvgNcoll->SetBinError(icent + 1, hProjectionNcoll->GetMeanError());

    if (printNpartNcoll) {
      std::cout
        << centMin << "-" << centMax << "% | "
        << hProjectionNpart->GetMean() << " ± "
        << hProjectionNpart->GetMeanError() << " | "
        << hProjectionNcoll->GetMean() << " ± "
        << hProjectionNcoll->GetMeanError()
        << std::endl;
    }

    delete hProjectionNpart;
    delete hProjectionNcoll;
  }

  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 800);
  canvNpart->SetLogz();
  canvNpart->SetTopMargin(0.03);
  canvNpart->SetLeftMargin(0.1);
  canvNpart->SetRightMargin(0.1);
  canvNpart->SetBottomMargin(0.13);
  h2dNpart->GetYaxis()->SetRangeUser(0, 35);
  h2dNpart->Draw("colz");
  hAvgNpart->Draw("pe same");
  canvNpart->SaveAs(Form("h2dNpart_%sOO.pdf", name));

  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 800);
  canvNcoll->SetLogz();
  canvNcoll->SetTopMargin(0.03);
  canvNcoll->SetLeftMargin(0.1);
  canvNcoll->SetRightMargin(0.1);
  canvNcoll->SetBottomMargin(0.13);
  h2dNcoll->GetYaxis()->SetRangeUser(0, 100);
  h2dNcoll->Draw("colz");
  hAvgNcoll->Draw("pe same");
  canvNcoll->SaveAs(Form("h2dNcoll_%sOO.pdf", name));

  delete canvNpart;
  delete canvNcoll;
  delete h2dNpart;
  delete h2dNcoll;
  infile->Close();
  delete infile;
  std::cout << std::endl;
  return { hAvgNpart, hAvgNcoll };
}

struct Histos {
  Histos(std::vector<TH1F*> vec, const char* c)
  {
    histos = vec;
    name = c;
  }
  void setColor(int color)
  {
    for (auto& hist : histos) {
      hist->SetMarkerColor(color);
      hist->SetLineColor(color);
    }
  }
  void setMarker(int marker)
  {
    for (auto& hist : histos) {
      hist->SetMarkerStyle(marker);
    }
  }

  const char* name;
  std::vector<TH1F*> histos;
};

void drawLightIonNpartNcollVsAnchorPoint()
{
  int cfgDraw = 0; // 0: BCs, 1: Collisions
  TH1F* hAvgNpart = new TH1F("hAvgNpart_global", "", 11, 0, 100);
  hAvgNpart->SetLineColor(kBlack);
  hAvgNpart->SetMarkerColor(kBlack);
  hAvgNpart->SetMarkerStyle(kFullCircle);
  hAvgNpart->GetXaxis()->SetLabelSize(0);
  hAvgNpart->GetYaxis()->SetLabelSize(0.045);
  hAvgNpart->GetYaxis()->SetTitleSize(0.05);
  hAvgNpart->GetYaxis()->SetTitle("Oxygen#minusOxygen <N_{part}>");
  hAvgNpart->GetYaxis()->SetTitleOffset(0);
  for (int bin = 0; bin < ooNpart.size(); ++bin) {
    hAvgNpart->SetBinContent(bin + 1, ooNpart[bin].value);
    hAvgNpart->SetBinError(bin + 1, ooNpart[bin].error);
  }

  TH1F* hAvgNcoll = new TH1F("hAvgNcoll_global", "", 11, 0, 100);
  hAvgNcoll->SetLineColor(kBlack);
  hAvgNcoll->SetMarkerColor(kBlack);
  hAvgNcoll->SetMarkerStyle(kFullCircle);
  hAvgNcoll->GetXaxis()->SetLabelSize(0);
  hAvgNcoll->GetYaxis()->SetLabelSize(0.045);
  hAvgNcoll->GetYaxis()->SetTitleSize(0.05);
  hAvgNcoll->GetYaxis()->SetTitle("Oxygen#minusOxygen <N_{coll}>");
  hAvgNcoll->GetYaxis()->SetTitleOffset(0);
  for (int bin = 0; bin < ooNpart.size(); ++bin) {
    hAvgNcoll->SetBinContent(bin + 1, ooNcoll[bin].value);
    hAvgNcoll->SetBinError(bin + 1, ooNcoll[bin].error);
  }

  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  std::vector<Histos*> bcs, col;
  std::vector<int> anchorpointPercentage = { 95, 90, 80, 70, 60, 50 };

  for (const auto& ap : anchorpointPercentage) {
    bcs.push_back(new Histos(doNpartNcoll("hFT0M_BCs", ap), Form("%d", ap)));
    col.push_back(new Histos(doNpartNcoll("hFT0M_Collisions", ap), Form("%d", ap)));
  }

  ColorManager cm(anchorpointPercentage.size());
  for (int ii{ 0 }; ii < anchorpointPercentage.size(); ++ii) {
    bcs[ii]->setColor(cm.getColor(ii));
    col[ii]->setColor(cm.getColor(ii));
    bcs[ii]->setMarker(kFullCircle);
    col[ii]->setMarker(kFullCircle);
  }

  // === Canvas: Npart with ratio ===
  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 800);
  canvNpart->SetTopMargin(0.0);
  canvNpart->SetBottomMargin(0.0);

  TPad* padNpartTop = new TPad("padNpartTop", "", 0, 0.3, 1, 1);
  padNpartTop->SetTicks(1, 1);
  padNpartTop->SetTopMargin(0.04);
  padNpartTop->SetBottomMargin(0.02);
  padNpartTop->SetLeftMargin(0.1);
  padNpartTop->SetRightMargin(0.03);
  padNpartTop->Draw();
  padNpartTop->cd();
  hAvgNpart->Draw("pe");
  for (int ii{ 0 }; ii < anchorpointPercentage.size(); ++ii) {
    if (cfgDraw == 0) bcs[ii]->histos[kNpart]->Draw("pe same");
    if (cfgDraw == 1) col[ii]->histos[kNpart]->Draw("pe same");
  }

  TLegend* legAnchor = new TLegend(0.8, 0.48, 1.12, 0.94);
  legAnchor->SetTextFont(42);
  legAnchor->SetBorderSize(0);
  legAnchor->SetFillColorAlpha(0, 0);
  legAnchor->AddEntry(hAvgNpart, "Table", "pl");
  for (int ii{ 0 }; ii < anchorpointPercentage.size(); ++ii) {
    legAnchor->AddEntry(bcs[ii]->histos[kNpart], Form("%d%%", anchorpointPercentage[ii]), "pl");
  }

  TLegend* legHist = new TLegend(0.7, 0.77, 0.91, 0.94);
  TH1F* hFullCircle = new TH1F();
  TH1F* hFullSquare = new TH1F();
  hFullCircle->SetMarkerStyle(kFullCircle);
  hFullCircle->SetMarkerSize(2);
  hFullCircle->SetMarkerColor(kBlack);
  hFullSquare->SetMarkerStyle(kFullSquare);
  hFullSquare->SetMarkerSize(2);
  hFullSquare->SetMarkerColor(kBlack);
  legHist->SetTextFont(42);
  legHist->SetBorderSize(0);
  legHist->SetFillColorAlpha(0, 0);
  legHist->AddEntry(hFullCircle, "Coll", "p");
  legHist->AddEntry(hFullSquare, "BCs", "p");

  legAnchor->Draw();
  // legHist->Draw();
  canvNpart->cd();
  TPad* padNpartBot = new TPad("padNpartBot", "", 0, 0, 1, 0.3);
  padNpartBot->SetTicks(1, 1);
  padNpartBot->SetTopMargin(0.02);
  padNpartBot->SetBottomMargin(0.35);
  padNpartBot->SetLeftMargin(0.1);
  padNpartBot->SetRightMargin(0.03);
  padNpartBot->Draw();
  padNpartBot->cd();

  bool first = true;
  for (int ii{ 0 }; ii < anchorpointPercentage.size(); ++ii) {
    TH1F* hRatioBcs = (TH1F*)bcs[ii]->histos[kNpart]->Clone(Form("hPartRatioBcs_%d", ii));
    TH1F* hRatioCol = (TH1F*)col[ii]->histos[kNpart]->Clone(Form("hPartRatioCol_%d", ii));
    hRatioBcs->Divide(hRatioBcs, hAvgNpart, 1, 1, "B");
    hRatioCol->Divide(hRatioCol, hAvgNpart, 1, 1, "B");

    if (first) {
      hRatioBcs->SetTitle("");
      hRatioBcs->GetYaxis()->SetTitle("");
      hRatioBcs->GetYaxis()->SetNdivisions(505);
      hRatioBcs->GetYaxis()->SetTitleSize(0.12);
      hRatioBcs->GetYaxis()->SetRangeUser(0.89, 1.11);
      hRatioBcs->GetYaxis()->SetTitleOffset(0.4);
      hRatioBcs->GetXaxis()->SetTitle("FT0M Amplitude");
      hRatioBcs->GetYaxis()->SetLabelSize(0.10);
      hRatioBcs->GetXaxis()->SetTitleSize(0.13);
      hRatioBcs->GetXaxis()->SetLabelSize(0.11);
      hRatioCol->SetTitle("");
      hRatioCol->GetYaxis()->SetTitle("");
      hRatioCol->GetYaxis()->SetNdivisions(505);
      hRatioCol->GetYaxis()->SetTitleSize(0.12);
      hRatioCol->GetYaxis()->SetRangeUser(0.89, 1.11);
      hRatioCol->GetYaxis()->SetTitleOffset(0.4);
      hRatioCol->GetXaxis()->SetTitle("FT0M Amplitude");
      hRatioCol->GetYaxis()->SetLabelSize(0.10);
      hRatioCol->GetXaxis()->SetTitleSize(0.13);
      hRatioCol->GetXaxis()->SetLabelSize(0.11);
      if (cfgDraw == 0) hRatioBcs->Draw("pe");
      if (cfgDraw == 1) hRatioCol->Draw("pe");
      first = false;
    } else {
      if (cfgDraw == 0) hRatioBcs->Draw("pe same");
      if (cfgDraw == 1) hRatioCol->Draw("pe same");
    }
  }
  TLine* line = new TLine(0, 1, 100, 1);
  line->SetLineColor(kGray+2);
  line->SetLineStyle(7);
  line->Draw();

  // === Canvas: Ncoll with ratio ===
  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 800);
  canvNcoll->SetTopMargin(0.0);
  canvNcoll->SetBottomMargin(0.0);

  TPad* padNcollTop = new TPad("padNcollTop", "", 0, 0.3, 1, 1);
  padNcollTop->SetTicks(1, 1);
  padNcollTop->SetTopMargin(0.04);
  padNcollTop->SetBottomMargin(0.02);
  padNcollTop->SetLeftMargin(0.1);
  padNcollTop->SetRightMargin(0.03);
  padNcollTop->Draw();
  padNcollTop->cd();
  hAvgNcoll->Draw("pe");
  for (int ii{ 0 }; ii < anchorpointPercentage.size(); ++ii) {
    if (cfgDraw == 0) bcs[ii]->histos[kNcoll]->Draw("pe same");
    if (cfgDraw == 1) col[ii]->histos[kNcoll]->Draw("pe same");
  }
  legAnchor->Draw();

  canvNcoll->cd();
  TPad* padNcollBot = new TPad("padNcollBot", "", 0, 0, 1, 0.3);
  padNcollBot->SetTicks(1, 1);
  padNcollBot->SetTopMargin(0.02);
  padNcollBot->SetBottomMargin(0.35);
  padNcollBot->SetLeftMargin(0.1);
  padNcollBot->SetRightMargin(0.03);
  padNcollBot->Draw();
  padNcollBot->cd();

  first = true;
  for (int ii{ 0 }; ii < anchorpointPercentage.size(); ++ii) {
    TH1F* hRatioBcs = (TH1F*)bcs[ii]->histos[kNcoll]->Clone(Form("hCollRatioBcs_%d", ii));
    TH1F* hRatioCol = (TH1F*)col[ii]->histos[kNcoll]->Clone(Form("hCollRatioCol_%d", ii));
    hRatioBcs->Divide(hRatioBcs, hAvgNcoll, 1, 1, "B");
    hRatioCol->Divide(hRatioCol, hAvgNcoll, 1, 1, "B");

    if (first) {
      hRatioBcs->SetTitle("");
      hRatioBcs->GetYaxis()->SetTitle("");
      hRatioBcs->GetYaxis()->SetNdivisions(505);
      hRatioBcs->GetYaxis()->SetTitleSize(0.12);
      hRatioBcs->GetYaxis()->SetRangeUser(0.81, 1.19);
      hRatioBcs->GetYaxis()->SetTitleOffset(0.4);
      hRatioBcs->GetXaxis()->SetTitle("FT0M Amplitude");
      hRatioBcs->GetYaxis()->SetLabelSize(0.10);
      hRatioBcs->GetXaxis()->SetTitleSize(0.13);
      hRatioBcs->GetXaxis()->SetLabelSize(0.11);
      hRatioCol->SetTitle("");
      hRatioCol->GetYaxis()->SetTitle("");
      hRatioCol->GetYaxis()->SetNdivisions(505);
      hRatioCol->GetYaxis()->SetTitleSize(0.12);
      hRatioCol->GetYaxis()->SetRangeUser(0.81, 1.19);
      hRatioCol->GetYaxis()->SetTitleOffset(0.4);
      hRatioCol->GetXaxis()->SetTitle("FT0M Amplitude");
      hRatioCol->GetYaxis()->SetLabelSize(0.10);
      hRatioCol->GetXaxis()->SetTitleSize(0.13);
      hRatioCol->GetXaxis()->SetLabelSize(0.11);
      if (cfgDraw == 0) hRatioBcs->Draw("pe");
      if (cfgDraw == 1) hRatioCol->Draw("pe");
      first = false;
    } else {
      if (cfgDraw == 0) hRatioBcs->Draw("pe same");
      if (cfgDraw == 1) hRatioCol->Draw("pe same");
    }
  }
  line->Draw();
  if (cfgDraw == 0) {
    canvNpart->SaveAs("hBCsAnchoredNpartOO.pdf");
    canvNcoll->SaveAs("hBCsAnchoredNcollOO.pdf");
  }
  if (cfgDraw == 1) {
    canvNpart->SaveAs("hColAnchoredNpartOO.pdf");
    canvNcoll->SaveAs("hColAnchoredNcollOO.pdf");
  }
}