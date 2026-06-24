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
#include <map>
#include <set>

#include "colorManager.h"

static std::map <int, const char*> calibs = {
  {567905, "AR_567905_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567911, "AR_567911_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567912, "AR_567912_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567925, "AR_567925_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567927, "AR_567927_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567939, "AR_567939_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567966, "AR_567966_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567967, "AR_567967_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567986, "AR_567986_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {567999, "AR_567999_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568012, "AR_568012_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568026, "AR_568026_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568039, "AR_568039_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568040, "AR_568040_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568041, "AR_568041_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568067, "AR_568067_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568080, "AR_568080_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568117, "AR_568117_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568146, "AR_568146_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568191, "AR_568191_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568209, "AR_568209_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568210, "AR_568210_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568211, "AR_568211_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568227, "AR_568227_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568228, "AR_568228_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568229, "AR_568229_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568242, "AR_568242_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568256, "AR_568256_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568269, "AR_568269_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568282, "AR_568282_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568294, "AR_568294_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568321, "AR_568321_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568322, "AR_568322_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568365, "AR_568365_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568367, "AR_568367_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568380, "AR_568380_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568393, "AR_568393_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568409, "AR_568409_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568421, "AR_568421_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568422, "AR_568422_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568451, "AR_568451_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568452, "AR_568452_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568467, "AR_568467_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568468, "AR_568468_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568484, "AR_568484_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568486, "AR_568486_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568500, "AR_568500_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568521, "AR_568521_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568522, "AR_568522_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568524, "AR_568524_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568525, "AR_568525_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568526, "AR_568526_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568538, "AR_568538_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568539, "AR_568539_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568568, "AR_568568_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568586, "AR_568586_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568587, "AR_568587_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568600, "AR_568600_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568601, "AR_568601_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568602, "AR_568602_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568616, "AR_568616_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568617, "AR_568617_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568630, "AR_568630_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568631, "AR_568631_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568646, "AR_568646_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568647, "AR_568647_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568648, "AR_568648_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568663, "AR_568663_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568664, "AR_568664_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568678, "AR_568678_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568694, "AR_568694_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568695, "AR_568695_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568709, "AR_568709_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"},
  {568721, "AR_568721_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root"}
};

static std::map<int, int> runToFill = {
  {568500, 11361},
  {568486, 11359},
  {568484, 11359},
  {568468, 11357},
  {568467, 11357},
  {568452, 11356},
  {568451, 11356},
  {568422, 11353},
  {568421, 11353},
  {568539, 11366},
  {568538, 11366},
  {568526, 11365},
  {568525, 11365},
  {568524, 11365},
  {568522, 11365},
  {568521, 11365},
  {568322, 11343},
  {568321, 11343},
  {568294, 11338},
  {568282, 11335},
  {568269, 11334},
  {568256, 11333},
  {568242, 11332},
  {568229, 11331},
  {568228, 11331},
  {568227, 11331},
  {568211, 11330},
  {568210, 11330},
  {568209, 11330},
  {568191, 11328},
  {568146, 11323},
  {568117, 11318},
  {568080, 11316},
  {568067, 11314},
  {568041, 11310},
  {568040, 11310},
  {568039, 11310},
  {568026, 11307},
  {568012, 11305},
  {567999, 11304},
  {567986, 11303},
  {567967, 11297},
  {567966, 11297},
  {567939, 11292},
  {567927, 11290},
  {567926, 11290},
  {567925, 11290},
  {567912, 11289},
  {567911, 11289},
  {567905, 11289},
  {567685, 11273},
  {567684, 11273},
  {567683, 11273},
  {567682, 11273},
  {567681, 11273},
  {567680, 11273},
  {567679, 11273},
  {567678, 11273},
  {567677, 11273},
  {567676, 11273},
  {568721, 11391},
  {568709, 11390},
  {568695, 11389},
  {568694, 11389},
  {568678, 11388},
  {568664, 11387},
  {568663, 11387},
  {568648, 11384},
  {568647, 11384},
  {568646, 11384},
  {568631, 11383},
  {568630, 11383},
  {568617, 11382},
  {568616, 11382},
  {568602, 11381},
  {568601, 11381},
  {568600, 11381},
  {568587, 11380},
  {568586, 11380},
  {568568, 11377},
  {568409, 11352},
  {568393, 11350},
  {568380, 11349},
  {568367, 11348},
  {568366, 11348},
  {568365, 11348}
};

std::array<Double_t, 2> calculateAvgNpart(const char* path, Int_t bin)
{
  const char* filePath = Form("/home/jesgum/alice/Centrality/AnalysisResults/LHC25_pass1/%s", path);
  std::unique_ptr<TFile> infile(TFile::Open(filePath, "READ"));
  if (!infile || infile->IsZombie()) {
    std::cerr << "Could not open " << filePath << std::endl;
    return {};
  }

  TH2F* h2dNpart = dynamic_cast<TH2F*>(infile->Get("h2dNpart"));
  if (!h2dNpart) {
    return {};
  }

  TH1D* proj = h2dNpart->ProjectionY(Form("hProjectionNpart_%d", bin), bin, bin);
  return {proj->GetMean(), proj->GetMeanError()};
}

std::array<Double_t, 2> calculateAvgNcoll(const char* path, Int_t bin)
{
  const char* filePath = Form("/home/jesgum/alice/Centrality/AnalysisResults/LHC25_pass1/%s", path);
  std::unique_ptr<TFile> infile(TFile::Open(filePath, "READ"));
  if (!infile || infile->IsZombie()) {
    std::cerr << "Could not open " << filePath << std::endl;
    return {};
  }

  TH2F* h2dNcoll = dynamic_cast<TH2F*>(infile->Get("h2dNcoll"));
  if (!h2dNcoll) {
    return {};
  }

  TH1D* proj = h2dNcoll->ProjectionY(Form("hProjectionNcoll_%d", bin), bin, bin);
  return {proj->GetMean(), proj->GetMeanError()};
}

void drawNpartNcollVsRun()
{
  gStyle->SetOptStat(0);
  std::vector<int> fillOrder;
  {
    std::set<int> seen;
    for (auto& kv : runToFill) {
      if (seen.insert(kv.second).second) {
        fillOrder.push_back(kv.second);
      }
    }
  }

  std::vector<int> palette = {
    kRed+1, kBlue+1, kGreen+2, kOrange+1, kMagenta+1,
    kCyan+2, kYellow+2, kViolet+1, kTeal+2, kPink+1,
    kSpring+5, kAzure+2, kOrange-3, kRed-7, kBlue-7
  };

  std::map<int, int> fillColor;
  for (size_t i = 0; i < fillOrder.size(); ++i) {
    fillColor[fillOrder[i]] = palette[i % palette.size()];
  }

  for (Int_t binCent{0}; binCent < 100; ++binCent) {
    TH1F* hNpart = new TH1F(Form("hNpart%d", binCent), "", calibs.size(), 0, calibs.size() - 1);
    hNpart->SetMarkerStyle(kFullCircle);
    hNpart->SetMarkerColor(kBlack);
    hNpart->SetLineColor(kBlack);
    TH1F* hNcoll = new TH1F(Form("hNcoll%d", binCent), "", calibs.size(), 0, calibs.size() - 1);
    hNcoll->SetMarkerStyle(kFullCircle);
    hNcoll->SetMarkerColor(kBlack);
    hNcoll->SetLineColor(kBlack);

    TLegend* legend = new TLegend();
    Int_t binRun{1};
    Double_t avgNpart{0}, avgNColl{0};
    std::vector<TBox*> fillBoxes;
    for (const auto& [runNumber, path] : calibs) {
      std::array<Double_t, 2> runNpart = calculateAvgNpart(path, binCent + 1);
      std::array<Double_t, 2> runNcoll = calculateAvgNcoll(path, binCent + 1);
      hNpart->GetXaxis()->SetBinLabel(binRun, Form("%d", runNumber));
      hNcoll->GetXaxis()->SetBinLabel(binRun, Form("%d", runNumber));
      hNpart->SetBinContent(binRun, runNpart[0]);
      hNcoll->SetBinContent(binRun, runNcoll[0]);
      hNpart->SetBinError(binRun, runNpart[1]);
      hNcoll->SetBinError(binRun, runNcoll[1]);
      avgNpart += runNpart[0];
      avgNColl += runNcoll[0];
      TBox* box = new TBox(hNcoll->GetXaxis()->GetBinLowEdge(binRun), 0, hNcoll->GetXaxis()->GetBinUpEdge(binRun), 1);
      box->SetFillColorAlpha(fillColor[runToFill[runNumber]], 0.2);
      box->SetLineWidth(0);
      fillBoxes.push_back(box);
      ++binRun;
    }

    avgNpart /= calibs.size();
    avgNColl /= calibs.size();

    TLine* lineNpart = new TLine(0, avgNpart, calibs.size() - 1, avgNpart);
    TLine* lineNcoll = new TLine(0, avgNColl, calibs.size() - 1, avgNColl);

    TLatex* latex = new TLatex();
    latex->SetTextAlign(33);
    TCanvas* canvNpart = new TCanvas(Form("canvNpart_%d", binCent), "", 1600, 1000);
    canvNpart->SetTicks(1, 1);
    canvNpart->SetLogz();
    canvNpart->SetTopMargin(0.05);
    canvNpart->SetLeftMargin(0.1);
    canvNpart->SetRightMargin(0.03);
    canvNpart->SetBottomMargin(0.13);
    hNpart->GetYaxis()->SetTitle("<N_{part}>");
    hNpart->SetMaximum(avgNpart + 5);
    hNpart->SetMinimum(avgNpart - 5);
    // hNpart->SetMaximum(hNpart->GetMaximum() * 1.001);
    // hNpart->SetMinimum(hNpart->GetMinimum() * 0.999);
    hNpart->Draw();
    for (auto& box : fillBoxes) {
      box->SetY1(avgNpart - 5);
      box->SetY2(avgNpart + 5);
      box->Draw("same");
    }
    latex->DrawTextNDC(0.97, 0.995, Form("%d-%d %%", binCent, binCent + 1));
    canvNpart->SaveAs(Form("FiguresPerCentralityPercentile/hNpart_%d.pdf", binCent));
    canvNpart->SaveAs(Form("FiguresPerCentralityPercentile/hNpart_%d.png", binCent));

    TCanvas* canvNcoll = new TCanvas(Form("canvNcoll_%d", binCent), "", 1600, 1000);
    canvNcoll->SetTicks(1, 1);
    canvNcoll->SetLogz();
    canvNcoll->SetTopMargin(0.05);
    canvNcoll->SetLeftMargin(0.1);
    canvNcoll->SetRightMargin(0.03);
    canvNcoll->SetBottomMargin(0.13);
    hNcoll->GetYaxis()->SetTitle("<N_{coll}>");
    // hNcoll->SetMaximum(hNcoll->GetMaximum() * 1.001);
    // hNcoll->SetMinimum(hNcoll->GetMinimum() * 0.999);
    hNcoll->SetMaximum(avgNColl + 5);
    hNcoll->SetMinimum(avgNColl - 5);
    hNcoll->Draw();
    for (auto& box : fillBoxes) {
      box->SetY1(avgNColl - 5);
      box->SetY2(avgNColl + 5);
      box->Draw("same");
    }
    latex->DrawTextNDC(0.97, 0.995, Form("%d-%d %%", binCent, binCent + 1));
    canvNcoll->SaveAs(Form("FiguresPerCentralityPercentile/hNcoll_%d.pdf", binCent));    
    canvNcoll->SaveAs(Form("FiguresPerCentralityPercentile/hNcoll_%d.png", binCent));    
    delete canvNpart; delete canvNcoll;
  }
}