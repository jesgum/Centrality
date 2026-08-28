#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <string>

#include "Includes/colorManager.h"

static bool doNeNe = false;

struct GlauberInfo {
  TFile* file = nullptr;
  GlauberInfo(const char* hist, const int anchorPoint, const int color)
  {
    const char* pathOO = Form("../AnalysisResults/LHC25ae_pass2/AR_564374_calibration_%s_Anchor%d.root", hist, anchorPoint);
    const char* pathNeNe = Form("../AnalysisResults/LHC25af_pass2_systematics/AR_564468_calibration_%s_Anchor%d.root", hist, anchorPoint);
    const char* path = doNeNe ? pathNeNe : pathOO;
    file = TFile::Open(path, "READ");
    if (!file || file->IsZombie()) {
      throw std::runtime_error(std::string("Cannot open file: ") + path);
    }

    hData = (TH1F*)file->Get("hData");
    hGlauber = (TH1F*)file->Get("hGlauberFine");
    TH1F* hAnchorPoint = (TH1F*)file->Get("hAnchorPoint");
    if (!hAnchorPoint) {
      throw std::runtime_error(Form("Cannot open histogram: hAnchorPoint for file: %s", path));
    }

    mAnchorPoint = hAnchorPoint->GetBinContent(1);
    delete hAnchorPoint;

    mColor = color;

    if (!hData) {
      throw std::runtime_error(std::string("Cannot open histogram: hData"));
    }

    if (!hGlauber) {
      throw std::runtime_error(std::string("Cannot open histogram: hGlauberFine"));
    }

    hGlauber->SetLineColor(color);

    hInsetData = (TH1F*)hData->Clone(Form("hInsetData_%s", path));
    hInsetGlauber = (TH1F*)hGlauber->Clone(Form("hInsetGlauber_%s", path));
    hRatio = (TH1F*)hGlauber->Clone(Form("hRatio_%s", path));
    hInsetRatio = (TH1F*)hInsetGlauber->Clone(Form("hInsetRatio_%s", path));

    hRatio->SetDirectory(nullptr);
    hRatio->SetTitle("");
    hInsetData->SetDirectory(nullptr);
    hInsetData->SetTitle("");
    hInsetGlauber->SetDirectory(nullptr);
    hInsetGlauber->SetTitle("");
    hInsetRatio->SetDirectory(nullptr);
    hInsetRatio->SetTitle("");

    TH1F* hDataClone = (TH1F*)hData->Clone("hDataClone");
    TH1F* hInsetGlauberClone = (TH1F*)hInsetGlauber->Clone("hInsetGlauberClone");

    const int rebin = 20;
    hRatio->Rebin(rebin);
    hDataClone->Rebin(rebin);
    hInsetGlauberClone->Rebin(rebin);

    hRatio->Divide(hDataClone, hInsetGlauberClone, 1, 1, "B");
    hInsetRatio->Divide(hInsetData, hInsetGlauber, 1, 1, "B");
    delete hDataClone;
    delete hInsetGlauberClone;
  }

  ~GlauberInfo()
  {
    if (file) {
      file->Close();
      delete file;
      hData = nullptr;
      hGlauber = nullptr;
    }
  }

  GlauberInfo(GlauberInfo&& other) noexcept
    : file(other.file),
      hData(other.hData),
      hGlauber(other.hGlauber),
      hRatio(other.hRatio),
      hInsetData(other.hInsetData),
      hInsetGlauber(other.hInsetGlauber),
      hInsetRatio(other.hInsetRatio)
  {
    other.file = nullptr;
    other.hData = nullptr;
    other.hGlauber = nullptr;
    other.hRatio = nullptr;
    other.hInsetData = nullptr;
    other.hInsetGlauber = nullptr;
    other.hInsetRatio = nullptr;
  }

  GlauberInfo& operator=(GlauberInfo&&) = default;

  int mAnchorPoint = 0;
  int mColor = 0;

  TH1F* hData = nullptr;
  TH1F* hGlauber = nullptr;
  TH1F* hRatio = nullptr;

  TH1F* hInsetData = nullptr;
  TH1F* hInsetGlauber = nullptr;
  TH1F* hInsetRatio = nullptr;
};

void drawGlauberSet(TCanvas* canv, std::vector<GlauberInfo>& info, const char* title, std::vector<int> anchor)
{
  const float xmax = doNeNe ? 32000 : 26000;
  const float ymax = 5e+6;
  const float ratio_ymin = 0.6;
  const float ratio_ymax = 2.5;

  canv->cd(1);
  TLegend* leg = new TLegend(0.76, 0.5, 0.98, 0.96);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0, 0);

  bool first = true;
  for (int ii{ 0 }; ii < info.size(); ++ii) {
    auto& g = info[ii];
    if (first) {
      g.hData->SetLineColor(kBlack);
      g.hData->GetXaxis()->SetRangeUser(0, xmax);
      g.hData->GetYaxis()->SetRangeUser(0.25, ymax);
      g.hData->SetTitle(title);
      g.hData->GetYaxis()->SetTitle("Count");
      g.hData->GetYaxis()->SetTitleSize(0.065);
      g.hData->GetYaxis()->SetTitleOffset(0.5);
      g.hData->GetYaxis()->SetLabelSize(0.055);
      g.hData->Draw("hist");
      leg->AddEntry(g.hData, "Data", "l");
      first = false;
    }
    g.hGlauber->Draw("HIST SAME");
    leg->AddEntry(g.hGlauber, Form("Anchor %d", anchor[ii]), "l");
  }

  leg->Draw();


  float width = 0.6 - 0.25;
  float height = 0.9 - 0.4;
  float xpos = 0.1;
  float ypos = 0.1;

  TPad* padTop = new TPad("padTop", "", xpos, ypos, xpos+width, ypos+height);
  padTop->SetLogy();
  padTop->SetTopMargin(0.0);
  padTop->SetRightMargin(0.03);
  padTop->SetFillStyle(4000);
  padTop->SetFrameFillStyle(4000);
  padTop->SetLeftMargin(0.08);
  padTop->SetBottomMargin(0.15);
  padTop->SetTicks(1, 1);
  padTop->Draw();
  padTop->cd();

  first = true;
  for (int ii{ 0 }; ii < info.size(); ++ii) {
    auto& g = info[ii];
    if (first) {
      g.hInsetData->SetLineColor(kBlack);
      g.hInsetData->GetXaxis()->SetRangeUser(0, 4000);
      g.hInsetData->GetYaxis()->SetRangeUser(3e+4, 5e+5);
      g.hInsetData->SetTitle(title);
      g.hInsetData->GetYaxis()->SetTitle("Count");
      g.hInsetData->GetYaxis()->SetTitleSize(0.065);
      g.hInsetData->GetYaxis()->SetTitleOffset(0.5);
      g.hInsetData->GetYaxis()->SetLabelSize(0.055);
      g.hInsetData->Draw("hist");
      first = false;
    }
    g.hInsetGlauber->Draw("HIST SAME");
  }

  canv->Update();
  canv->cd(2);

  first = true;
  for (int ii{ 0 }; ii < info.size(); ++ii) {
    auto& g = info[ii];
    g.hRatio->GetXaxis()->SetRangeUser(0, xmax);
    g.hRatio->GetYaxis()->SetRangeUser(ratio_ymin, ratio_ymax);

    g.hRatio->GetXaxis()->SetLabelSize(0.05);
    g.hRatio->GetYaxis()->SetLabelSize(0.05);

    g.hRatio->GetXaxis()->SetTitle("FT0M Amplitude");
    g.hRatio->GetYaxis()->SetTitle("Data / Fit");

    g.hRatio->GetXaxis()->SetTitleSize(0.06);
    g.hRatio->GetYaxis()->SetTitleSize(0.06);

    g.hRatio->GetXaxis()->SetTitleOffset(1);
    g.hRatio->GetYaxis()->SetTitleOffset(0.5);

    if (first) {
      g.hRatio->Draw("HIST");
      first = false;
    } else {
      g.hRatio->Draw("HIST SAME");
    }
  }

  for (int ii{ 0 }; ii < info.size(); ++ii) {
    auto& g = info[ii];
    TLine* anchorLine = new TLine(g.mAnchorPoint, ratio_ymin, g.mAnchorPoint, ratio_ymax);
    anchorLine->SetLineStyle(7);
    anchorLine->SetLineColorAlpha(g.mColor, 0.6);
    anchorLine->Draw();
  }

  TLine* line = new TLine(0, 1, xmax, 1);
  line->SetLineColor(kGray + 2);
  line->SetLineStyle(7);
  line->Draw();

  TLegend* leg2 = new TLegend(0.2, 0.2, 0.45, 0.26);
  leg2->SetBorderSize(0);
  leg2->SetFillColorAlpha(0, 0);
  TH1F* hAnchorLine = new TH1F();
  hAnchorLine->SetLineColor(kBlack);
  hAnchorLine->SetLineWidth(2);
  hAnchorLine->SetLineStyle(7);
  leg2->AddEntry(hAnchorLine, "Anchor point", "l");
  leg2->Draw();

  TPad* padBot = new TPad("padBot", "", 0.23, 0.4, 0.58, 0.9);
  padBot->SetTopMargin(0.0);
  padBot->SetRightMargin(0.03);
  padBot->SetFillStyle(4000);
  padBot->SetFrameFillStyle(4000);
  padBot->SetLeftMargin(0.08);
  padBot->SetBottomMargin(0.15);
  padBot->SetTicks(1, 1);
  padBot->Draw();
  padBot->cd();

  first = true;
  for (int ii{ 0 }; ii < info.size(); ++ii) {
    auto& g = info[ii];
    g.hInsetRatio->GetXaxis()->SetRangeUser(0, 4000);
    g.hInsetRatio->GetYaxis()->SetRangeUser(ratio_ymin, 1.5);

    g.hInsetRatio->GetXaxis()->SetLabelSize(0.05);
    g.hInsetRatio->GetYaxis()->SetLabelSize(0.05);

    g.hInsetRatio->GetXaxis()->SetTitle("FT0M Amplitude");
    g.hInsetRatio->GetYaxis()->SetTitle("Data / Fit");

    g.hInsetRatio->GetXaxis()->SetTitleSize(0.06);
    g.hInsetRatio->GetYaxis()->SetTitleSize(0.06);

    g.hInsetRatio->GetXaxis()->SetTitleOffset(1);
    g.hInsetRatio->GetYaxis()->SetTitleOffset(0.5);

    if (first) {
      first = false;
      g.hInsetRatio->Draw("HIST");
    } else {
      g.hInsetRatio->Draw("HIST SAME");
    }
  }

  for (int ii{ 0 }; ii < info.size(); ++ii) {
    auto& g = info[ii];
    TLine* anchorLine = new TLine(g.mAnchorPoint, ratio_ymin, g.mAnchorPoint, ratio_ymax);
    anchorLine->SetLineStyle(7);
    anchorLine->SetLineColorAlpha(g.mColor, 0.6);
    anchorLine->Draw();
  }

  TLine* line2 = new TLine(0, 1, 4000, 1);
  line2->SetLineColor(kGray + 2);
  line2->SetLineStyle(7);
  line2->Draw();
}

TCanvas* makeCanvas(const char* title)
{
  TCanvas* canv = new TCanvas(title, "", 1600, 1600);

  canv->SetRightMargin(0.0);
  canv->SetLeftMargin(0.06);
  canv->Divide(1, 2, 0, 0);

  // Top pad
  canv->cd(1)->SetLogy();
  canv->cd(1)->SetTopMargin(0.0);
  canv->cd(1)->SetRightMargin(0.03);
  canv->cd(1)->SetBottomMargin(0.0);
  canv->cd(1)->SetTicks(1, 1);

  // Bottom pad
  canv->cd(2)->SetTopMargin(0.0);
  canv->cd(2)->SetRightMargin(0.03);
  canv->cd(2)->SetBottomMargin(0.15);
  canv->cd(2)->SetTicks(1, 1);
  return canv;
}

void drawLightIonGlauberAnchorPoints(bool isNeNe = false)
{
  doNeNe = isNeNe;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // std::vector<int> anchorpointPercentage = {
  //   100, 90, 85, 80, 75,
  //   70, 65, 60, 55, 50, 45,
  //   40, 30, 25, 20, 15, 10, 5
  // };

  // std::vector<int> anchorpointPercentage = {
  //   50,
  //   60, 65, 70, 75,
  //   80, 85, 90, 95
  // };

  std::vector<int> anchorpointPercentage = {
    90, 80, 70, 60, 50
  };
  
  ColorManager cm(anchorpointPercentage.size());

  std::vector<GlauberInfo> bcs, col;
  bcs.reserve(anchorpointPercentage.size());
  col.reserve(anchorpointPercentage.size());

  for (size_t i = 0; i < anchorpointPercentage.size(); ++i) {
    int anchor = anchorpointPercentage[i];
    int color = cm.getColor(i);

    bcs.emplace_back("hFT0M_BCs", anchor, color);
    col.emplace_back("hFT0M_Collisions", anchor, color);
  }

  TCanvas* canvCol = makeCanvas("canvCol");
  drawGlauberSet(canvCol, col, "Collisions", anchorpointPercentage);
  canvCol->cd();
  TLatex* latCol = new TLatex();
  latCol->SetNDC();
  latCol->SetTextSize(0.04);
  latCol->SetTextAlign(32);
  latCol->DrawLatexNDC(0.97, 0.975, "Collisions");
  if (doNeNe) canvCol->SaveAs("FiguresSystematics/hGlauberCol_NeNe.pdf");
  if (!doNeNe) canvCol->SaveAs("FiguresSystematics/hGlauberCol_OO.pdf");

  TCanvas* canvBcs = makeCanvas("canvBcs");
  drawGlauberSet(canvBcs, bcs, "hFT0M_BCs", anchorpointPercentage);
  canvBcs->cd();
  TLatex* latBc = new TLatex();
  latBc->SetNDC();
  latBc->SetTextSize(0.04);
  latBc->SetTextAlign(32);
  latBc->DrawLatexNDC(0.97, 0.975, "Bunch crossings");
  if (doNeNe) canvBcs->SaveAs("FiguresSystematics/hGlauberBcs_NeNe.pdf");
  if (!doNeNe) canvBcs->SaveAs("FiguresSystematics/hGlauberBcs_OO.pdf");
}