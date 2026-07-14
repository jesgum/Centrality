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

std::vector<TH1F*> doNpartNcoll(TString dataset, TString ar, const char* name)
{
  Bool_t fixedMu = true;
  Bool_t fixedK = true;
  Int_t ancesterMode = 2;
  Bool_t printNpartNcoll = true;

  TString filePath = Form("../AnalysisResults/%s/%s_calibration_%s_Anchor90.root", dataset.Data(), ar.Data(), name);
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

  TH1F* hAvgNpart = new TH1F("hAvgNpart_local", "", 100, 0, 100);
  TH1F* hAvgNcoll = new TH1F("hAvgNcoll_local", "", 100, 0, 100);

  hAvgNpart->SetDirectory(0);
  hAvgNcoll->SetDirectory(0);
  hAvgNpart->SetMarkerStyle(kFullCircle);
  hAvgNpart->SetMarkerColor(kBlack);
  hAvgNcoll->SetMarkerStyle(kFullCircle);
  hAvgNcoll->SetMarkerColor(kBlack);

  if (printNpartNcoll) {
    std::cout << " - - - cent - - - <Npart> - - - <Ncoll> - - - " << std::endl;
  }

  for (Int_t icent{0}; icent < h2dNpart->GetNbinsX(); ++icent) {
    TH1D* hProjectionNpart = dynamic_cast<TH1D*>(h2dNpart->ProjectionY(Form("hProjectionNpart_%d", icent), icent + 1, icent + 1));
    TH1D* hProjectionNcoll = dynamic_cast<TH1D*>(h2dNcoll->ProjectionY(Form("hProjectionNcoll_%d", icent), icent + 1, icent + 1));
    hAvgNpart->SetBinContent(icent + 1, hProjectionNpart->GetMean());
    hAvgNpart->SetBinError(icent + 1, hProjectionNpart->GetMeanError());
    hAvgNcoll->SetBinContent(icent + 1, hProjectionNcoll->GetMean());
    hAvgNcoll->SetBinError(icent + 1, hProjectionNcoll->GetMeanError());

    if (printNpartNcoll) {
      std::cout << Form("%d-%d%% | ", icent, icent+1) 
                << hProjectionNpart->GetMean() << " | " 
                << hProjectionNcoll->GetMean() << std::endl;
    }    
    delete hProjectionNpart; delete hProjectionNcoll;
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

  delete canvNpart; delete canvNcoll;
  delete h2dNpart; delete h2dNcoll;
  infile->Close(); delete infile;
  std::cout << std::endl;
  return { hAvgNpart, hAvgNcoll };
}

struct Histos {
  Histos(std::vector<TH1F*> vec, const char* run)
  {
    histos = vec;
    runNumber = run;
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

  const char* runNumber;
  std::vector<TH1F*> histos;
};

void drawNpartNcollLightIon()
{
  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  // 564374
  // 564356
  // 564359
  // 564373

  Histos* ooBCs_564374 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564374", "hFT0M_BCs"), "BCs 564374");
  Histos* ooBCs_564356 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564356", "hFT0M_BCs"), "BCs 564356");
  Histos* ooBCs_564359 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564359", "hFT0M_BCs"), "BCs 564359");
  Histos* ooBCs_564373 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564373", "hFT0M_BCs"), "BCs 564373");
  Histos* ooColls_564374 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564374", "hFT0M_Collisions"), "Coll 564374");
  Histos* ooColls_564356 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564356", "hFT0M_Collisions"), "Coll 564356");
  Histos* ooColls_564359 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564359", "hFT0M_Collisions"), "Coll 564359");
  Histos* ooColls_564373 = new Histos(doNpartNcoll("LHC25ae_pass2", "AR_564373", "hFT0M_Collisions"), "Coll 564373");

  std::vector<Histos*> registry = { ooBCs_564374,
                                    ooBCs_564356,
                                    ooBCs_564359,
                                    ooBCs_564373,
                                    ooColls_564374,
                                    ooColls_564356,
                                    ooColls_564359,
                                    ooColls_564373};
                                    

  ColorManager cm(4);
  ooBCs_564374->setColor(cm.getColor(0));
  ooBCs_564374->setMarker(kFullSquare);
  ooBCs_564356->setColor(cm.getColor(1));
  ooBCs_564356->setMarker(kFullSquare);
  ooBCs_564359->setColor(cm.getColor(2));
  ooBCs_564359->setMarker(kFullSquare);
  ooBCs_564373->setColor(cm.getColor(3));
  ooBCs_564373->setMarker(kFullSquare);

  ooColls_564374->setColor(cm.getColor(0));
  ooColls_564374->setMarker(kFullCircle);
  ooColls_564356->setColor(cm.getColor(1));
  ooColls_564356->setMarker(kFullCircle);
  ooColls_564359->setColor(cm.getColor(2));
  ooColls_564359->setMarker(kFullCircle);
  ooColls_564373->setColor(cm.getColor(3));
  ooColls_564373->setMarker(kFullCircle);


  TH1F* h564374 = new TH1F();
  TH1F* h564356 = new TH1F();
  TH1F* h564359 = new TH1F();
  TH1F* h564373 = new TH1F();
  h564374->SetLineColor(cm.getColor(0));
  h564356->SetLineColor(cm.getColor(1));
  h564359->SetLineColor(cm.getColor(2));
  h564373->SetLineColor(cm.getColor(3));
  h564374->SetLineWidth(2);
  h564356->SetLineWidth(2);
  h564359->SetLineWidth(2);
  h564373->SetLineWidth(2);
  TLegend* legRun = new TLegend(0.77, 0.6, 1.1, 0.94);
  legRun->SetTextFont(42);
  legRun->SetBorderSize(0);
  legRun->SetFillColorAlpha(0, 0);
  legRun->AddEntry(h564374, "564374", "l");
  legRun->AddEntry(h564356, "564356", "l");
  legRun->AddEntry(h564359, "564359", "l");
  legRun->AddEntry(h564373, "564373", "l");

  TH1F* hFullCircle = new TH1F();
  TH1F* hFullSquare = new TH1F();
  hFullCircle->SetMarkerStyle(kFullCircle);
  hFullCircle->SetMarkerSize(2);
  hFullCircle->SetMarkerColor(kBlack);
  hFullSquare->SetMarkerStyle(kFullSquare);
  hFullSquare->SetMarkerSize(2);
  hFullSquare->SetMarkerColor(kBlack);
  TLegend* legHist = new TLegend(0.66, 0.77, 0.87, 0.94);
  legHist->SetTextFont(42);
  legHist->SetBorderSize(0);
  legHist->SetFillColorAlpha(0, 0);
  legHist->AddEntry(hFullCircle, "Coll", "p");
  legHist->AddEntry(hFullSquare, "BCs", "p");

  for (auto& hist : registry) {
    hist->histos[kNpart]->GetXaxis()->SetLabelSize(0);
    hist->histos[kNcoll]->GetXaxis()->SetLabelSize(0);
    hist->histos[kNpart]->GetYaxis()->SetLabelSize(0.045);
    hist->histos[kNcoll]->GetYaxis()->SetLabelSize(0.045);
    hist->histos[kNpart]->GetYaxis()->SetTitleSize(0.05);
    hist->histos[kNcoll]->GetYaxis()->SetTitleSize(0.05);
    hist->histos[kNpart]->GetYaxis()->SetTitle("Oxygen#minusOxygen <N_{part}>");
    hist->histos[kNcoll]->GetYaxis()->SetTitle("Oxygen#minusOxygen <N_{coll}>");
    hist->histos[kNpart]->GetYaxis()->SetTitleOffset(0);
    hist->histos[kNcoll]->GetYaxis()->SetTitleOffset(0);
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
  ooBCs_564374->histos[kNpart]->Draw("pe");
  ooBCs_564356->histos[kNpart]->Draw("pe same");
  ooBCs_564359->histos[kNpart]->Draw("pe same");
  ooBCs_564373->histos[kNpart]->Draw("pe same");
  ooColls_564374->histos[kNpart]->Draw("pe same");
  ooColls_564356->histos[kNpart]->Draw("pe same");
  ooColls_564359->histos[kNpart]->Draw("pe same");
  ooColls_564373->histos[kNpart]->Draw("pe same");
  legRun->Draw();
  legHist->Draw();

  canvNpart->cd();
  TPad* padNpartBot = new TPad("padNpartBot", "", 0, 0, 1, 0.3);
  padNpartBot->SetTicks(1, 1);
  padNpartBot->SetTopMargin(0.02);
  padNpartBot->SetBottomMargin(0.35);
  padNpartBot->SetLeftMargin(0.1);
  padNpartBot->SetRightMargin(0.03);
  padNpartBot->Draw();
  padNpartBot->cd();

  std::vector<std::pair<Histos*, Histos*>> runPairsNpart = {
    {ooBCs_564374, ooColls_564374},
    {ooBCs_564356, ooColls_564356},
    {ooBCs_564359, ooColls_564359},
    {ooBCs_564373, ooColls_564373}
  };

  std::vector<TH1F*> hRatiosNpart;
  for (int i = 0; i < 4; ++i) {
    TH1F* hRatio = (TH1F*)runPairsNpart[i].first->histos[kNpart]->Clone(Form("hRatioNpart_%d", i));
    hRatio->Divide(runPairsNpart[i].second->histos[kNpart]);
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("BCs / Colls");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatiosNpart.push_back(hRatio);
  }
  hRatiosNpart[0]->Draw("pe");
  for (int i = 1; i < 4; ++i) hRatiosNpart[i]->Draw("pe same");

  TLine* lineNpart = new TLine(hRatiosNpart[0]->GetXaxis()->GetXmin(), 1,
                                hRatiosNpart[0]->GetXaxis()->GetXmax(), 1);
  lineNpart->SetLineStyle(2);
  lineNpart->SetLineColor(kGray + 1);
  lineNpart->Draw("same");
  canvNpart->SaveAs("hNpartOO.pdf");

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
  ooBCs_564374->histos[kNcoll]->Draw("pe");
  ooBCs_564356->histos[kNcoll]->Draw("pe same");
  ooBCs_564359->histos[kNcoll]->Draw("pe same");
  ooBCs_564373->histos[kNcoll]->Draw("pe same");
  ooColls_564374->histos[kNcoll]->Draw("pe same");
  ooColls_564356->histos[kNcoll]->Draw("pe same");
  ooColls_564359->histos[kNcoll]->Draw("pe same");
  ooColls_564373->histos[kNcoll]->Draw("pe same");
  legRun->Draw();
  legHist->Draw();

  canvNcoll->cd();
  TPad* padNcollBot = new TPad("padNcollBot", "", 0, 0, 1, 0.3);
  padNcollBot->SetTicks(1, 1);
  padNcollBot->SetTopMargin(0.02);
  padNcollBot->SetBottomMargin(0.35);
  padNcollBot->SetLeftMargin(0.1);
  padNcollBot->SetRightMargin(0.03);
  padNcollBot->Draw();
  padNcollBot->cd();

  std::vector<std::pair<Histos*, Histos*>> runPairsNcoll = {
    {ooBCs_564374, ooColls_564374},
    {ooBCs_564356, ooColls_564356},
    {ooBCs_564359, ooColls_564359},
    {ooBCs_564373, ooColls_564373}
  };

  std::vector<TH1F*> hRatiosNcoll;
  for (int i = 0; i < 4; ++i) {
    TH1F* hRatio = (TH1F*)runPairsNcoll[i].first->histos[kNcoll]->Clone(Form("hRatioNcoll_%d", i));
    hRatio->Divide(runPairsNcoll[i].second->histos[kNcoll]);
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("BCs / Colls");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatiosNcoll.push_back(hRatio);
  }
  hRatiosNcoll[0]->Draw("pe");
  for (int i = 1; i < 4; ++i) hRatiosNcoll[i]->Draw("pe same");

  TLine* lineNcoll = new TLine(hRatiosNcoll[0]->GetXaxis()->GetXmin(), 1,
                                hRatiosNcoll[0]->GetXaxis()->GetXmax(), 1);
  lineNcoll->SetLineStyle(2);
  lineNcoll->SetLineColor(kGray + 1);
  lineNcoll->Draw("same");
  canvNcoll->SaveAs("hNcollOO.pdf");
}