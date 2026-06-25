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

std::vector<TH1F*> doNpartNcoll(TString dataset, TString ar, TString basehistos, const char* name)
{
  Bool_t fixedMu = true;
  Bool_t fixedK = true;
  Int_t ancesterMode = 2;
  Bool_t printNpartNcoll = true;

  TString mu = (fixedMu) ? "fixedMu" : "freeMu";
  TString k = (fixedK) ? "fixedK" : "freeK";
  TString filePath = Form("../AnalysisResults/%s/%s_calibration_ancestorMode%d_%s_%s_hFT0C_BCs%s.root", dataset.Data(), ar.Data(), ancesterMode, k.Data(), mu.Data(), basehistos.Data());
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

  h2dNpart->GetXaxis()->SetTitle("Centrality percentile (FT0C)");
  h2dNcoll->GetXaxis()->SetTitle("Centrality percentile (FT0C)");
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
  h2dNpart->Draw("colz");
  hAvgNpart->Draw("pe same");
  // canvNpart->SaveAs(Form("h2dNpart_%s.pdf", name));

  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 800);
  canvNcoll->SetLogz();
  canvNcoll->SetTopMargin(0.03);
  canvNcoll->SetLeftMargin(0.1);
  canvNcoll->SetRightMargin(0.1);
  canvNcoll->SetBottomMargin(0.13);
  h2dNcoll->Draw("colz");
  hAvgNcoll->Draw("pe same");
  // canvNcoll->SaveAs(Form("h2dNcoll_%s.pdf", name));

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

  const char* runNumber;
  std::vector<TH1F*> histos;
};

void drawTGlauberMC()
{
  gStyle->SetOptStat(0);
  enum GlauberQuantity {
    kNpart = 0,
    kNcoll,
    kNanc,
    kEcc,
    kB
  };

  Histos* pbpb25_568117 = new Histos(doNpartNcoll("LHC25_pass1", "AR_568117", "", "PbPb25"), "568117");
  Histos* pbpb25_568117_Pb = new Histos(doNpartNcoll("LHC25_pass1", "AR_568117", "_basehistos_Pb", "PbPb25"), "568117");
  Histos* pbpb25_568117_PbHN = new Histos(doNpartNcoll("LHC25_pass1", "AR_568117", "_basehistos_PbHN", "PbPb25"), "568117");
  Histos* pbpb25_568117_Pbpn = new Histos(doNpartNcoll("LHC25_pass1", "AR_568117", "_basehistos_Pbpn", "PbPb25"), "568117");
  Histos* pbpb25_568117_Pbstar = new Histos(doNpartNcoll("LHC25_pass1", "AR_568117", "_basehistos_Pbstar", "PbPb25"), "568117");
  Histos* pbpb25_568117_Pbrw = new Histos(doNpartNcoll("LHC25_pass1", "AR_568117", "_basehistos_Pbrw", "PbPb25"), "568117");
  Histos* pbpb25_568117_Pbpnrw = new Histos(doNpartNcoll("LHC25_pass1", "AR_568117", "_basehistos_Pbpnrw", "PbPb25"), "568117");

  ColorManager cm(6);
  pbpb25_568117->setColor(cm.getColor(-1));
  pbpb25_568117_Pb->setColor(cm.getColor(0));
  pbpb25_568117_PbHN->setColor(cm.getColor(1));
  pbpb25_568117_Pbpn->setColor(cm.getColor(2));
  pbpb25_568117_Pbstar->setColor(cm.getColor(3));
  pbpb25_568117_Pbrw->setColor(cm.getColor(4));
  pbpb25_568117_Pbpnrw->setColor(cm.getColor(5));

  TLegend* legend = new TLegend(0.67, 0.7, 0.87, 0.95);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0, 0);
  legend->AddEntry(pbpb25_568117->histos[kNpart], "default", "pl");
  legend->AddEntry(pbpb25_568117_Pb->histos[kNpart], "Pb", "lp");
  legend->AddEntry(pbpb25_568117_PbHN->histos[kNpart], "PbHN", "pl");
  legend->AddEntry(pbpb25_568117_Pbpn->histos[kNpart], "Pbpn", "lp");
  legend->AddEntry(pbpb25_568117_Pbstar->histos[kNpart], "Pb*", "pl");
  legend->AddEntry(pbpb25_568117_Pbrw->histos[kNpart], "Pbrw", "pl");
  legend->AddEntry(pbpb25_568117_Pbpnrw->histos[kNpart], "Pbpnrw*", "pl");

  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 800);
  canvNpart->SetTopMargin(0.03);
  canvNpart->SetLeftMargin(0.1);
  canvNpart->SetRightMargin(0.1);
  canvNpart->SetBottomMargin(0.13);
  // hAvgNpart->Draw();
  // pbpb25_567905->histos[kNpart]->Draw("pe");
  pbpb25_568117->histos[kNpart]->Draw("pe");
  pbpb25_568117_Pb->histos[kNpart]->Draw("pe same");
  pbpb25_568117_PbHN->histos[kNpart]->Draw("pe same");
  pbpb25_568117_Pbpn->histos[kNpart]->Draw("pe same");
  pbpb25_568117_Pbstar->histos[kNpart]->Draw("pe same");
  pbpb25_568117_Pbrw->histos[kNpart]->Draw("pe same");
  pbpb25_568117_Pbpnrw->histos[kNpart]->Draw("pe same");
  legend->Draw();
  canvNpart->SaveAs("hNpartTGlauberMC.pdf");

  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 800);
  canvNcoll->SetTopMargin(0.03);
  canvNcoll->SetLeftMargin(0.1);
  canvNcoll->SetRightMargin(0.1);
  canvNcoll->SetBottomMargin(0.13);
  // hAvgNcoll->Draw();
  // pbpb25_567905->histos[kNcoll]->Draw("pe");
  pbpb25_568117->histos[kNcoll]->Draw("pe");
  pbpb25_568117_Pb->histos[kNcoll]->Draw("pe same");
  pbpb25_568117_PbHN->histos[kNcoll]->Draw("pe same");
  pbpb25_568117_Pbpn->histos[kNcoll]->Draw("pe same");
  pbpb25_568117_Pbstar->histos[kNcoll]->Draw("pe same");
  pbpb25_568117_Pbrw->histos[kNcoll]->Draw("pe same");
  pbpb25_568117_Pbpnrw->histos[kNcoll]->Draw("pe same");
  legend->Draw();
  canvNcoll->SaveAs("hNcollTGlauberMC.pdf");

  std::vector<Histos*> registry = { pbpb25_568117_Pb,
                                    pbpb25_568117_PbHN,
                                    pbpb25_568117_Pbpn,
                                    pbpb25_568117_Pbstar,
                                    pbpb25_568117_Pbrw,
                                    pbpb25_568117_Pbpnrw };
  
  TLegend* legend2 = new TLegend(0.17, 0.7, 0.37, 0.95);
  legend2->SetBorderSize(0);
  legend2->SetFillColorAlpha(0, 0);
  legend2->AddEntry(pbpb25_568117_Pb->histos[kNpart], "Pb", "lp");
  legend2->AddEntry(pbpb25_568117_PbHN->histos[kNpart], "PbHN", "pl");
  legend2->AddEntry(pbpb25_568117_Pbpn->histos[kNpart], "Pbpn", "lp");
  legend2->AddEntry(pbpb25_568117_Pbstar->histos[kNpart], "Pb*", "pl");
  legend2->AddEntry(pbpb25_568117_Pbrw->histos[kNpart], "Pbrw", "pl");
  legend2->AddEntry(pbpb25_568117_Pbpnrw->histos[kNpart], "Pbpnrw*", "pl");

  // Ratio histograms Npart
  TCanvas* canvNpartRatio = new TCanvas("canvNpartRatio", "", 1200, 800);
  canvNpartRatio->SetTopMargin(0.03);
  canvNpartRatio->SetLeftMargin(0.13);
  canvNpartRatio->SetRightMargin(0.03);
  canvNpartRatio->SetBottomMargin(0.13);

  TLine* line = new TLine(0, 1, 100, 1);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kGray);

  bool firstNpart = true;
  for (auto& hist : registry) {
    TH1F* hRatio = (TH1F*)hist->histos[kNpart]->Clone();
    hRatio->Divide(hist->histos[kNpart], pbpb25_568117->histos[kNpart], 1, 1, "");
    hRatio->GetYaxis()->SetRangeUser(0.985, 1.03);
    hRatio->GetYaxis()->SetTitle("<N_{part}> ratio to basehistos");
    hRatio->GetXaxis()->SetTitle("Centrality percentile (FT0C)");
    hRatio->Draw(firstNpart ? "hist" : "hist same");
    firstNpart = false;
  }

  legend2->Draw();
  line->Draw();
  canvNpartRatio->SaveAs("hNpartRatioTGlauberMC.pdf");

  // Ratio histograms Ncoll
  TCanvas* canvNcollRatio = new TCanvas("canvNcollRatio", "", 1200, 800);
  canvNcollRatio->SetTopMargin(0.03);
  canvNcollRatio->SetLeftMargin(0.13);
  canvNcollRatio->SetRightMargin(0.03);
  canvNcollRatio->SetBottomMargin(0.13);

  bool firstNcoll = true;
  for (auto& hist : registry) {
    TH1F* hRatio = (TH1F*)hist->histos[kNcoll]->Clone();
    hRatio->Divide(hist->histos[kNcoll], pbpb25_568117->histos[kNcoll], 1, 1, "B");

    hRatio->GetYaxis()->SetRangeUser(0.985, 1.09);
    hRatio->GetYaxis()->SetTitle("<N_{coll}> ratio to basehistos");
    hRatio->GetXaxis()->SetTitle("Centrality percentile (FT0C)");
    hRatio->Draw(firstNcoll ? "hist" : "hist same");
    firstNcoll = false;
  }
  legend2->Draw();
  line->Draw();
  canvNcollRatio->SaveAs("hNcollRatioTGlauberMC.pdf");                
}