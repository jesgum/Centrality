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

  TString filePath = Form("../AnalysisResults/%s/%s_calibration_%s.root", dataset.Data(), ar.Data(), name);
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

  Histos* ooBCs = new Histos(doNpartNcoll("LHC25ae_pass2_small", "AR_564374", "hFT0M_BCs"), "BCs");
  Histos* ooColls = new Histos(doNpartNcoll("LHC25ae_pass2_small", "AR_564374", "hFT0M_Collisions"), "Collisions");

  std::vector<Histos*> registry = { ooBCs,
                                    ooColls };
                                    

  ooBCs->setColor(kBlue + 1);
  ooColls->setColor(kRed + 1);

  TLegend* legend = new TLegend(0.67, 0.64, 0.97, 0.94);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0, 0);

  for (auto& hist : registry) {
    legend->AddEntry(hist->histos[kNpart], hist->runNumber, "pl");
    hist->histos[kNpart]->GetXaxis()->SetLabelSize(0);
    hist->histos[kNcoll]->GetXaxis()->SetLabelSize(0);
    hist->histos[kNpart]->GetYaxis()->SetLabelSize(0.045);
    hist->histos[kNcoll]->GetYaxis()->SetLabelSize(0.045);
    hist->histos[kNpart]->GetYaxis()->SetTitleSize(0.05);
    hist->histos[kNcoll]->GetYaxis()->SetTitleSize(0.05);
    hist->histos[kNpart]->GetYaxis()->SetTitle("Oxygen#minusOxygen <N_{part}>");
    hist->histos[kNcoll]->GetYaxis()->SetTitle("Oxygen#minusOxygen <N_{coll}>");
  }

// === Canvas: Npart with ratio ===
TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 800);
canvNpart->SetTopMargin(0.0);
canvNpart->SetBottomMargin(0.0);

// Upper pad (main plot)
TPad* padNpartTop = new TPad("padNpartTop", "", 0, 0.3, 1, 1);
padNpartTop->SetTicks(1, 1);
padNpartTop->SetTopMargin(0.04);
padNpartTop->SetBottomMargin(0.02);
padNpartTop->SetLeftMargin(0.1);
padNpartTop->SetRightMargin(0.03);
padNpartTop->Draw();
padNpartTop->cd();
ooBCs->histos[kNpart]->Draw("pe");
ooColls->histos[kNpart]->Draw("pe same");
legend->Draw();

// Lower pad (ratio)
canvNpart->cd();
TPad* padNpartBot = new TPad("padNpartBot", "", 0, 0, 1, 0.3);
padNpartBot->SetTicks(1, 1);
padNpartBot->SetTopMargin(0.02);
padNpartBot->SetBottomMargin(0.35);
padNpartBot->SetLeftMargin(0.1);
padNpartBot->SetRightMargin(0.03);
padNpartBot->Draw();
padNpartBot->cd();

TH1* hRatioNpart = (TH1*)ooBCs->histos[kNpart]->Clone("hRatioNpart");
hRatioNpart->Divide(ooColls->histos[kNpart]);
hRatioNpart->SetMarkerColor(kBlack);
hRatioNpart->SetLineColor(kBlack);
hRatioNpart->SetTitle("");
hRatioNpart->GetYaxis()->SetTitle("BCs / Colls");
hRatioNpart->GetYaxis()->SetNdivisions(505);
hRatioNpart->GetYaxis()->SetTitleSize(0.12);
hRatioNpart->GetYaxis()->SetTitleOffset(0.4);
hRatioNpart->GetYaxis()->SetLabelSize(0.10);
hRatioNpart->GetXaxis()->SetTitleSize(0.13);
hRatioNpart->GetXaxis()->SetLabelSize(0.11);

TLine* lineNpart = new TLine(hRatioNpart->GetXaxis()->GetXmin(), 1,
                              hRatioNpart->GetXaxis()->GetXmax(), 1);
lineNpart->SetLineStyle(2);
lineNpart->SetLineColor(kGray + 1);
hRatioNpart->Draw("pe");
lineNpart->Draw("same");
canvNpart->SaveAs("hNpartOO.pdf");

// === Canvas: Ncoll with ratio ===
TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 800);
canvNcoll->SetTopMargin(0.0);
canvNcoll->SetBottomMargin(0.0);

// Upper pad (main plot)
TPad* padNcollTop = new TPad("padNcollTop", "", 0, 0.3, 1, 1);
padNcollTop->SetTicks(1, 1);
padNcollTop->SetTopMargin(0.04);
padNcollTop->SetBottomMargin(0.02);
padNcollTop->SetLeftMargin(0.1);
padNcollTop->SetRightMargin(0.03);
padNcollTop->Draw();
padNcollTop->cd();
ooBCs->histos[kNcoll]->Draw("pe");
ooColls->histos[kNcoll]->Draw("pe same");
legend->Draw();

// Lower pad (ratio)
canvNcoll->cd();
TPad* padNcollBot = new TPad("padNcollBot", "", 0, 0, 1, 0.3);
padNcollBot->SetTicks(1, 1);
padNcollBot->SetTopMargin(0.02);
padNcollBot->SetBottomMargin(0.35);
padNcollBot->SetLeftMargin(0.1);
padNcollBot->SetRightMargin(0.03);
padNcollBot->Draw();
padNcollBot->cd();

TH1* hRatioNcoll = (TH1*)ooBCs->histos[kNcoll]->Clone("hRatioNcoll");
hRatioNcoll->Divide(ooColls->histos[kNcoll]);
hRatioNcoll->SetMarkerColor(kBlack);
hRatioNcoll->SetLineColor(kBlack);
hRatioNcoll->SetTitle("");
hRatioNcoll->GetYaxis()->SetTitle("BCs / Colls");
hRatioNcoll->GetYaxis()->SetNdivisions(505);
hRatioNcoll->GetYaxis()->SetTitleSize(0.12);
hRatioNcoll->GetYaxis()->SetTitleOffset(0.4);
hRatioNcoll->GetYaxis()->SetLabelSize(0.10);
hRatioNcoll->GetXaxis()->SetTitleSize(0.13);
hRatioNcoll->GetXaxis()->SetLabelSize(0.11);

TLine* lineNcoll = new TLine(hRatioNcoll->GetXaxis()->GetXmin(), 1,
                              hRatioNcoll->GetXaxis()->GetXmax(), 1);
lineNcoll->SetLineStyle(2);
lineNcoll->SetLineColor(kGray + 1);
hRatioNcoll->Draw("pe");
lineNcoll->Draw("same");
canvNcoll->SaveAs("hNcollOO.pdf");
}