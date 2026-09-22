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

#include "Includes/colorManager.h"
#include "Includes/centralityHelper.h"
#include "Includes/systematicsHelper.h"

void drawNpartNcollTrento()
{
  gStyle->SetOptStat(0);
  syst::GlauberParameters base("LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.root");
  syst::GlauberParameters<TH2D> trento("../Trento/trentofit_calib.root", true);
  
  TCanvas* canvNpart = new TCanvas("canvNpart", "", 1200, 800);
  canvNpart->SetLogz();
  canvNpart->SetTopMargin(0.03);
  canvNpart->SetLeftMargin(0.1);
  canvNpart->SetRightMargin(0.1);
  canvNpart->SetBottomMargin(0.13);
  trento.h2dNpart->Draw("colz");
  trento.hNpart->Draw("same");
  canvNpart->SaveAs("hNpart.pdf");

  TCanvas* canvNcoll = new TCanvas("canvNcoll", "", 1200, 800);
  canvNcoll->SetLogz();
  canvNcoll->SetTopMargin(0.03);
  canvNcoll->SetLeftMargin(0.1);
  canvNcoll->SetRightMargin(0.1);
  canvNcoll->SetBottomMargin(0.13);
  trento.h2dNcoll->Draw("colz");
  trento.hNcoll->Draw("same");
  canvNcoll->SaveAs("hNcoll.pdf");

  trento.hNpart->SetMarkerColor(kRed + 1);
  trento.hNpart->SetLineColor(kRed + 1);
  trento.hNcoll->SetMarkerColor(kRed + 1);
  trento.hNcoll->SetLineColor(kRed + 1);
  base.hNpart->SetMarkerColor(kBlue + 1);
  base.hNpart->SetLineColor(kBlue + 1);
  base.hNcoll->SetMarkerColor(kBlue + 1);
  base.hNcoll->SetLineColor(kBlue + 1);

  // === Canvas: Npart, Trento vs base with ratio ===
  TCanvas* canvNpartCompare = new TCanvas("canvNpartCompare", "", 1200, 800);
  canvNpartCompare->SetTopMargin(0.0);
  canvNpartCompare->SetBottomMargin(0.0);

  TPad* padNpartTop = new TPad("padNpartTop", "", 0, 0.3, 1, 1);
  padNpartTop->SetTicks(1, 1);
  padNpartTop->SetTopMargin(0.04);
  padNpartTop->SetBottomMargin(0.02);
  padNpartTop->SetLeftMargin(0.12);
  padNpartTop->SetRightMargin(0.03);
  padNpartTop->Draw();
  padNpartTop->cd();
  base.hNpart->SetMinimum(0);
  base.hNpart->GetXaxis()->SetLabelSize(0);
  base.hNpart->GetYaxis()->SetTitle("<N_{part}>");
  base.hNpart->GetYaxis()->SetTitleSize(0.05);
  base.hNpart->GetYaxis()->SetLabelSize(0.045);
  base.hNpart->GetYaxis()->SetTitleOffset(1.1);
  base.hNpart->Draw("pe");
  trento.hNpart->Draw("pe same");

  TLegend* legNpartCompare = new TLegend(0.65, 0.78, 0.87, 0.94);
  legNpartCompare->SetBorderSize(0);
  legNpartCompare->SetFillColorAlpha(0, 0);
  legNpartCompare->AddEntry(base.hNpart, "Base", "p");
  legNpartCompare->AddEntry(trento.hNpart, "Trento", "p");
  legNpartCompare->Draw();

  canvNpartCompare->cd();
  TPad* padNpartBot = new TPad("padNpartBot", "", 0, 0, 1, 0.3);
  padNpartBot->SetTicks(1, 1);
  padNpartBot->SetTopMargin(0.02);
  padNpartBot->SetBottomMargin(0.35);
  padNpartBot->SetLeftMargin(0.12);
  padNpartBot->SetRightMargin(0.03);
  padNpartBot->Draw();
  padNpartBot->cd();

  TH1F* hRatioNpart = dynamic_cast<TH1F*>(trento.hNpart->Clone("hRatioNpart"));
  hRatioNpart->Divide(base.hNpart);
  hRatioNpart->SetTitle("");
  hRatioNpart->GetYaxis()->SetTitle("Trento / Base");
  hRatioNpart->GetYaxis()->SetNdivisions(505);
  hRatioNpart->GetYaxis()->SetTitleSize(0.12);
  hRatioNpart->GetYaxis()->SetTitleOffset(0.45);
  hRatioNpart->GetYaxis()->SetLabelSize(0.10);
  hRatioNpart->GetXaxis()->SetTitle("Centrality (%)");
  hRatioNpart->GetXaxis()->SetTitleSize(0.13);
  hRatioNpart->GetXaxis()->SetLabelSize(0.11);
  hRatioNpart->Draw("pe");

  TLine* lineNpartCompare = new TLine(hRatioNpart->GetXaxis()->GetXmin(), 1,
                                       hRatioNpart->GetXaxis()->GetXmax(), 1);
  lineNpartCompare->SetLineStyle(2);
  lineNpartCompare->SetLineColor(kGray + 1);
  lineNpartCompare->Draw("same");

  canvNpartCompare->SaveAs("hNpartCompareTrentoBase.pdf");

  // === Canvas: Ncoll, Trento vs base with ratio ===
  TCanvas* canvNcollCompare = new TCanvas("canvNcollCompare", "", 1200, 800);
  canvNcollCompare->SetTopMargin(0.0);
  canvNcollCompare->SetBottomMargin(0.0);

  TPad* padNcollTop = new TPad("padNcollTop", "", 0, 0.3, 1, 1);
  padNcollTop->SetTicks(1, 1);
  padNcollTop->SetTopMargin(0.04);
  padNcollTop->SetBottomMargin(0.02);
  padNcollTop->SetLeftMargin(0.12);
  padNcollTop->SetRightMargin(0.03);
  padNcollTop->Draw();
  padNcollTop->cd();
  base.hNcoll->SetMinimum(0);
  base.hNcoll->GetXaxis()->SetLabelSize(0);
  base.hNcoll->GetYaxis()->SetTitle("<N_{coll}>");
  base.hNcoll->GetYaxis()->SetTitleSize(0.05);
  base.hNcoll->GetYaxis()->SetLabelSize(0.045);
  base.hNcoll->GetYaxis()->SetTitleOffset(1.1);
  base.hNcoll->Draw("pe");
  trento.hNcoll->Draw("pe same");

  TLegend* legNcollCompare = new TLegend(0.65, 0.78, 0.87, 0.94);
  legNcollCompare->SetBorderSize(0);
  legNcollCompare->SetFillColorAlpha(0, 0);
  legNcollCompare->AddEntry(base.hNcoll, "Base", "p");
  legNcollCompare->AddEntry(trento.hNcoll, "Trento", "p");
  legNcollCompare->Draw();

  canvNcollCompare->cd();
  TPad* padNcollBot = new TPad("padNcollBot", "", 0, 0, 1, 0.3);
  padNcollBot->SetTicks(1, 1);
  padNcollBot->SetTopMargin(0.02);
  padNcollBot->SetBottomMargin(0.35);
  padNcollBot->SetLeftMargin(0.12);
  padNcollBot->SetRightMargin(0.03);
  padNcollBot->Draw();
  padNcollBot->cd();

  TH1F* hRatioNcoll = dynamic_cast<TH1F*>(trento.hNcoll->Clone("hRatioNcoll"));
  hRatioNcoll->Divide(base.hNcoll);
  hRatioNcoll->SetTitle("");
  hRatioNcoll->GetYaxis()->SetTitle("Trento / Base");
  hRatioNcoll->GetYaxis()->SetNdivisions(505);
  hRatioNcoll->GetYaxis()->SetTitleSize(0.12);
  hRatioNcoll->GetYaxis()->SetTitleOffset(0.45);
  hRatioNcoll->GetYaxis()->SetLabelSize(0.10);
  hRatioNcoll->GetXaxis()->SetTitle("Centrality (%)");
  hRatioNcoll->GetXaxis()->SetTitleSize(0.13);
  hRatioNcoll->GetXaxis()->SetLabelSize(0.11);
  hRatioNcoll->Draw("pe");

  TLine* lineNcollCompare = new TLine(hRatioNcoll->GetXaxis()->GetXmin(), 1,
                                       hRatioNcoll->GetXaxis()->GetXmax(), 1);
  lineNcollCompare->SetLineStyle(2);
  lineNcollCompare->SetLineColor(kGray + 1);
  lineNcollCompare->Draw("same");

  canvNcollCompare->SaveAs("hNcollCompareTrentoBase.pdf");
}