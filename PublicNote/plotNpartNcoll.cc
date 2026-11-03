#include <cstdio>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLatex.h"
#include "TString.h"

#include "plotHelpers.h"

struct SystemInfo {
  const char* name;       // used in the output file name
  const char* fileTag;    // Files/h2d<var>_<fileTag>.root
  const char* canvSuffix; // canvas in the file is canv<var><canvSuffix>
};

void plotNpartNcoll()
{
  plot::init();

  const SystemInfo systems[] = {
    { "OO", "hFT0M_BCsOO", "" },
    { "PbPb", "hFT0C_PbPb", "PbPb" },
  };

  for (const SystemInfo& sys : systems) {
    for (const char* var : { "Npart", "Ncoll" }) {
      TFile* f = TFile::Open(Form("Files/h2d%s_%s.root", var, sys.fileTag));
      TCanvas* c = (TCanvas*)f->Get(Form("canv%s%s", var, sys.canvSuffix));
      c->SetName(Form("canv%s_%s", var, sys.name));
      c->Draw();
      plot::removeStats(c);

      // Canvas size in pixels (originals: OO 1200x800, PbPb 800x800)
      c->SetCanvasSize(800, 800);
      plot::setMargins(c, 0.15, 0.15, 0.15, 0.15);  // left, right, bottom, top

      printf("\n=== %s %s ===\n", sys.name, var);
      plot::printCanvasContents(c);

      // The histograms in the plot, in drawing order
      std::vector<TH1*> hists = plot::getHistograms(c);
      TH2* h2d = (TH2*)hists.at(0);  // the 2D Npart/Ncoll vs centrality (also the frame)
      h2d->GetXaxis()->SetTitle("FT0C Centrality (%)");
      TH1* hAvg = hists.at(1);       // average markers per centrality bin

      // Which plot is this?
      const bool isOO = TString(sys.name) == "OO";
      const bool isPbPb = TString(sys.name) == "PbPb";
      const bool isNpart = TString(var) == "Npart";
      const bool isNcoll = TString(var) == "Ncoll";

      if (isNpart) {
        h2d->GetYaxis()->SetTitle("#it{N}_{part}");
      }
      if (isNcoll) {
        h2d->GetYaxis()->SetTitle("#it{N}_{coll}");
      }
      
      
      TLatex* lat = new TLatex();
      lat->SetTextSize(0.04);
      lat->SetTextAlign(31);
      const float latPosX = 0.84;
      if (isOO && isNpart) {
        lat->DrawLatexNDC(latPosX, 0.8, "ALICE OO #sqrt{#it{s}_{NN}} = 5.36 TeV");
        h2d->GetYaxis()->SetRangeUser(0, 37);
        plot::setMinimum(c, 1e-2);
        plot::setMaximum(c, 2e+4);
      }
      if (isOO && isNcoll) {
        lat->DrawLatexNDC(latPosX, 0.8, "ALICE OO #sqrt{#it{s}_{NN}} = 5.36 TeV");
        h2d->RebinY(3);
        h2d->GetYaxis()->SetRangeUser(0, 112);
        plot::setMinimum(c, 1e-5);
        plot::setMaximum(c, 2e+4);
      }
      if (isPbPb && isNpart) {
        lat->DrawLatexNDC(latPosX, 0.8, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");
        // h2d->GetYaxis()->SetRangeUser(0, 37);
        plot::setMinimum(c, 1e-2);
        plot::setMaximum(c, 2e+4);
      }
      if (isPbPb && isNcoll) {
        lat->DrawLatexNDC(latPosX, 0.8, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");
        h2d->GetYaxis()->SetRangeUser(0, 2400);
        plot::setMinimum(c, 1e-2);
        plot::setMaximum(c, 2e+4);
      }

      c->SaveAs(Form("Figures/h2d%s_%s.pdf", var, sys.name));
    }
  }
}
