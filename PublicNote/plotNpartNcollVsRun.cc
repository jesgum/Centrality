#include <cstdio>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"

#include "plotHelpers.h"

// <Npart> and <Ncoll> vs centrality, one curve per run, for each Pb-Pb dataset.
// Input: Files/<var>VsRun_<dataset>.root, canvas canv<var>_<dataset> with sub-pads
//   padLeg  legend with one entry per run
//   padTop  reference (LHC23_pass5, run 544122, is the frame) + one histogram per run
//   padBot  hRatioFrame ("Var / Base"), a line at 1, one ratio per run, dataset label

struct Dataset {
  const char* name; // as in the file name
  int year;
};

void plotNpartNcollVsRun()
{
  plot::init();

  const Dataset datasets[] = {
    { "LHC23_pass5", 2023 },
    { "LHC24_pass3", 2024 },
    { "LHC25_pass1", 2025 },
    { "LHC26_pass1", 2026 },
  };

  double titleSize = -1; // fractions of the canvas width, set from the first canvas
  double labelSize = -1;

  for (const Dataset& ds : datasets) {
    for (const auto& var : { "Npart", "Ncoll" }) {
      TFile* f = TFile::Open(Form("Files/%sVsRun_%s.root", var, ds.name));
      if (!f || f->IsZombie()) {
        printf("cannot open Files/%sVsRun_%s.root, skipping\n", var, ds.name);
        continue;
      }
      TCanvas* c = (TCanvas*)f->Get(Form("canv%s_%s", var, ds.name));
      c->SetName(Form("canv%sVsRun_%d", var, ds.year));
      c->Draw();
      plot::removeStats(c);

      printf("\n=== %s %d ===\n", var, ds.year);
      plot::printCanvasContents(c);

      const bool isNpart = TString(var) == "Npart";
      const bool isNcoll = TString(var) == "Ncoll";

      // Upper pad: frame is the reference histogram, followed by one histogram per run
      TPad* padTop = plot::get<TPad>(c, "padTop");
      std::vector<TH1*> hists = plot::getHistograms(padTop);
      TH1* hFrame = hists.at(0);
      if (isNpart) {
        hFrame->GetYaxis()->SetTitle("#LT#it{N}_{part}#GT");
      }
      if (isNcoll) {
        hFrame->GetYaxis()->SetTitle("#LT#it{N}_{coll}#GT");
      }
      // hFrame->GetXaxis()->SetRangeUser(0, 90);
      // hFrame->GetYaxis()->SetRangeUser(0, 420);
      padTop->Modified();

      // Lower pad: ratio to the reference
      TPad* padBot = plot::get<TPad>(c, "padBot");
      TH1* hRatio = plot::getFrame(c, "padBot"); // hRatioFrame
      hRatio->GetYaxis()->SetTitle("Run / reference");
      hRatio->GetXaxis()->SetTitle("FT0C Centrality (%)");
      padBot->Modified();

      // Legend with the run numbers (its own pad)
      TPad* padLeg = plot::get<TPad>(c, "padLeg");
      TLegend* leg = plot::get<TLegend>(padLeg, "TPave");
      // leg->SetNColumns(4);
      // leg->SetTextSize(0.03);
      padLeg->Modified();

      // Replace the dataset label by the year and add the ALICE label
      plot::removeLatex(c, ds.name);
      // Same absolute text size in all canvases: sizes as a fraction of the canvas
      // width, taken from the first (2023) canvas so it keeps its current look
      if (titleSize < 0) {
        titleSize = plot::toCanvasWidthFraction(padTop, 0.08);
        labelSize = plot::toCanvasWidthFraction(padTop, 0.07);
      }
      padTop->cd();
      TLatex* lat = new TLatex();
      lat->SetTextAlign(31);
      plot::setTextSizeOfCanvasWidth(lat, c, titleSize, 6);
      lat->DrawLatexNDC(0.92, 0.85, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");
      plot::setTextSizeOfCanvasWidth(lat, c, labelSize, 4);
      lat->DrawLatexNDC(0.92, 0.75, Form("Data-taking year: %d", ds.year));

      // Year-specific tweaks, e.g.
      // if (ds.year == 2026) { ... }

      c->Modified();
      c->Update();
      c->SaveAs(Form("Figures/%sVsRun_PbPb_%d.pdf", var, ds.year));
    }
  }
}
