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

// O-O <Npart> and <Ncoll> vs centrality, one curve per run (LHC25ae_pass2).
// Input: Files/h<var>OOVsRun.root (from Macros/drawLightIonNpartNcollVsRun.cc),
// canvas "canv" (1200x800) with sub-pads
//   padTop  reference run 564374 (frame) + one histogram per run + legend
//   padBot  ratio of each run to the reference (frame is the first ratio)

void plotNpartNcollOOVsRun()
{
  plot::init();

  for (const auto& var : { "Npart", "Ncoll" }) {
    const TString fileName = Form("Files/h%sOOVsRun.root", var);
    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
      printf("cannot open %s, skipping\n", fileName.Data());
      continue;
    }
    TCanvas* c = (TCanvas*)f->Get("canv");
    c->SetName(Form("canv%sVsRun_OO", var));
    c->Draw();
    plot::removeStats(c);

    // Canvas size in pixels (original 1200x800)
    // c->SetCanvasSize(1200, 800);

    printf("\n=== %s OO ===\n", var);
    plot::printCanvasContents(c);

    const bool isNpart = TString(var) == "Npart";

    TPad* padTop = plot::get<TPad>(c, "padTop");
    TPad* padBot = plot::get<TPad>(c, "padBot");

    // Margins of the sub-pads (left, right, bottom, top); originals from the file.
    // Keep left/right equal in both pads so the x-axes stay aligned.
    const double leftMargin = 0.1;
    const double rightMargin = 0.03;
    plot::setMargins(padTop, leftMargin, rightMargin, 0.02, 0.04);
    plot::setMargins(padBot, leftMargin, rightMargin, 0.35, 0.02);

    // Upper pad: frame is the reference run, followed by one histogram per run
    TH1* hFrame = plot::getFrame(padTop);
    hFrame->GetYaxis()->SetTitle(isNpart ? "#LT#it{N}_{part}#GT" : "#LT#it{N}_{coll}#GT");
    hFrame->SetMaximum(hFrame->GetMaximum() * 1.12);
    // hFrame->GetYaxis()->SetRangeUser(0, 30);
    padTop->Modified();

    // Lower pad: ratio to the reference run
    TH1* hRatio = plot::getFrame(padBot);
    hRatio->GetXaxis()->SetTitle("FT0M Centrality (%)");
    hRatio->GetYaxis()->SetTitle("Run / reference");
    // hRatio->GetYaxis()->SetRangeUser(0.95, 1.05);
    plot::drawHorizontalLine(padBot, 1.); // line at unity, behind the ratios
    padBot->Modified();

    // Legend with the run numbers: position in NDC of padTop (x1, y1, x2, y2),
    // originals from the file
    TLegend* leg = plot::get<TLegend>(padTop, "TPave");
    leg->SetX1NDC(0.75);
    leg->SetY1NDC(0.35);
    leg->SetX2NDC(0.95);
    leg->SetY2NDC(0.9);
    // leg->SetTextSize(0.04);
    padTop->Modified();

    // ALICE label, left of the legend above the curves
    padTop->cd();
    TLatex* lat = new TLatex();
    const float latPosX = 0.22;
    lat->SetTextSize(0.09);
    lat->DrawLatexNDC(latPosX, 0.84, "ALICE OO #sqrt{#it{s}_{NN}} = 5.36 TeV");

    c->Modified();
    c->Update();
    c->SaveAs(Form("Figures/%sVsRun_OO.pdf", var));
  }
}
